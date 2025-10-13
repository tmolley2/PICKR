#08_score_pairs.py
import math, glob, os, re
import pandas as pd
from scipy.stats import norm
from pathlib import Path

# ------------------------------------------------------------------
#   Configuration
# ------------------------------------------------------------------
SPECIES_LIST = ["frog", "cow", "human", "monkey", "rat", "pig", "mouse", "zebrafish", "worm", "fly"]

# Base Directories
BASE_TEMP_DIR = Path('temp')
BASE_OUTPUT_DIR = Path('scored') # Final output will go here

# ------------------------------------------------------------------
#   Scoring hyper-parameters
# ------------------------------------------------------------------
N_MAX, C_MAX  = 100, 5000       # ceilings for evidence scaling
WE, WB, WS    = 0.5, 0.30, 0.20   # additive weights (Evidence, Bio, Synergy)
DIMER_HARD_CUTOFF = 6           # pair_comp ≥ this ⇒ PICKR = 0
P_MAX, R_MAX, T_MAX   = 25, 100, 500    # ceilings for log‐scaling
W_SHARED = 0.5
W_REG = 0.475
W_TOT = 0.425        # each tier = 0.5

# Parameters for percentile calculation & filtering
PICKR_MEAN = 54.6
PICKR_STDEV = 16.11
TOP_N_PAIRS = 30

# Derived constants & setup
SCORE_COLS = ["evidence_score", "biophysics_score", "synergy_score", "pickr_score"]
SAFE_FLOOR = 1e-3             # avoids log10(0)

# ------------------------------------------------------------------
#   Utility functions
# ------------------------------------------------------------------

def safe_log10(x: float) -> float:
    """log10 with small floor to prevent −inf."""
    return math.log10(max(x, SAFE_FLOOR))

# ✅ MODIFIED: Now correctly handles Xenopus .L/.S suffixes
def remove_pseudogenes(row: pd.Series) -> pd.Series:
    """
    Cleans every *_cross_specificity_[0-4] and its *_seq column in-place.
    Handles Xenopus .L/.S suffixes correctly.
    """
    # Get the base gene name, stripping .L or .S if present
    tgt_full = str(row.get("f_gene", "")).upper()
    tgt_base = tgt_full.rsplit('.', 1)[0] if tgt_full.endswith(('.L', '.S')) else tgt_full

    for side in ("f", "r"):
        for d in range(0, 5):
            g_col, s_col = f"{side}_cross_specificity_{d}", f"{side}_cross_specificity_{d}_seq"
            if g_col not in row or s_col not in row: continue

            genes_raw = str(row[g_col]) if pd.notna(row[g_col]) else ""
            seqs_raw  = str(row[s_col]) if pd.notna(row[s_col]) else ""
            genes = [g.strip() for g in genes_raw.split(",") if g.strip()]
            seqs  = [s.strip() for s in seqs_raw.split(";") if s.strip()]

            keep_genes, keep_seqs = [], []
            for g, aln in zip(genes, seqs):
                g_up = g.upper()
                if g_up.startswith(("LOC", "GM")) or g_up.endswith('RIK'): continue
                # Don't filter out the other homolog (e.g., keep .S when target is .L)
                if tgt_base in g_up and g_up != tgt_full and g_up.rsplit('.', 1)[0] != tgt_base:
                    continue
                keep_genes.append(g)
                keep_seqs.append(aln.split(":", 1)[-1])

            row[g_col] = ",".join(keep_genes)
            row[s_col] = ";".join(keep_seqs)
    return row

def cross_specificity_score(aln: str, primer_len: int) -> float:
    """Calculates the specificity factor for a single off-target alignment."""
    if ":" in aln: aln = aln.split(":", 1)[1]
    mismatches = [i+1 for i, c in enumerate(aln) if c != '.']
    d = len(mismatches)
    if d == 0: return 1.0

    BASE_FACTOR = {1: 0.10, 2: 0.25, 3: 0.65, 4: 0.9}
    NEAR_3_CUTOFF, POS_WEIGHT = 3, 0.60
    if d > 4: return 1.0

    dist_from_3 = [primer_len - i for i in mismatches]
    n_near = sum(1 for d3 in dist_from_3 if d3 <= NEAR_3_CUTOFF)
    frac_near = n_near / d
    rescue = (1.0 - BASE_FACTOR[d]) * POS_WEIGHT * frac_near
    return max(0.0, min(BASE_FACTOR[d] + rescue, 1.0))

def primer_quality(length_nt: int, gc_pct: float, tm_c: float) -> float:
    """Average compliance score ∈[0,1] across length, GC%, Tm."""
    if 19 <= length_nt <= 25: length_s = 1.0
    else: length_s = max(0.0, 1 - 0.2 * abs(length_nt - 22))
    if 45 <= gc_pct <= 55: gc_s = 1.0
    elif 40 <= gc_pct < 45 or 55 < gc_pct <= 60: gc_s = 0.75
    else: gc_s = max(0.0, 0.75 - (abs(gc_pct - 50) - 10) * 0.075)
    tm_s = max(0.0, 1 - abs(tm_c - 60) / 3)
    return (length_s + gc_s + tm_s) / 3

def evidence_score(row):
    """Calculates the evidence component of the PICKR score."""
    shared_norm = min(safe_log10(row.pair_shared_pmcid_count + 1) / math.log10(P_MAX + 1), 1.0)
    reg_f_norm = min(safe_log10(row.f_reg_pmcid_count + 1) / math.log10(R_MAX + 1), 1.0)
    reg_r_norm = min(safe_log10(row.r_reg_pmcid_count + 1) / math.log10(R_MAX + 1), 1.0)
    tot_f_norm = min(safe_log10(row.f_reg_pmcid_count + row.f_inv_pmcid_count + 1) / math.log10(T_MAX + 1), 1.0)
    tot_r_norm = min(safe_log10(row.r_reg_pmcid_count + row.r_inv_pmcid_count + 1) / math.log10(T_MAX + 1), 1.0)
    return min(shared_norm*W_SHARED + ((reg_f_norm+reg_r_norm)/2)*W_REG + ((tot_f_norm+tot_r_norm)/2)*W_TOT, 1.0)

def synergy_score(tm_f: float, tm_r: float, pair_comp: float) -> float:
    """Calculates the synergy component of the PICKR score."""
    t_penalty = abs(tm_f - tm_r) / 3
    dimer_pen = max(0.0, (pair_comp - 6) / 6)
    return max(0.0, 1 - t_penalty - dimer_pen)

def pickr_components(row):
    """Calculates all components of the PICKR score for a single primer pair."""
    def has_real_perfect_off(side: str) -> bool:
        raw = row.get(f"{side}_cross_specificity_0")
        if pd.isna(raw) or not raw: return False
        tgt_full = str(row.get('f_gene', '')).upper()
        tgt_base = tgt_full.rsplit('.', 1)[0] if tgt_full.endswith(('.L', '.S')) else tgt_full
        for g in re.split(r'[;|\n,]+', str(raw)):
            g_clean = g.strip().upper()
            if not g_clean or g_clean.startswith(('LOC', 'GM')) or g_clean.endswith('RIK'): continue
            if tgt_base in g_clean and g_clean != tgt_full and g_clean.rsplit('.', 1)[0] != tgt_base: continue
            return True
        return False

    if has_real_perfect_off('f') or has_real_perfect_off('r') or row.pair_comp >= DIMER_HARD_CUTOFF:
        return {col: 0.0 for col in SCORE_COLS}

    e = evidence_score(row)
    qual_f = primer_quality(len(str(row.f_sequence)), row.f_GC_pct, row.f_Tm_C)
    qual_r = primer_quality(len(str(row.r_sequence)), row.r_GC_pct, row.r_Tm_C)
    b_uncorr = 0.5 * (qual_f + qual_r)

    def compound_specificity(side: str) -> float:
        score = 1.0
        tgt_full = str(row.get('f_gene', '')).upper()
        tgt_base = tgt_full.rsplit('.', 1)[0] if tgt_full.endswith(('.L', '.S')) else tgt_full
        for d in range(1, 5):
            genes_raw = row.get(f"{side}_cross_specificity_{d}")
            seqs_raw = row.get(f"{side}_cross_specificity_{d}_seq")
            if pd.isna(genes_raw) or pd.isna(seqs_raw): continue
            for g, aln in zip(str(genes_raw).split(','), str(seqs_raw).split(';')):
                g_clean = g.strip().upper()
                if not g_clean or g_clean.startswith(('LOC', 'GM')) or g_clean.endswith('RIK'): continue
                if tgt_base in g_clean and g_clean.rsplit('.', 1)[0] == tgt_base: continue
                score *= cross_specificity_score(aln, len(str(row[f"{side}_sequence"])))
        return score

    b = b_uncorr * compound_specificity('f') * compound_specificity('r')
    s = synergy_score(row.f_Tm_C, row.r_Tm_C, row.pair_comp)
    pickr = 100 * (WE * e + WB * b + WS * s)
    return dict(evidence_score=e, biophysics_score=b, synergy_score=s, pickr_score=pickr)

# ------------------------------------------------------------------
# Main loop
# ------------------------------------------------------------------
for species in SPECIES_LIST:
    print(f"\n{'='*20} SCORING PRIMERS FOR: {species.upper()} {'='*20}")

    INPUT_FOLDER = BASE_TEMP_DIR / f'primer_pairings_shared_counts_{species}'
    OUTPUT_FOLDER = BASE_OUTPUT_DIR / f'csv_data_{species}'

    if not INPUT_FOLDER.is_dir():
        print(f"Warning: Input directory '{INPUT_FOLDER}' not found. Skipping {species}.")
        continue
    
    OUTPUT_FOLDER.mkdir(exist_ok=True)

    file_list = glob.glob(os.path.join(INPUT_FOLDER, "*.csv"))

    for csv_path in file_list:
        gene = os.path.splitext(os.path.basename(csv_path))[0]
        print(f"Processing: {gene}")

        try:
            df = pd.read_csv(csv_path)
            # Ensure numeric columns are treated as numbers, filling errors with 0
            numeric_cols = ['f_GC_pct', 'f_Tm_C', 'r_GC_pct', 'r_Tm_C', 'pair_comp', 
                            'pair_shared_pmcid_count', 'f_reg_pmcid_count', 'f_inv_pmcid_count',
                            'r_reg_pmcid_count', 'r_inv_pmcid_count']
            for col in numeric_cols:
                df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0)
        except pd.errors.EmptyDataError:
            print(f"  - Skipping empty file: {gene}.csv")
            continue
        
        # 1) Filter out pseudogenes from cross-specificity results
        df = df.apply(remove_pseudogenes, axis=1)

        # 2) Score each row
        comp_df = df.apply(pickr_components, axis=1, result_type="expand")

        # 3) Attach the new score columns
        df.drop(columns=[c for c in SCORE_COLS if c in df.columns], inplace=True, errors="ignore")
        df_final = pd.concat([df, comp_df], axis=1)

        # 4) Calculate Percentile based on the pickr_score
        z_scores = (df_final['pickr_score'] - PICKR_MEAN) / PICKR_STDEV
        df_final['percentile'] = norm.cdf(z_scores) * 100

        # 5) Filter out any pairs that scored 0
        original_count = len(df_final)
        df_final = df_final[df_final['pickr_score'] > 0].copy()
        
        if original_count > len(df_final):
            print(f"  - Dropped {original_count - len(df_final)} pairs with a score of 0.")
        
        # 6) Sort by score and keep only the top N pairs
        df_final.sort_values(by='pickr_score', ascending=False, inplace=True)
        if len(df_final) > TOP_N_PAIRS:
            print(f"  - Keeping top {TOP_N_PAIRS} of {len(df_final)} pairs.")
            df_final = df_final.head(TOP_N_PAIRS)
        
        final_count = len(df_final)

        # 7) Save the final scored and filtered file
        if final_count > 0:
            outfile = os.path.join(OUTPUT_FOLDER, f"{gene}.csv")
            df_final.to_csv(outfile, index=False)
            print(f"✓ {gene}: scored and saved → {outfile} ({final_count} pairs)")
        else:
            print(f"  - No valid pairs remained for {gene} after scoring.")

print("\n--- All species have been scored. ---")