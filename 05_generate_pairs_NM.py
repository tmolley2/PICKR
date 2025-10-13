# 05_generate_pairs.py


import os
import pandas as pd
from pathlib import Path
from Bio.Seq import Seq

# ----------------------------------------------------------------------------------------------------------------------------------
# 1) Configuration
# ----------------------------------------------------------------------------------------------------------------------------------
# List of species to process
SPECIES_LIST = ["frog", "zebrafish", "worm", "fly"]

# Base Input/Output Directories
BASE_INPUT_DIR = Path("temp")
BASE_OUTPUT_DIR = Path("temp")

# ----------------------------------------------------------------------------------------------------------------------------------
# 2) Utility Functions
# ----------------------------------------------------------------------------------------------------------------------------------

def max_contiguous_matches(s1: str, s2: str) -> int:
    """
    Align s1 against s2 at all offsets and return the maximum contiguous matching run.
    """
    s1, s2 = s1.upper(), s2.upper()
    max_run = 0
    for shift in range(-len(s2) + 1, len(s1)):
        run = 0
        for i in range(len(s1)):
            j = i - shift
            if 0 <= j < len(s2) and s1[i] == s2[j]:
                run += 1
                max_run = max(max_run, run)
            else:
                run = 0
    return max_run

def self_complementarity(primer_seq: str) -> int:
    """Calculates self-complementarity: primer vs its own reverse-complement."""
    rc = str(Seq(primer_seq).reverse_complement())
    return max_contiguous_matches(primer_seq, rc)

def pair_complementarity(fwd_seq: str, rev_seq: str) -> int:
    """Calculates pair complementarity: forward primer vs reverse complement of reverse primer."""
    rev_rc = str(Seq(rev_seq).reverse_complement())
    return max_contiguous_matches(fwd_seq, rev_rc)

# ----------------------------------------------------------------------------------------------------------------------------------
# 3) Main Execution Block
# ----------------------------------------------------------------------------------------------------------------------------------

# ✅ NEW: Main loop to iterate over each species
for species in SPECIES_LIST:
    print(f"\n{'='*20} GENERATING PAIRS FOR SPECIES: {species.upper()} {'='*20}")

    # ✅ NEW: Dynamically set input and output directories for the current species
    ANNOTATED_DIR = BASE_INPUT_DIR / f'gene_csvs/{species}'
    PAIRINGS_DIR = BASE_OUTPUT_DIR / f'primer_pairings_{species}'

    # Ensure the species-specific directories exist
    if not ANNOTATED_DIR.is_dir():
        print(f"Warning: Input directory '{ANNOTATED_DIR}' not found. Skipping {species}.")
        continue
    
    PAIRINGS_DIR.mkdir(exist_ok=True)

    # Loop through each annotated CSV in the species-specific folder
    for fn in os.listdir(ANNOTATED_DIR):
        if not fn.endswith('.csv'):
            continue
        
        gene = fn.replace('.csv', '')
        input_path = ANNOTATED_DIR / fn
        
        try:
            df = pd.read_csv(input_path)
        except pd.errors.EmptyDataError:
            print(f"Skipping empty file: {fn}")
            continue

        # Separate forward vs reverse_complement primers
        df_fwd = df[df['orientation'] == 'forward'].copy()
        df_rev = df[df['orientation'] == 'reverse_complement'].copy()

        if df_fwd.empty or df_rev.empty:
            print(f"  - Skipping {gene}: not enough forward/reverse primers to form pairs.")
            continue

        pair_records = []
        # Generate all possible forward-reverse pairs
        for f in df_fwd.itertuples(index=False):
            for r in df_rev.itertuples(index=False):
                # Ensure start/end columns are numeric before calculation
                f_start = pd.to_numeric(f.start, errors='coerce')
                r_end = pd.to_numeric(r.end, errors='coerce')
                
                if pd.isna(f_start) or pd.isna(r_end):
                    continue

                gap = r_end - f_start
                
                # Check for valid amplicon size
                if 50 < gap < 300:
                    record = {}
                    # Add all columns with 'f_' and 'r_' prefixes
                    for col in df.columns:
                        record[f'f_{col}'] = getattr(f, col)
                        record[f'r_{col}'] = getattr(r, col)
                    
                    record['f_self_comp'] = self_complementarity(f.sequence)
                    record['r_self_comp'] = self_complementarity(r.sequence)
                    comp_score = pair_complementarity(f.sequence, r.sequence)
                    record['pair_comp'] = comp_score
                    
                    # Only include pairs with low complementarity (<=6)
                    if comp_score <= 6:
                        pair_records.append(record)

        # Write out pairings for the current gene
        out_path = PAIRINGS_DIR / f"{gene}.csv"
        if pair_records:
            pd.DataFrame(pair_records).to_csv(out_path, index=False)
            print(f"✅ Wrote {len(pair_records)} primer pairs for {gene} to {out_path}")
        else:
            print(f"  - No valid primer pairs found for {gene}")

print("\nAll species have been processed.")
