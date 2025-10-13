# 02_blast_primers.py
import subprocess, os, pandas as pd, re, textwrap, tempfile
from typing import Dict, Any, Optional, List, Tuple, Set
from tqdm import tqdm
import re, gzip, requests, pathlib
from collections import Counter
import numpy as np

# --- SETTINGS ---
BLAST_DB = "blast_dbs/human/refseq_rna"
THREADS = 24
FIELD_SPEC = "6 qseqid sacc stitle staxids bitscore evalue qcovs"
PRIMER_SUMMARY_FILE = "primer_3.csv"
BLAST_RESULTS_FILE = "primer_3_blast.csv"
# -----------------

STOP_WORDS = {"RT", "PCR", "DNA", "RNA"}
DNA_LETTERS = set("ACGT")
TOKEN_RX = re.compile(r"[A-Za-z0-9\-]{3,10}")

# --- Helper Functions (No changes needed here) ---

def load_hgnc_symbol_set(cache_dir="hgnc_cache") -> Set[str]:
    cache_path = pathlib.Path(cache_dir) / "homo_sapiens.gene_info.gz"
    if not cache_path.exists():
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        url = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
        print("⏬ downloading HGNC symbol list (Homo sapiens)...")
        r = requests.get(url, timeout=60)
        r.raise_for_status()
        cache_path.write_bytes(r.content)
    symbols = set()
    with gzip.open(cache_path, "rt", encoding='utf-8') as fh:
        for line in fh:
            if line.startswith("#"): continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) > 2:
                symbols.add(cols[2].upper())
    print(f"Loaded {len(symbols)} symbols into HGNC set.")
    return symbols

def extract_symbols_from_title(title: str, hgnc_set: Set[str]) -> List[str]:
    symbols_found = []
    common_species_pattern = r"\b(Homo sapiens|Mus musculus|Rattus norvegicus)\b"
    cleaned_title = re.sub(common_species_pattern, "", title, flags=re.IGNORECASE)
    m = re.search(r"\(([A-Z0-9\-]{3,10})\)", cleaned_title, re.IGNORECASE)
    if m:
        tok = m.group(1).upper()
        if tok in hgnc_set and tok not in STOP_WORDS and not tok[0].isdigit() and not (set(tok) <= DNA_LETTERS):
            symbols_found.append(tok)
    for raw in TOKEN_RX.findall(cleaned_title):
        tok = raw.upper()
        if tok in STOP_WORDS or tok[0].isdigit() or tok not in hgnc_set or set(tok) <= DNA_LETTERS:
            continue
        symbols_found.append(tok)
    return symbols_found

def determine_best_gene_from_blast_hits(qid_group: pd.DataFrame, hgnc_set: Set[str]) -> Tuple[Optional[str], Optional[float]]:
    if qid_group.empty: return None, None
    qid_group["qcovs"] = pd.to_numeric(qid_group["qcovs"], errors='coerce')
    hits_100_qcov = qid_group[qid_group["qcovs"] == 100].copy()
    if hits_100_qcov.empty: return None, qid_group["evalue"].min()
    all_symbols = [symbol for title in hits_100_qcov["stitle"] for symbol in extract_symbols_from_title(title, hgnc_set)]
    best_evalue_100_qcov = hits_100_qcov["evalue"].min()
    if not all_symbols: return None, best_evalue_100_qcov
    most_common_symbol = Counter(all_symbols).most_common(1)[0][0]
    return most_common_symbol, best_evalue_100_qcov

def find_first_hgnc_in_author_list(author_genes_str: Optional[str], hgnc_set: Set[str]) -> Optional[str]:
    if pd.isna(author_genes_str) or not isinstance(author_genes_str, str): return None
    for entry in author_genes_str.split(','):
        match = re.match(r"([A-Za-z0-9\-]+)", entry.strip())
        if match:
            symbol_candidate_upper = match.group(1).upper()
            if symbol_candidate_upper == "18S": return "18S"
            if symbol_candidate_upper in hgnc_set: return symbol_candidate_upper
    return None

# --- Main Logic ---
HGNC = load_hgnc_symbol_set()

try:
    summary_df = pd.read_csv(PRIMER_SUMMARY_FILE)
    print(f"Loaded {len(summary_df)} total entries from '{PRIMER_SUMMARY_FILE}'.")
except FileNotFoundError:
    print(f"Error: '{PRIMER_SUMMARY_FILE}' not found. Please ensure this file exists.")
    exit()

if os.path.exists(BLAST_RESULTS_FILE):
    print(f"Loading existing BLAST results from '{BLAST_RESULTS_FILE}'.")
    blast_cols_to_load = ['sequence', 'gene', 'evalue']
    existing_blast_df = pd.read_csv(BLAST_RESULTS_FILE, usecols=lambda c: c in blast_cols_to_load)
    blasted_sequences = set(existing_blast_df['sequence'])
    print(f"Found {len(blasted_sequences)} previously blasted sequences.")
else:
    print("No existing BLAST results file found. All sequences will be blasted.")
    existing_blast_df = pd.DataFrame(columns=['sequence', 'gene', 'evalue'])
    blasted_sequences = set()

all_sequences = set(summary_df['sequence'])
new_sequences_to_blast = all_sequences - blasted_sequences

if not new_sequences_to_blast:
    print("✅ No new sequences to BLAST. Output file is already up to date.")
    print("Re-merging data to ensure counts and other metadata are current...")
    final_df = summary_df.merge(existing_blast_df, on="sequence", how="left")
else:
    print(f"Found {len(new_sequences_to_blast)} new unique sequences to BLAST.")
    id2seq = {}
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".fa", encoding='utf-8') as temp_fasta_file:
        for i, seq in enumerate(new_sequences_to_blast):
            qid = f"q{i}"
            id2seq[qid] = seq
            temp_fasta_file.write(f">{qid}\n{seq}\n")
    
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".b6") as temp_blast_out_file:
        pass

    cmd = [
        "blastn", "-task", "blastn-short", "-query", temp_fasta_file.name,
        "-db", BLAST_DB, "-out", temp_blast_out_file.name, "-outfmt", FIELD_SPEC,
        "-max_target_seqs", "25", "-num_threads", str(THREADS), "-dust", "no",
        "-soft_masking", "false",
    ]
    print("Running BLAST on new sequences...")
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True, encoding='utf-8')
        print("BLAST done ✔")
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"BLAST command failed: {e}")
        os.remove(temp_fasta_file.name)
        os.remove(temp_blast_out_file.name)
        exit()

    try:
        hits_raw = pd.read_csv(temp_blast_out_file.name, sep="\t", names=["qid", "sacc", "stitle", "staxids", "bitscore", "evalue", "qcovs"])
    except pd.errors.EmptyDataError:
        hits_raw = pd.DataFrame()
    
    new_blast_results = []
    qids_with_hits = set(hits_raw['qid']) if not hits_raw.empty else set()

    if not hits_raw.empty:
        for qid, group in tqdm(hits_raw.groupby("qid"), desc="Processing new BLAST hits", total=len(qids_with_hits)):
            gene_symbol, best_evalue = determine_best_gene_from_blast_hits(group, HGNC)
            new_blast_results.append({"sequence": id2seq[qid], "gene": gene_symbol, "evalue": best_evalue})
    
    for qid, seq in id2seq.items():
        if qid not in qids_with_hits:
            new_blast_results.append({"sequence": seq, "gene": None, "evalue": np.nan})

    new_hits_df = pd.DataFrame(new_blast_results)
    
    combined_blast_df = pd.concat([existing_blast_df, new_hits_df], ignore_index=True).drop_duplicates(subset=['sequence'], keep='last')
    final_df = summary_df.merge(combined_blast_df, on="sequence", how="left")

    os.remove(temp_fasta_file.name)
    os.remove(temp_blast_out_file.name)

if 'gene' not in final_df.columns: final_df['gene'] = pd.NA
final_df['gene'] = final_df['gene'].fillna('')

for index, row in tqdm(final_df.iterrows(), total=final_df.shape[0], desc="Applying fallback gene assignment"):
    if not row["gene"]:
        fallback_gene = find_first_hgnc_in_author_list(row.get("Source_Genes"), HGNC)
        if fallback_gene:
            final_df.loc[index, "gene"] = fallback_gene

# --- 9. Sort and save the final, complete file ---
final_df = final_df.sort_values(["Occurrence_Count", "evalue"], ascending=[False, True]).reset_index(drop=True)

# ✅ FIX: The order of operations is corrected here.
# First, rename the columns to their final names.
final_df.rename(columns={'Source_Genes': 'author_gene_list', 'Source_PMCIDs': 'pmcid_list'}, inplace=True)

# Second, define the list of columns to keep, using the NEW names.
final_columns = ["sequence", "gene", "author_gene_list", "evalue", "Occurrence_Count", "pmcid_list", "Citation_Count"]
        
# Third, ensure all desired columns exist in the DataFrame.
for col in final_columns:
    if col not in final_df.columns:
        final_df[col] = np.nan

# Finally, select and reorder the columns for the output file.
final_df = final_df[final_columns]

final_df.to_csv(BLAST_RESULTS_FILE, index=False)
print(f"✅ Saved updated results to '{BLAST_RESULTS_FILE}'")