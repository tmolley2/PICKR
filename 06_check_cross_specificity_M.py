# 06_check_cross_specificity.py

import os
import subprocess
import pandas as pd
import re
from tqdm import tqdm
import tempfile
from pathlib import Path

# --- Configuration ---
SPECIES_LIST = ["cow", "monkey", "rat", "pig", "mouse"]

# Base Directories
BASE_TEMP_DIR = Path('temp')
BASE_BLAST_DIR = Path('blast_dbs/blast_cross')
BASE_SCORED_DIR = Path('scored') # A base directory for scored files

# --- Performance Settings ---
# Set this to the number of CPU cores you want to use for the BLAST search
NUM_THREADS = 24

# ----------------------------------------------------------------------------------------------------------------------------------
# Utility Function (No changes needed)
# ----------------------------------------------------------------------------------------------------------------------------------
def extract_gene_symbol(desc: str) -> str:
    """Tries multiple patterns to find the gene symbol in the header."""
    match = re.search(r'\[Gene=([^\]]+)\]', desc, re.IGNORECASE)
    if match: return match.group(1)
    match = re.search(r'gene[="\']([^"\']+)["\']?', desc, re.IGNORECASE)
    if match: return match.group(1)
    match = re.search(r'\(([^,)]+)\)', desc)
    if match: return match.group(1).split()[0]
    return None

# ----------------------------------------------------------------------------------------------------------------------------------
# Main Execution Loop
# ----------------------------------------------------------------------------------------------------------------------------------

# ✅ NEW: Main loop to iterate over each species
for species in SPECIES_LIST:
    print(f"\n{'='*20} PROCESSING CROSS-SPECIFICITY FOR: {species.upper()} {'='*20}")

    # --- Dynamically set paths for the current species ---
    PAIRINGS_DIR = BASE_TEMP_DIR / f'primer_pairings_{species}'
    OUTPUT_DIR = BASE_TEMP_DIR / f'cross_specificity_{species}'
    EXISTING_SCORES_DIR = BASE_SCORED_DIR / f'scored_{species}' # Using a consistent naming scheme
    
    # Paths for the species-specific BLAST database
    SPECIES_BLAST_DIR = BASE_BLAST_DIR / species
    FASTA_PATH = SPECIES_BLAST_DIR / f'{species}.rna.fna'
    DB_PREFIX = SPECIES_BLAST_DIR / f'{species}_refseq_rna'

    # --- Create directories if they don't exist ---
    if not PAIRINGS_DIR.is_dir():
        print(f"Warning: Primer pairings directory not found at '{PAIRINGS_DIR}'. Skipping {species}.")
        continue
    OUTPUT_DIR.mkdir(exist_ok=True)
    EXISTING_SCORES_DIR.mkdir(exist_ok=True)
    SPECIES_BLAST_DIR.mkdir(exist_ok=True)

    # --- STEP 1: BUILD BLAST DB IF NEEDED ---
    if not FASTA_PATH.is_file():
        print(f"Error: FASTA file for BLAST not found at '{FASTA_PATH}'. Skipping {species}.")
        continue
    
    if not (DB_PREFIX.with_suffix('.nsq')).is_file():
        print(f"Building BLAST database for {species}...")
        subprocess.run(['makeblastdb', '-in', str(FASTA_PATH), '-dbtype', 'nucl', '-out', str(DB_PREFIX)], check=True)
    else:
        print(f"BLAST database for {species} already exists.")

    # --- STEP 2: INDEX ACCESSION → GENE SYMBOL ---
    print(f"Building gene map for {species}...")
    tx2gene = {}
    with open(FASTA_PATH) as fasta:
        for line in fasta:
            if not line.startswith('>'): continue
            header = line[1:].strip()
            acc = header.split()[0].split('.')[0]
            tx2gene[acc] = extract_gene_symbol(header) or acc
    print(f"Gene map for {species} created with {len(tx2gene)} entries.")

    # --- STEP 3A: PRE-SCAN - Identify ALL unique primers for the current species ---
    all_primers_to_check = set()
    files_to_process = [fn for fn in os.listdir(PAIRINGS_DIR) if fn.endswith('.csv')]

    for fn in tqdm(files_to_process, desc=f"Pre-scanning files for {species}"):
        pairings_df = pd.read_csv(PAIRINGS_DIR / fn, dtype=str).fillna('')
        all_primers_to_check.update(pairings_df['f_sequence'].str.upper().str.strip().dropna())
        all_primers_to_check.update(pairings_df['r_sequence'].str.upper().str.strip().dropna())

    print(f"Found {len(all_primers_to_check)} unique primers for {species}.")
    if not all_primers_to_check:
        print(f"No primers to check for {species}. Moving to next species.")
        continue

    # --- STEP 3B: BATCH BLAST - Run BLAST once for all unique primers ---
    blast_results_cache = {}
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".fa") as tmp_fasta, \
         tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".tsv") as tmp_blast_out:
        
        for seq in all_primers_to_check:
            tmp_fasta.write(f">{seq}\n{seq}\n")
        
        fasta_path_str = tmp_fasta.name
        blast_out_path_str = tmp_blast_out.name

    cmd = [
        'blastn', '-task', 'blastn-short', '-db', str(DB_PREFIX),
        '-query', fasta_path_str, '-out', blast_out_path_str,
        '-qcov_hsp_perc', '100', '-outfmt', '6 qseqid sseqid mismatch sseq',
        '-num_threads', str(NUM_THREADS)
    ]
    
    print(f"🚀 Running batch BLAST for {species}...")
    subprocess.run(cmd, check=True)
    print("✅ Batch BLAST complete.")

    # --- STEP 3C: CACHE RESULTS ---
    print("Caching BLAST results...")
    blast_output_df = pd.read_csv(blast_out_path_str, sep='\t', names=['qseqid', 'sseqid', 'mismatch', 'sseq'])
    
    for seq, group in tqdm(blast_output_df.groupby('qseqid'), desc="Parsing results"):
        mm_hits = {i: set() for i in range(5)}
        mm_pats = {i: {} for i in range(5)}
        for _, row in group.iterrows():
            mm = int(row['mismatch'])
            if mm > 4: continue
            gene = tx2gene.get(row['sseqid'].split('.')[0], row['sseqid'])
            mm_hits[mm].add(gene)
            if mm >= 1:
                diffs = [i for i, (a, b) in enumerate(zip(row['qseqid'], row['sseq'])) if a != b]
                if len(diffs) == mm:
                    pat = ['.'] * len(row['qseqid'])
                    for d in diffs: pat[d] = row['sseq'][d]
                    mm_pats[mm][gene] = ''.join(pat)
        blast_results_cache[seq] = (mm_hits, mm_pats)

    os.remove(fasta_path_str)
    os.remove(blast_out_path_str)

    # --- STEP 4: PROCESS PRIMER PAIRS USING THE CACHE ---
    for fn in tqdm(files_to_process, desc=f"Assembling files for {species}"):
        pairings_df = pd.read_csv(PAIRINGS_DIR / fn, dtype=str).fillna('')
        pairings_df['f_sequence'] = pairings_df['f_sequence'].str.upper().str.strip()
        pairings_df['r_sequence'] = pairings_df['r_sequence'].str.upper().str.strip()

        gene_target = fn.replace('.csv', '').upper()
        # Handle Xenopus .L/.S suffixes in the target name
        if '.L' in gene_target or '.S' in gene_target:
            gene_target = gene_target.rsplit('.', 1)[0]

        for index, row in pairings_df.iterrows():
            fseq, rseq = row['f_sequence'], row['r_sequence']
            
            f_mm_hits, f_mm_pats = blast_results_cache.get(fseq, ({i:set() for i in range(5)}, {i:{} for i in range(5)}))
            r_mm_hits, r_mm_pats = blast_results_cache.get(rseq, ({i:set() for i in range(5)}, {i:{} for i in range(5)}))

            for mm in range(5):
                # Forward primer
                f_hits_current = {g for g in f_mm_hits[mm] if g and g.upper().split('.')[0] != gene_target}
                f_pats_current = {g: p for g, p in f_mm_pats[mm].items() if g and g.upper().split('.')[0] != gene_target}
                pairings_df.loc[index, f'f_cross_specificity_{mm}'] = ','.join(sorted(f_hits_current))
                pairings_df.loc[index, f'f_cross_specificity_{mm}_seq'] = ';'.join(f_pats_current.values()) if mm >= 1 else ''
                
                # Reverse primer
                r_hits_current = {g for g in r_mm_hits[mm] if g and g.upper().split('.')[0] != gene_target}
                r_pats_current = {g: p for g, p in r_mm_pats[mm].items() if g and g.upper().split('.')[0] != gene_target}
                pairings_df.loc[index, f'r_cross_specificity_{mm}'] = ','.join(sorted(r_hits_current))
                pairings_df.loc[index, f'r_cross_specificity_{mm}_seq'] = ';'.join(r_pats_current.values()) if mm >= 1 else ''

        pairings_df.to_csv(OUTPUT_DIR / fn, index=False)

print("\n--- All species cross-specificity checks complete ---")
