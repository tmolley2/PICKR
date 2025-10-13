# 07_map_pmcids.py (Corrected)
import os
import glob
import csv
from collections import defaultdict
import pandas as pd
from tqdm import tqdm
from pathlib import Path

# --------------------------------------------------------------------------------------
# 1) Configuration
# --------------------------------------------------------------------------------------
SPECIES_LIST = ["frog", "cow", "human", "monkey", "rat", "pig", "mouse", "zebrafish", "worm", "fly"]

# Base Directories
BASE_TEMP_DIR = Path('temp')

# --------------------------------------------------------------------------------------
# 2) Build mapping from sequence → set of PMCIDs (Loaded ONCE)
# --------------------------------------------------------------------------------------
print("Loading raw primer data...")
seq_to_pmcids = defaultdict(set)
raw_primer_file = 'raw_primers.csv'

if not os.path.exists(raw_primer_file):
    print(f"Error: '{raw_primer_file}' not found. This file is required to proceed.")
    exit()

with open(raw_primer_file, newline='', encoding='utf-8', errors='replace') as raw_file:
    reader = csv.DictReader(raw_file)
    for row in tqdm(reader, desc="Loading raw primers"):
        seq = row.get('Sequence', '').strip().upper()
        pmcid = row.get('PMCID', '').strip()
        if seq and pmcid:
            seq_to_pmcids[seq].add(pmcid)
print("Raw primer data loaded.")

# --------------------------------------------------------------------------------------
# 3) Helper function
# --------------------------------------------------------------------------------------
_complement = str.maketrans('ACGT', 'TGCA')
def reverse_complement(seq):
    """Computes the reverse-complement of a DNA sequence."""
    if not isinstance(seq, str): return ""
    return seq.translate(_complement)[::-1]

# --------------------------------------------------------------------------------------
# 4) Main Processing Loop for Each Species
# --------------------------------------------------------------------------------------
for species in SPECIES_LIST:
    print(f"\n{'='*20} PROCESSING ORIENTATIONS FOR: {species.upper()} {'='*20}")

    input_folder = BASE_TEMP_DIR / f'cross_specificity_{species}'
    output_folder = BASE_TEMP_DIR / f'primer_pairings_shared_counts_{species}'

    if not input_folder.is_dir():
        print(f"Warning: Input directory '{input_folder}' not found. Skipping {species}.")
        continue
    
    output_folder.mkdir(exist_ok=True)

    file_list = glob.glob(os.path.join(input_folder, '*.csv'))

    for file_path in tqdm(file_list, desc=f"Processing {species} files"):
        try:
            df = pd.read_csv(file_path, dtype=str).fillna('')
        except pd.errors.EmptyDataError:
            print(f"Skipping empty file: {os.path.basename(file_path)}")
            continue
        
        # Prepare lists for all new columns
        f_reg_pmcid_list, f_inv_pmcid_list = [], []
        r_reg_pmcid_list, r_inv_pmcid_list = [], []
        f_reg_pmcid_count, f_inv_pmcid_count = [], []
        r_reg_pmcid_count, r_inv_pmcid_count = [], []
        pair_shared_pmcid_list, pair_shared_pmcid_count = [], []

        for _, row in df.iterrows():
            f_seq = row.get('f_sequence', '').strip().upper()
            r_seq = row.get('r_sequence', '').strip().upper()

            # --- ✅ CORRECTED LOGIC ---
            # Directly use the global seq_to_pmcids map as the source of truth.
            # Do not rely on `f_pmcid_list` or `r_pmcid_list` from the input file for this step.
            
            f_rc = reverse_complement(f_seq)
            r_rc = reverse_complement(r_seq)
            
            # Get the complete sets of PMCIDs for regular and inverted sequences
            f_reg_set = seq_to_pmcids.get(f_seq, set())
            f_inv_set = seq_to_pmcids.get(f_rc, set())
            r_reg_set = seq_to_pmcids.get(r_seq, set())
            r_inv_set = seq_to_pmcids.get(r_rc, set())

            # Convert to sorted lists for consistent output
            f_reg = sorted(list(f_reg_set))
            f_inv = sorted(list(f_inv_set))
            r_reg = sorted(list(r_reg_set))
            r_inv = sorted(list(r_inv_set))
            # --- END OF CORRECTION ---

            f_reg_pmcid_list.append(','.join(f_reg))
            f_inv_pmcid_list.append(','.join(f_inv))
            r_reg_pmcid_list.append(','.join(r_reg))
            r_inv_pmcid_list.append(','.join(r_inv))
            
            f_reg_pmcid_count.append(len(f_reg))
            f_inv_pmcid_count.append(len(f_inv))
            r_reg_pmcid_count.append(len(r_reg))
            r_inv_pmcid_count.append(len(r_inv))
            
            # This logic was already correct, using the global map
            shared = seq_to_pmcids.get(f_seq, set()) & seq_to_pmcids.get(r_seq, set())
            pair_shared_pmcid_list.append(','.join(sorted(shared)))
            pair_shared_pmcid_count.append(len(shared))

        # Attach all new columns to the DataFrame
        df['f_reg_pmcid_list'] = f_reg_pmcid_list
        df['f_inv_pmcid_list'] = f_inv_pmcid_list
        df['r_reg_pmcid_list'] = r_reg_pmcid_list
        df['r_inv_pmcid_list'] = r_inv_pmcid_list
        
        df['f_reg_pmcid_count'] = f_reg_pmcid_count
        df['f_inv_pmcid_count'] = f_inv_pmcid_count
        df['r_reg_pmcid_count'] = r_reg_pmcid_count
        df['r_inv_pmcid_count'] = r_inv_pmcid_count
        
        df['pair_shared_pmcid_list'] = pair_shared_pmcid_list
        df['pair_shared_pmcid_count'] = pair_shared_pmcid_count

        # Save the fully processed file to the final output folder
        out_path = output_folder / os.path.basename(file_path)
        df.to_csv(out_path, index=False)

print(f"\n✅ All species have been processed.")
