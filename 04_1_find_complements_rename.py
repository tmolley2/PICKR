import pandas as pd
from pathlib import Path
from tqdm import tqdm

# --- Configuration ---

# 1. List of species corresponding to your BLAST output files.
SPECIES_LIST = ["frog", "fly", "worm", "zebrafish"]

# 2. The directory where your 'primer_3_blast_*.csv' files are located.
#    Assumes the script is in the same directory as the files.
INPUT_DIR = Path(".") 

# 3. The name of the NCBI annotation file you downloaded.
NCBI_ANNOTATION_FILE = "gene2refseq.gz"

# 4. The name of the column in your CSVs that contains the accession number (e.g., 'XM_001234.1').
ACCESSION_COLUMN = 'gene1'

# --- Main Script ---

def create_accession_to_symbol_map(annotation_file: Path) -> dict:
    """
    Reads the NCBI gene2refseq.gz file and creates a mapping dictionary
    from RNA accession number to official gene symbol.
    """
    print(f"✔️ Loading NCBI annotation file: {annotation_file}")
    try:
        # Define the columns we need to speed up loading
        columns_to_use = ['RNA_nucleotide_accession.version', 'Symbol']
        
        df_ncbi = pd.read_csv(
            annotation_file,
            sep='\t',
            compression='gzip',
            usecols=columns_to_use,
            na_values=['-'] # Treat '-' as Not a Number (NaN)
        )
        
        # Drop rows where accession or symbol is missing
        df_ncbi.dropna(subset=columns_to_use, inplace=True)
        
        # Create the dictionary for mapping
        # This handles cases where one accession might appear multiple times by keeping the first symbol
        accession_map = df_ncbi.drop_duplicates(subset=['RNA_nucleotide_accession.version']).set_index('RNA_nucleotide_accession.version')['Symbol'].to_dict()
        
        print(f"✔️ Created mapping for {len(accession_map):,} unique accessions.")
        return accession_map

    except FileNotFoundError:
        print(f"❌ ERROR: Annotation file not found at '{annotation_file}'.")
        print("Please download 'gene2refseq.gz' from the NCBI FTP site.")
        return None
    except Exception as e:
        print(f"❌ An error occurred while reading the NCBI file: {e}")
        return None

def main():
    """
    Main function to process each species file and add gene symbols.
    """
    accession_map = create_accession_to_symbol_map(Path(NCBI_ANNOTATION_FILE))
    
    if accession_map is None:
        print("Exiting due to error.")
        return

    print("\nProcessing species files...")
    for species in tqdm(SPECIES_LIST, desc="Annotating Species"):
        input_file = INPUT_DIR / f"primer_3_blast_{species}.csv"
        output_file = INPUT_DIR / f"primer_3_blast_{species}_with_symbols.csv"

        if not input_file.exists():
            print(f"\n⚠️ Warning: Input file not found for {species}, skipping: {input_file}")
            continue

        # Load the BLAST results from the previous script
        df_species = pd.read_csv(input_file)

        # Create the new 'gene_symbol' column by mapping from the accession column
        df_species['gene_symbol'] = df_species[ACCESSION_COLUMN].map(accession_map)

        # For any accessions that weren't in the NCBI file, fill the 'gene_symbol'
        # with the original accession number so no information is lost.
        df_species['gene_symbol'].fillna(df_species[ACCESSION_COLUMN], inplace=True)
        
        # Save the newly annotated file
        df_species.to_csv(output_file, index=False)

    print(f"\n🎉 Success! All files have been processed. Look for the '_with_symbols.csv' files.")


if __name__ == "__main__":
    main()