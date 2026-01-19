# PICKR
Primer PICKR codebase
# Primer PICKR — Publication Integrations for Composite-Knowledge Ranking

**Primer PICKR** is a large-scale, open, continuously updated database that extracts and ranks >7 million published PCR primers from >400 000 peer-reviewed studies to generate experimentally informed primer pairs across ten model organisms.  
🌐 Public web access: [https://www.primerpickr.com](https://www.primerpickr.com)

This repository accompanies our manuscript: "Primer PICKR: Literature-Mined Scoring Platform for Robust RT-qPCR Primers".

It provides the source code, schema, and documentation for data parsing, alignment, and scoring workflows underlying the generation online database.

---

## 1. System Requirements
Primer PICKR was developed and tested on:

| Component | Recommended | Notes |
|------------|--------------|-------|
| OS | Linux (Ubuntu 20.04 +) / macOS 12 + / Windows 10 + | Tested on Ubuntu 22.04 LTS |
| Python | ≥ 3.9 | core scripts use `pandas`, `Biopython`, `requests`, `numpy`, `sqlalchemy` |
| Memory | ≥ 16 GB RAM | bulk sequence alignment requires several GB |
| Disk | ≥ 50 GB free space | for full dataset re-build |
| Optional tools | BLAST+ (2.13 +) | for transcript alignment |

## 2. Installation Guide

The access teh tool, www.primerpickr.com hosts the active database. To regerenate the database yourself, refere tot eh codes outlined int his github and use the following: 

bash
git clone https://github.com/your-org/primer-pickr.git
cd primer-pickr
python3 -m venv venv
source venv/bin/activate


iinstall the following packages:

pandas

os

glob

pathlib

tqdm

gzip

numpy

scipy

bio.seq

re
### Installing BLAST+ Locally

This project relies on NCBI BLAST+ for sequence alignment and primer validation. BLAST+ must be installed locally and accessible from your command line.

1. Check if BLAST is Already Installed

Open a terminal (or Anaconda Prompt on Windows) and run:

blastn -version


If BLAST is installed, you should see version information. If not, follow the instructions below.

macOS

### Option A: Install via Homebrew (recommended)

If you don’t have Homebrew, install it first:

/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"


Then install BLAST:

brew install blast


Verify installation:

blastn -version


### Option B: Manual Installation

Download the macOS BLAST+ binaries from NCBI:
https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

Extract the archive:

tar -xzf ncbi-blast-*-x64-macosx.tar.gz


Add BLAST to your PATH (replace path as needed):

export PATH=$PATH:/path/to/ncbi-blast-*/bin


Add this line to ~/.zshrc or ~/.bashrc to make it permanent.

## 3. Demo

To demo the databse, please go to www.primerpickr.com to access the tool.

## 4. Instructions for Use

a. Web Interface

1.) Visit primerpickr.com
2.) Enter a gene symbol or accession to retrieve ranked primer pairs.
3.) Filter by organism, amplicon length, and score threshold.
4.) Copy outputs into a CSV for integration into qPCR assay design pipelines.

Notes for Reviewers and Editors

The web instance (https://www.primerpickr.com) provides all functionalities for end users.
The GitHub repository contains reproducible scripts and demo data to satisfy transparency requirements without requiring reviewers to rebuild the full database.

## 5. Expected Results / Outputs

Running the PICKR pipeline produces a set of **gene-specific primer candidates**, **ranked primer pairs**, and **annotated CSV outputs** that can be used directly for downstream primer selection and validation. The exact number of hits per gene will vary depending on how frequently that gene appears in the literature and how stringent your filtering thresholds are.

### Summary of Outputs (What You Should Expect)

After completing the full workflow, you should expect:

- **One annotated primer file per gene**
  - Primers are mapped to a reference transcript
  - Orientation and reference coordinates are appended
  - Unmappable primers are filtered out and logged separately

- **One primer-pair file per gene**
  - All valid forward–reverse primer combinations are enumerated
  - Amplicon distance constraints (e.g., <300 bp) are enforced
  - High complementarity primer pairs (likely primer-dimers) are removed

- **Optional cross-species scoring outputs**
  - Each primer pair receives a score for predicted cross-reactivity in multiple species
  - Useful for flagging primers likely to amplify orthologs or off-target species templates

---

### Example Output Files

Depending on which scripts you run, you will typically generate files in the following folders:

#### 1) Annotated primers (per gene)
**Folder:** `annotated_with_complements/`  
**Example file:** `GAPDH_annotated.csv`

Expected columns include:
- `sequence`
- `match_start`, `match_end`
- `orientation` (`forward` or `reverse_complement`)
- (plus any original metadata columns from the literature extraction step)

Unmatched primers (no valid mapping) are saved separately as:
- `annotated_with_complements/<gene>_skipped_primers.csv`

Genes that could not be processed (e.g., missing reference mRNA or insufficient entries) are logged in:
- `skipped_genes/skipped_genes.csv`

---

#### 2) Primer pair enumeration + filtering (per gene)
**Folder:** `primer_pairings/`  
**Example file:** `GAPDH_pairs.csv`

Expected columns include:
- `f_sequence`, `r_sequence`
- `f_match_start`, `r_match_end`
- `f_self_comp`, `r_self_comp` (self-dimer complementarity scores)
- `pair_comp` (pairwise complementarity score)
- All original primer metadata columns preserved with prefixes:
  - `f_<metadata_column>`
  - `r_<metadata_column>`

Only primer pairs passing constraints (e.g., amplicon length <300 bp and pair complementarity ≤6) are retained.

---

#### 3) Cross-species reactivity scores (optional)
**Folder:** `species_crosschecked/`  
**Example file:** `GAPDH_species_scores.csv`

Expected columns include one new score per species, such as:
- `mouse`
- `rat`
- `zebrafish`
- `drosophila`
- `c_elegans`
- `e_coli`


Species scores are computed per primer pair by averaging the forward and reverse match scores.

---

### What “Good” Results Look Like

For commonly studied genes (e.g., housekeeping genes), you should typically observe:
- Multiple primer candidates mapping cleanly to the reference
- Dozens to hundreds of valid primer pairs per gene before filtering
- A smaller subset of primer pairs remaining after removing:
  - invalid orientations
  - excessive amplicon length
  - high primer–primer complementarity (primer dimer risk)

For less-studied genes, you may observe:
- Few or no primer candidates extracted from literature
- Missing or ambiguous reference sequences
- Limited valid pairings

These cases are expected and are logged in the `skipped_genes/` and `*_skipped_primers.csv` outputs for review.


### An example output file is provide in the folder as OUTPUT.csv
---

### Notes on Reproducibility

Output files are deterministic given the same:
- input CSVs
- reference sequences
- filtering thresholds

However, results may differ if:
- reference FASTA files are updated
- different transcript sets are used (e.g., RefSeq vs Ensembl)
- mismatch thresholds or amplicon constraints are changed

