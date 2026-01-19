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
