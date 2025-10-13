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

The access teh tool, www.primerpickr.com hosts the active database. To regerenate the databse yourself, refere tot eh codes outlined int his github and use the following: 

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
