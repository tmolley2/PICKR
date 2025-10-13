# 03_calculate_tm.py

import os, re
import pandas as pd
from Bio.SeqUtils.MeltingTemp import Tm_NN          # Biopython ≥1.83
from Bio.SeqUtils import gc_fraction
from pathlib import Path
from tqdm import tqdm

# -------- file paths --------
SRC  = Path("primer_3_blast.csv")   # adjust if needed
OUTDIR   = Path("temp/gene_csvs")            # all per-gene files go here
GENE_COL = "gene"                # rename if your header differs
ID_COL   = "pmcid_list"

OUTDIR.mkdir(exist_ok=True)

# -------- load & annotate --------
df = pd.read_csv(SRC)

# 2) group & count
counts = (
    df.groupby("gene")
      .size()
      .reset_index(name="n_sequences")
      .sort_values("n_sequences", ascending=False)
)
counts.to_csv("temp/sequence_count_per_gene.csv", index=False)

# ---3 calcualte melt temperature
def calc_tm_gc(seq: str):
    tm = Tm_NN(seq, Na=50, Mg=1.5, dNTPs=0.2)
    gc = gc_fraction(seq) * 100          # gc_fraction returns 0–1
    return round(tm, 2), round(gc, 1)

df[["Tm_C", "GC_pct"]] = df["sequence"].apply(
    lambda s: pd.Series(calc_tm_gc(s))
)

# ── 4. Compute “times used” = # papers for each primer row
df["pmcid_count"] = df[ID_COL].fillna("").apply(
    lambda cell: 0 if cell == "" else len(re.split(r"[;,]", cell))
)

# ── 4. Split by gene, keep ≤20 most-used, write CSV
def safe_name(gene: str) -> str:
    """Remove characters not allowed in filenames (Windows safe too)."""
    return re.sub(r"[^\w\-]", "_", gene)

for gene, sub in tqdm(df.groupby(GENE_COL), desc="Genes"):
    top = (
        sub.sort_values("pmcid_count", ascending=False)
           .head(40)                           # ≤ 20 rows
           .reset_index(drop=True)
    )
    out_file = OUTDIR / f"{safe_name(gene)}.csv"
    top.to_csv(out_file, index=False)

print(f"✓ Per-gene CSVs written to folder '{OUTDIR}'")