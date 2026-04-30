"""
generate_gene_positions.py

Run this script ONCE on your Mac to generate hg38_gene_positions.csv
which is required for the infercnvpy CNV inference notebook.

Usage:
    conda activate arc-breast
    pip install pybiomart
    python scripts/generate_gene_positions.py

Output:
    data/raw/hg38_gene_positions.csv (~20,000 genes with GRCh38 coordinates)

This file is gitignored (too large for GitHub).
Runtime: 2-5 minutes depending on Ensembl server speed.
"""

from pathlib import Path
import sys

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from src import config as cfg

OUTPUT_PATH = cfg.RAW_DIR / "hg38_gene_positions.csv"

if OUTPUT_PATH.exists():
    print(f"Already exists: {OUTPUT_PATH}")
    print("Delete it and re-run if you want to regenerate.")
    sys.exit(0)

print("Fetching GRCh38 gene coordinates from Ensembl biomart...")
print("This takes 2-5 minutes depending on server speed.")
print()

try:
    import pybiomart
except ImportError:
    print("pybiomart not installed. Run:")
    print("  pip install pybiomart")
    sys.exit(1)

import pandas as pd

# Connect to Ensembl
dataset = pybiomart.Dataset(
    name='hsapiens_gene_ensembl',
    host='http://www.ensembl.org',
)

print("Connected to Ensembl. Fetching gene positions...")

# Fetch HGNC symbol + genomic coordinates
coords = dataset.query(attributes=[
    'hgnc_symbol',
    'chromosome_name',
    'start_position',
    'end_position',
    'gene_biotype',
])
coords.columns = ['gene_name', 'chromosome', 'start', 'end', 'biotype']

print(f"Raw fetch: {len(coords)} entries")

# Filter to standard chromosomes only (1-22, X, Y)
standard_chroms = [str(i) for i in range(1, 23)] + ['X', 'Y']
coords = coords[coords['chromosome'].isin(standard_chroms)].copy()

# Add 'chr' prefix to match infercnvpy format
coords['chromosome'] = 'chr' + coords['chromosome'].astype(str)

# Keep only protein-coding and lncRNA genes (most informative for CNV)
keep_biotypes = ['protein_coding', 'lncRNA', 'antisense_RNA']
coords = coords[coords['biotype'].isin(keep_biotypes)].copy()

# Remove entries without gene symbol
coords = coords[coords['gene_name'].notna() & (coords['gene_name'] != '')].copy()

# Remove duplicates — keep first entry per gene symbol
coords = coords.drop_duplicates(subset='gene_name', keep='first')

# Sort by chromosome and position
chrom_order = {f'chr{i}': i for i in range(1, 23)}
chrom_order['chrX'] = 23
chrom_order['chrY'] = 24
coords['chrom_num'] = coords['chromosome'].map(chrom_order)
coords = coords.sort_values(['chrom_num', 'start']).drop(columns=['chrom_num', 'biotype'])

print(f"After filtering: {len(coords)} genes")
print(f"Chromosomes covered: {sorted(coords['chromosome'].unique())}")

# Verify ARC is present
arc = coords[coords['gene_name'] == 'ARC']
if len(arc) > 0:
    print(f"\nARC coordinates confirmed: {arc[['chromosome','start','end']].to_dict('records')[0]}")
else:
    print("\nWARNING: ARC not found in coordinate file — it may be missing from Ensembl HGNC symbols")

# Save
OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
coords.to_csv(OUTPUT_PATH, index=False)
print(f"\nSaved: {OUTPUT_PATH}")
print(f"File size: {OUTPUT_PATH.stat().st_size / 1e6:.1f} MB")
print("\nDone. You can now run notebook 06b.")
