"""
Central configuration for the breast cancer ARC scRNA-seq pipeline.

All paths, parameters, and mode flags live here so they can be modified
in one place without touching notebooks or modules.
"""

from pathlib import Path

# =============================================================================
# Execution mode
# =============================================================================
# DEV_MODE = True  -> Subsample to ~40K cells for local Mac execution (8 GB RAM)
# DEV_MODE = False -> Full dataset. Run on Colab / Kaggle / workstation (>=16 GB)
DEV_MODE = True

# Number of cells per patient when subsampling in DEV_MODE.
# 26 patients * 1500 cells = ~39,000 cells total
CELLS_PER_PATIENT_DEV = 1500

# Random seed for reproducibility across sampling, PCA init, clustering
RANDOM_SEED = 42

# =============================================================================
# Paths
# =============================================================================
# Resolve paths relative to the project root (this file's grandparent directory)
PROJECT_ROOT = Path(__file__).resolve().parent.parent

DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"
FIGURES_DIR = PROJECT_ROOT / "figures"
RESULTS_DIR = PROJECT_ROOT / "results"

# Ensure directories exist at import time
for d in (RAW_DIR, PROCESSED_DIR, FIGURES_DIR, RESULTS_DIR):
    d.mkdir(parents=True, exist_ok=True)

# =============================================================================
# GSE176078 file locations
# =============================================================================
# The Wu et al. 2021 paper deposits a single tar archive with the processed
# count matrix + metadata. The archive contains matrix.mtx, features.tsv,
# barcodes.tsv, and metadata.csv files.
GSE_ACCESSION = "GSE176078"
GSE_TAR_FILENAME = "GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz"
GSE_TAR_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176078/suppl/"
    f"{GSE_TAR_FILENAME}"
)

# After extraction, these are the expected files inside data/raw/
RAW_MATRIX = RAW_DIR / "Wu_etal_2021_BRCA_scRNASeq" / "count_matrix_sparse.mtx"
RAW_GENES = RAW_DIR / "Wu_etal_2021_BRCA_scRNASeq" / "count_matrix_genes.tsv"
RAW_BARCODES = RAW_DIR / "Wu_etal_2021_BRCA_scRNASeq" / "count_matrix_barcodes.tsv"
RAW_METADATA = RAW_DIR / "Wu_etal_2021_BRCA_scRNASeq" / "metadata.csv"

# =============================================================================
# Intermediate h5ad files (one per pipeline stage)
# =============================================================================
H5AD_RAW = PROCESSED_DIR / "01_raw.h5ad"
H5AD_QC = PROCESSED_DIR / "02_qc_filtered.h5ad"
H5AD_NORM = PROCESSED_DIR / "03_normalized_hvg.h5ad"
H5AD_INTEGRATED = PROCESSED_DIR / "04_integrated_clustered.h5ad"
H5AD_ANNOTATED = PROCESSED_DIR / "05_annotated.h5ad"

# =============================================================================
# QC thresholds
# =============================================================================
# IMPORTANT: The deposited GSE176078 matrix is already QC-filtered by Wu et al.
# (~100,064 cells passing their original QC). Our QC steps replicate Wu et al.'s
# Methods exactly (>200 genes, >250 UMIs, <20% mito) plus Scrublet doublet
# detection (which Wu did not perform). MAD-based filtering was considered but
# not applied to remain faithful to the original paper's QC.
#
# Wu et al. 2021 (Nature Genetics) Methods, verbatim:
#   "An additional cutoff was applied, filtering for cells with a gene and
#    UMI count greater than 200 and 250, respectively. All cells with a
#    mitochondrial UMI count percentage greater than 20% were removed."
#
# Tumor scRNA-seq tolerates higher mitochondrial % than healthy tissue
# because hypoxic, necrotic, and apoptotic regions are biologically
# meaningful (Luecken & Theis 2019, Mol Syst Biol).

MIN_GENES_PER_CELL = 200
# Matches Wu et al. 2021 exactly. Cells with fewer detected genes are
# typically empty droplets or extreme dropouts.

MIN_TOTAL_COUNTS_PER_CELL = 250
# Matches Wu et al. 2021 exactly. Wu required UMI count > 250 in addition
# to gene count > 200. 

MIN_CELLS_PER_GENE = 3
# Convention. Standard scanpy default; not specified by Wu et al. Genes
# detected in <3 cells out of ~100,000 are almost certainly noise.

MAX_PCT_MITO = 20.0
# Matches Wu et al. 2021 exactly. Wu's Methods state: "All cells with a
# mitochondrial UMI count percentage greater than 20% were removed."
# Independently confirmed by Yu et al. 2025 Nat Comms re-analysis of
# GSE176078 (doi:10.1038/s41467-025-58511-0).

MAX_TOTAL_COUNTS_MAD = 4
# MAD-based outlier filtering — NOT USED in current pipeline.
# Following Wu et al. 2021's exact QC (gene >200, UMI >250, mito <20%
# only), we do not apply MAD filtering. Doublets are handled by Scrublet.
# Variable retained here for reference only.
# MAX_TOTAL_COUNTS_MAD = 4

# Scrublet doublet detection — NOT used by Wu et al.
EXPECTED_DOUBLET_RATE = 0.06
# Our addition. Wu et al. used EmptyDrops (DropletUtils) which removes
# empty droplets but NOT doublets. Adding Scrublet here is a defensible
# improvement — Wu's deposited data may contain undetected doublets.
# 6% is the typical 10x Chromium prior at ~10,000-cell loading. Since
# Wu's data is already partially cleaned, true remaining doublet rate is
# probably lower (2-4%); Scrublet may slightly over-flag. Acceptable for
# our purposes.

# =============================================================================
# Normalization + HVG
# =============================================================================
NORMALIZE_TOTAL = 1e4
# Convention. scanpy/Seurat default. The exact value is arbitrary - what
# matters is that all cells are scaled to the same total before log1p.
# 10,000 places log-transformed values in a numerically convenient range.

N_HVGS = 2000
# Matches Wu et al. 2021 exactly. Wu's Methods: "A total of 2000 features
# for anchoring (FindIntegrationAnchors step)... were used." Note: Wu used
# 5000 features for sub-clustering immune/mesenchymal lineages; we keep
# 2000 for the main analysis, consistent with their primary integration.

# Genes we always keep regardless of HVG rank
GENE_OF_INTEREST = "ARC"
PROTECTED_GENES = [GENE_OF_INTEREST]

# =============================================================================
# Integration + clustering
# =============================================================================
N_PCS = 50
# Our choice. Wu et al. 2021 used 30 PCs for main integration and up to 100
# for non-batch-corrected clustering. 50 is a middle-ground default that
# captures most biological variance for ~10-20 cell types in breast tumor
# microenvironment. Could justify either 30 (Wu-strict) or 100 (more granular).

BATCH_KEY = "orig.ident"
# Wu et al.'s deposited metadata uses Seurat's default 'orig.ident' column
# to encode patient/sample ID (one ID per tumor donor in this dataset).
# Used as the batch variable for Harmony integration.

HARMONY_THETA = 2.0
# Default Harmony strength (higher = stronger batch correction).
# Wu et al. used Seurat's CCA-based integration (FindIntegrationAnchors +
# IntegrateData), not Harmony. We use Harmony as a faster CPU alternative
# that produces qualitatively similar results on this dataset.

N_NEIGHBORS = 15
# Scanpy default. Wu's Methods don't specify Seurat's k.param explicitly
# (they used FindNeighbors with defaults, which is k=20 in Seurat).
# 15 is the equivalent scanpy default. Effects on clustering are minor.

LEIDEN_RESOLUTIONS = [0.3, 0.5, 0.8, 1.0]
LEIDEN_DEFAULT = 0.8
# Wu et al. used resolution 0.8 (Seurat default at the time) with Louvain
# clustering. We test multiple resolutions and use 0.5 as a slightly more
# conservative default. Resolution 0.8 is included in the sweep so you can
# directly compare against Wu's clustering granularity if desired.

# =============================================================================
# Cell type marker genes
# =============================================================================
# Canonical lineage markers for independent validation of authors' annotations.
# Order matters for dot plots (groups related markers together).
CELL_TYPE_MARKERS = {
    "Epithelial_luminal":   ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "Epithelial_basal":     ["KRT5", "KRT14", "KRT17"],
    "T_cells":              ["CD3D", "CD3E", "CD8A", "CD4"],
    "B_cells":              ["MS4A1", "CD79A", "CD19"],
    "Plasma":               ["MZB1", "JCHAIN", "IGHG1"],
    "Myeloid":              ["CD68", "CD163", "LYZ", "CSF1R"],
    "Endothelial":          ["PECAM1", "VWF", "CDH5"],
    "Fibroblasts_CAFs":     ["COL1A1", "DCN", "PDGFRA"],
    "Smooth_muscle_like":   ["ACTA2", "MYH11"],
    "Perivascular":         ["RGS5", "MCAM"],
}

# =============================================================================
# ARC analysis
# =============================================================================
# Subtype labels as they appear in Wu et al. metadata.
# Clinical IHC-based subtype column: "subtype"
# PAM50 intrinsic molecular subtype column: "subtype_by_PAM50"
SUBTYPE_COLUMN = "subtype"               # Primary grouping variable
PAM50_COLUMN = "subtype_by_PAM50"        # Secondary grouping
CELLTYPE_COLUMN_AUTHORS = "celltype_major"  # Authors' major cell type labels
CELLTYPE_COLUMN_MINOR = "celltype_minor"    # Finer subdivisions

# Pseudo-bulk aggregation: only keep (patient, cell_type) groups with >= N cells
MIN_CELLS_FOR_PSEUDOBULK = 10
