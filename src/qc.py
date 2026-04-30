"""
Quality control utilities for scRNA-seq data.

Functions:
    annotate_qc_metrics   - Add mitochondrial/ribosomal % and total-count metrics
    filter_cells_genes    - Hard thresholds on min genes / min cells
    filter_by_mito        - Remove cells above mitochondrial % cutoff
    filter_by_counts_mad  - Remove cells with outlier total counts (MAD-based)
    run_scrublet          - Per-sample doublet detection
    summarize_qc          - Print a compact QC report
"""

from __future__ import annotations

import warnings
from typing import Optional

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from src import config as cfg


def annotate_qc_metrics(adata: AnnData) -> None:
    """
    Add per-cell QC metrics to adata.obs in-place.

    Adds:
        pct_counts_mt      - percent of counts from mitochondrial genes (MT-*)
        pct_counts_ribo    - percent of counts from ribosomal genes (RPS*/RPL*)
        total_counts, n_genes_by_counts - total UMI and number of detected genes

    Note: this uses gene SYMBOLS. Wu et al. data uses symbols like "MT-ND1".
    """
    # Flag mitochondrial and ribosomal genes in adata.var
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))

    # scanpy computes total_counts, n_genes_by_counts, pct_counts_mt, pct_counts_ribo
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )


def filter_cells_genes(
    adata: AnnData,
    min_genes: int = 200,
    min_cells: int = 3,
) -> AnnData:
    """
    Apply hard QC thresholds.

    Parameters
    ----------
    min_genes : int
        Minimum genes a cell must express to be kept.
    min_cells : int
        Minimum cells a gene must be detected in to be kept.

    Returns a new AnnData (not in place). Important: scanpy's filter_* functions
    modify in-place, but returning a view allows chaining in notebooks.
    """
    n_cells_before, n_genes_before = adata.shape
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    print(
        f"filter_cells_genes: "
        f"{n_cells_before} -> {adata.n_obs} cells, "
        f"{n_genes_before} -> {adata.n_vars} genes"
    )
    return adata


def filter_by_mito(adata: AnnData, max_pct: float = 20.0) -> AnnData:
    """
    Remove cells above the mitochondrial percentage threshold.

    Why 20% for tumors: dying/hypoxic tumor cells have elevated mito content
    but are biologically meaningful. Wu et al. 2021 used 20% for this dataset.
    For healthy tissue, 10% is more typical.
    """
    if "pct_counts_mt" not in adata.obs.columns:
        raise ValueError(
            "pct_counts_mt not found. Run annotate_qc_metrics first."
        )
    n_before = adata.n_obs
    adata = adata[adata.obs["pct_counts_mt"] < max_pct].copy()
    print(f"filter_by_mito (<{max_pct}%): {n_before} -> {adata.n_obs} cells")
    return adata


def filter_by_counts_mad(
    adata: AnnData,
    batch_key: Optional[str] = None,
    n_mads: float = 4.0,
) -> AnnData:
    """
    Remove cells with total-count values above n_mads MADs from the median.

    Applied per batch (usually patient) if batch_key is given, because
    different samples have different sequencing depths and we don't want
    to penalize deeper-sequenced batches.

    This flags likely doublets and over-captured cells. Not a substitute
    for Scrublet, but a cheap complementary filter.
    """
    if "total_counts" not in adata.obs.columns:
        raise ValueError(
            "total_counts not found. Run annotate_qc_metrics first."
        )

    def _is_outlier(values: np.ndarray) -> np.ndarray:
        log_values = np.log1p(values)  # Use log-scale: counts are right-skewed
        median = np.median(log_values)
        mad = np.median(np.abs(log_values - median))
        if mad == 0:
            return np.zeros_like(values, dtype=bool)
        return log_values > median + n_mads * mad

    n_before = adata.n_obs
    if batch_key and batch_key in adata.obs.columns:
        outlier_mask = np.zeros(adata.n_obs, dtype=bool)
        for batch in adata.obs[batch_key].unique():
            idx = adata.obs[batch_key] == batch
            outlier_mask[idx] = _is_outlier(
                adata.obs.loc[idx, "total_counts"].values
            )
    else:
        outlier_mask = _is_outlier(adata.obs["total_counts"].values)

    adata = adata[~outlier_mask].copy()
    print(
        f"filter_by_counts_mad ({n_mads} MADs): "
        f"{n_before} -> {adata.n_obs} cells"
    )
    return adata


def run_scrublet(
    adata: AnnData,
    batch_key: str = cfg.BATCH_KEY,
    expected_doublet_rate: float = 0.06,
    random_state: int = 42,
    fallback_threshold: float = 0.25,
    verbose: bool = True,
) -> AnnData:
    """
    Run Scrublet per-batch to detect doublets, then filter them out.

    Per-sample is critical: Scrublet's simulated doublets assume cells
    come from the same droplet population. Pooling samples first creates
    spurious "doublets" from cross-sample expression differences.

    fallback_threshold: when Scrublet cannot auto-detect a threshold from
    the score distribution (common on pre-filtered data with few obvious
    doublets), use this manual value. 0.25 is a defensible default that
    flags clear outliers without over-filtering.
    """
    try:
        import scrublet as scr
    except ImportError as err:
        raise ImportError(
            "scrublet is required for doublet detection. "
            "Install with: pip install scrublet"
        ) from err

    doublet_scores = np.zeros(adata.n_obs)
    predicted_doublets = np.zeros(adata.n_obs, dtype=bool)
    per_patient_summary = []

    for batch in adata.obs[batch_key].unique():
        mask = (adata.obs[batch_key] == batch).values
        n_cells = int(mask.sum())
        if n_cells < 30:
            warnings.warn(f"Skipping Scrublet for '{batch}': only {n_cells} cells")
            per_patient_summary.append((batch, n_cells, 0, "skipped: <30 cells"))
            continue

        # Use AnnData slicing then access .X — safer than direct matrix subset
        sub_adata = adata[mask].copy()
        sub_counts = sub_adata.X

        try:
            scrub = scr.Scrublet(
                sub_counts,
                expected_doublet_rate=expected_doublet_rate,
                random_state=random_state,
            )
            scores, preds = scrub.scrub_doublets(
                min_counts=2,
                min_cells=3,
                min_gene_variability_pctl=85,
                n_prin_comps=30,
                verbose=False,
            )
        except Exception as err:
            warnings.warn(f"Scrublet failed for '{batch}': {err}")
            per_patient_summary.append((batch, n_cells, 0, f"failed: {err}"))
            continue

        if scores is None:
            per_patient_summary.append((batch, n_cells, 0, "no scores returned"))
            continue

        doublet_scores[mask] = scores

        # Fall back to manual threshold if Scrublet couldn't auto-detect
        if preds is None or preds.sum() == 0:
            preds = scores >= fallback_threshold
            note = f"auto-threshold failed; used fallback >={fallback_threshold}"
        else:
            note = "auto-threshold OK"

        predicted_doublets[mask] = preds
        n_doublets = int(preds.sum())
        per_patient_summary.append((batch, n_cells, n_doublets, note))

    adata.obs["doublet_score"] = doublet_scores
    adata.obs["predicted_doublet"] = predicted_doublets

    if verbose:
        print(f"\n{'Patient':<12} {'n_cells':>8} {'doublets':>10} {'%':>6}  Note")
        print("-" * 80)
        for batch, n_cells, n_db, note in per_patient_summary:
            pct = 100 * n_db / n_cells if n_cells else 0
            print(f"{str(batch):<12} {n_cells:>8} {n_db:>10} {pct:>5.1f}%  {note}")

    n_before = adata.n_obs
    adata = adata[~adata.obs["predicted_doublet"]].copy()
    n_removed = n_before - adata.n_obs
    print(
        f"\nrun_scrublet: removed {n_removed} doublets "
        f"({100*n_removed/n_before:.2f}%)"
    )
    return adata


def summarize_qc(adata: AnnData) -> pd.DataFrame:
    """
    Print and return a per-batch QC summary table.
    Useful for spotting outlier patients before integration.
    """
    df = adata.obs.groupby(cfg.BATCH_KEY).agg(
        n_cells=(cfg.BATCH_KEY, "size"),
        median_genes=("n_genes_by_counts", "median"),
        median_counts=("total_counts", "median"),
        median_pct_mito=("pct_counts_mt", "median"),
    ).round(2)
    print(df)
    return df