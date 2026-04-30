"""
Cell type annotation utilities.

We support two paths:
    (a) Use Wu et al.'s published annotations directly (`celltype_major`)
    (b) Score clusters against canonical marker lists, then assign types

Doing both gives cross-validation: if the authors' labels and independent
marker scoring agree, we're confident in the annotation.
"""

from __future__ import annotations

from typing import Dict, List

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData


def score_cell_type_markers(
    adata: AnnData,
    marker_dict: Dict[str, List[str]],
) -> AnnData:
    """
    Compute per-cell marker module scores for each cell type.

    Uses scanpy's `sc.tl.score_genes`, which calculates the average
    expression of a gene set minus a reference set of randomly selected
    genes with similar expression distribution. Stored in adata.obs
    as '{celltype}_score'.

    This is NOT a classifier — it produces a continuous score. Combine with
    clustering to assign cluster-level identities.
    """
    for cell_type, markers in marker_dict.items():
        present = [g for g in markers if g in adata.var_names]
        missing = [g for g in markers if g not in adata.var_names]
        if missing:
            print(f"[{cell_type}] missing markers: {missing}")
        if not present:
            print(f"[{cell_type}] SKIPPED — no markers in dataset")
            continue
        sc.tl.score_genes(
            adata,
            gene_list=present,
            score_name=f"{cell_type}_score",
            use_raw=False,
        )
    return adata


def assign_cluster_identities(
    adata: AnnData,
    marker_dict: Dict[str, List[str]],
    cluster_key: str = "leiden_r0.5",
) -> pd.DataFrame:
    """
    For each cluster, rank cell types by mean marker score.
    Return a DataFrame with the top assignment and the full score table.

    The top-scoring cell type per cluster is a first-pass annotation;
    always visualize the dotplot and adjust by hand for edge cases
    (e.g., proliferating T cells can score for T_cells + cycling markers).
    """
    score_cols = [f"{ct}_score" for ct in marker_dict if f"{ct}_score" in adata.obs.columns]
    if not score_cols:
        raise ValueError(
            "No marker scores found. Run score_cell_type_markers first."
        )

    # Mean score per cluster per cell type
    cluster_scores = adata.obs.groupby(cluster_key)[score_cols].mean()
    cluster_scores.columns = [c.replace("_score", "") for c in cluster_scores.columns]

    # Top assignment per cluster
    top_assignment = cluster_scores.idxmax(axis=1)
    top_score = cluster_scores.max(axis=1)
    result = pd.DataFrame({
        "top_cell_type": top_assignment,
        "top_score": top_score.round(3),
    })
    result = pd.concat([result, cluster_scores.round(3)], axis=1)
    return result


def apply_cluster_annotation(
    adata: AnnData,
    cluster_to_celltype: Dict[str, str],
    cluster_key: str = "leiden_r0.5",
    annotation_key: str = "celltype_predicted",
) -> AnnData:
    """
    Map cluster IDs to cell type labels, writing to adata.obs[annotation_key].

    cluster_to_celltype: e.g., {'0': 'T_cells', '1': 'Epithelial_luminal', ...}
    """
    adata.obs[annotation_key] = (
        adata.obs[cluster_key].astype(str).map(cluster_to_celltype)
    )
    # Fill unmapped clusters with 'Unknown' rather than NaN
    adata.obs[annotation_key] = adata.obs[annotation_key].fillna("Unknown")
    return adata


def compare_annotations(
    adata: AnnData,
    col_authors: str = "celltype_major",
    col_predicted: str = "celltype_predicted",
) -> pd.DataFrame:
    """
    Cross-tabulate authors' labels vs. our predicted labels.

    High diagonal values = good agreement. Off-diagonal entries point to
    cell types where our clustering disagrees with the paper.
    """
    if col_authors not in adata.obs.columns:
        raise ValueError(f"'{col_authors}' not in obs")
    if col_predicted not in adata.obs.columns:
        raise ValueError(f"'{col_predicted}' not in obs")
    ct = pd.crosstab(
        adata.obs[col_authors],
        adata.obs[col_predicted],
        margins=True,
    )
    return ct
