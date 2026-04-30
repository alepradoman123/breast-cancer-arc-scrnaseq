"""
ARC expression analysis across breast cancer subtypes and cell types.

Key analyses:
    1. Per-cell ARC expression summary (mean, median, % expressing)
       stratified by (subtype, cell_type).
    2. Pseudo-bulk aggregation: sum raw counts per (patient, cell_type)
       and compare across subtypes. This is the statistically valid approach
       because cells from the same patient are not independent.
    3. Rank cell types by ARC expression within each subtype.

Why pseudo-bulk matters:
    Per-cell Wilcoxon tests on scRNA-seq DE analysis dramatically inflate
    p-values because they ignore the patient-level structure of the data
    (Squair et al. 2021, Nature Communications). Aggregating to pseudo-bulk
    gives one observation per patient-celltype, which is the correct unit
    of replication for subtype comparisons.
"""

from __future__ import annotations

from typing import Optional, Tuple

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

from src import config as cfg


def _get_gene_vector(adata: AnnData, gene: str, layer: Optional[str] = None) -> np.ndarray:
    """Extract a dense 1D expression vector for one gene."""
    if gene not in adata.var_names:
        raise ValueError(f"Gene '{gene}' not in dataset")
    idx = adata.var_names.get_loc(gene)
    if layer is not None:
        mat = adata.layers[layer]
    else:
        mat = adata.X
    col = mat[:, idx]
    if sparse.issparse(col):
        col = col.toarray().ravel()
    else:
        col = np.asarray(col).ravel()
    return col


def summarize_gene_by_group(
    adata: AnnData,
    gene: str = "ARC",
    group_cols: Tuple[str, ...] = ("subtype", "celltype_major"),
    layer: Optional[str] = None,
) -> pd.DataFrame:
    """
    Summarize expression of one gene per group.

    Reports:
        n_cells        - total cells in group
        n_expressing   - cells with non-zero expression
        pct_expressing - percent of cells expressing
        mean_expr      - mean expression (log-normalized if layer is None)
        median_expr    - median expression
        mean_expr_pos  - mean across expressing cells only (avoids zero dilution)
    """
    expr = _get_gene_vector(adata, gene, layer=layer)
    df = adata.obs[list(group_cols)].copy()
    df["expr"] = expr
    df["expressing"] = expr > 0

    def agg(x: pd.DataFrame) -> pd.Series:
        pos = x["expr"][x["expressing"]]
        return pd.Series({
            "n_cells": len(x),
            "n_expressing": int(x["expressing"].sum()),
            "pct_expressing": 100 * x["expressing"].mean(),
            "mean_expr": x["expr"].mean(),
            "median_expr": x["expr"].median(),
            "mean_expr_pos": pos.mean() if len(pos) > 0 else np.nan,
        })

    summary = df.groupby(list(group_cols), observed=True).apply(agg).reset_index()
    numeric_cols = ["pct_expressing", "mean_expr", "median_expr", "mean_expr_pos"]
    summary[numeric_cols] = summary[numeric_cols].round(3)
    return summary


def top_cell_types_per_subtype(
    summary: pd.DataFrame,
    subtype_col: str = "subtype",
    celltype_col: str = "celltype_major",
    rank_by: str = "mean_expr",
    top_n: int = 3,
    min_cells: int = 20,
) -> pd.DataFrame:
    """
    For each subtype, return the top N cell types by `rank_by`.

    Filters out (subtype, cell_type) combinations with fewer than min_cells
    to avoid rankings driven by tiny groups.
    """
    df = summary[summary["n_cells"] >= min_cells].copy()
    df = df.sort_values([subtype_col, rank_by], ascending=[True, False])
    top = df.groupby(subtype_col, observed=True).head(top_n)
    return top.reset_index(drop=True)


def pseudobulk_counts(
    adata: AnnData,
    group_cols: Tuple[str, ...] = (cfg.BATCH_KEY, "celltype_major"),
    counts_layer: str = "counts",
    min_cells: int = 10,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Aggregate raw counts to pseudo-bulk per group.

    Returns:
        counts_df  - (n_groups, n_genes) summed raw counts
        metadata   - per-group metadata (n_cells, and all group_cols)

    Only groups with >= min_cells are kept. This is critical — pseudo-bulk
    from a handful of cells is noisy.
    """
    if counts_layer not in adata.layers:
        raise ValueError(
            f"Layer '{counts_layer}' not found. Store raw counts there first."
        )

    obs = adata.obs[list(group_cols)].copy()
    obs["_group"] = obs.apply(
        lambda r: "__".join(str(r[c]) for c in group_cols), axis=1
    )
    group_sizes = obs["_group"].value_counts()
    keep_groups = group_sizes[group_sizes >= min_cells].index.tolist()

    X = adata.layers[counts_layer]
    records = []
    meta_records = []
    for g in keep_groups:
        mask = (obs["_group"] == g).values
        sub = X[mask]
        if sparse.issparse(sub):
            summed = np.asarray(sub.sum(axis=0)).ravel()
        else:
            summed = sub.sum(axis=0)
        records.append(summed)

        group_vals = obs.loc[mask, list(group_cols)].iloc[0].to_dict()
        group_vals["n_cells"] = int(mask.sum())
        group_vals["_group"] = g
        meta_records.append(group_vals)

    counts_df = pd.DataFrame(
        records, index=[m["_group"] for m in meta_records], columns=adata.var_names
    )
    metadata = pd.DataFrame(meta_records).set_index("_group")
    return counts_df, metadata


def pseudobulk_expression(
    counts_df: pd.DataFrame,
    metadata: pd.DataFrame,
    gene: str = "ARC",
    normalize: bool = True,
) -> pd.DataFrame:
    """
    Extract one gene's pseudo-bulk expression. If normalize=True, returns
    counts per million (CPM) on log1p scale — suitable for cross-sample
    comparisons.
    """
    if gene not in counts_df.columns:
        raise ValueError(f"Gene '{gene}' not in counts matrix")

    if normalize:
        totals = counts_df.sum(axis=1)
        cpm = counts_df.div(totals, axis=0) * 1e6
        expr = np.log1p(cpm[gene])
    else:
        expr = counts_df[gene]

    out = metadata.copy()
    out[f"{gene}_expr"] = expr.values
    return out


def wilcoxon_per_celltype(
    pseudobulk_df: pd.DataFrame,
    gene: str = "ARC",
    subtype_col: str = "subtype",
    celltype_col: str = "celltype_major",
    min_samples: int = 3,
) -> pd.DataFrame:
    """
    Pairwise Wilcoxon rank-sum tests between subtypes WITHIN each cell type,
    using pseudo-bulk values (one per patient-celltype).

    With typically 3-10 patients per subtype, this is correctly powered
    for pseudo-bulk but very noisy — treat as exploratory. Report effect
    sizes (median difference) rather than relying on p-values alone.
    """
    from scipy.stats import mannwhitneyu
    from itertools import combinations

    expr_col = f"{gene}_expr"
    rows = []
    for ct in pseudobulk_df[celltype_col].unique():
        sub = pseudobulk_df[pseudobulk_df[celltype_col] == ct]
        subtypes = sub[subtype_col].unique()
        for s1, s2 in combinations(subtypes, 2):
            x = sub.loc[sub[subtype_col] == s1, expr_col].values
            y = sub.loc[sub[subtype_col] == s2, expr_col].values
            if len(x) < min_samples or len(y) < min_samples:
                continue
            try:
                stat, p = mannwhitneyu(x, y, alternative="two-sided")
            except ValueError:
                # All values identical -> test undefined
                stat, p = np.nan, np.nan
            rows.append({
                celltype_col: ct,
                "subtype_1": s1,
                "subtype_2": s2,
                "n1": len(x), "n2": len(y),
                "median_1": np.median(x),
                "median_2": np.median(y),
                "median_diff": np.median(x) - np.median(y),
                "u_stat": stat,
                "p_value": p,
            })
    df = pd.DataFrame(rows)
    if not df.empty:
        # FDR correction across all comparisons
        from statsmodels.stats.multitest import multipletests
        valid = df["p_value"].notna()
        df["fdr"] = np.nan
        if valid.any():
            df.loc[valid, "fdr"] = multipletests(
                df.loc[valid, "p_value"], method="fdr_bh"
            )[1]
    return df.sort_values("p_value")
