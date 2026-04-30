"""
Preprocessing utilities: normalization, HVG selection, PCA, Harmony integration,
neighborhood graph, UMAP, Leiden clustering.

Order of operations matters:
    1. Store raw counts in a layer ('counts')
    2. Normalize + log1p (operates on .X)
    3. HVG selection on normalized data (or raw with seurat_v3)
    4. Scale (optional; we keep unscaled for memory)
    5. PCA on HVGs
    6. Harmony on PCA embedding (batch correction)
    7. Neighbors + UMAP on Harmony embedding
    8. Leiden on neighbors graph
"""

from __future__ import annotations

from typing import Iterable, Optional

import numpy as np
import scanpy as sc
from anndata import AnnData

from src import config as cfg


def normalize_log(
    adata: AnnData,
    target_sum: float = 1e4,
    store_raw_counts_layer: bool = True,
) -> AnnData:
    """
    Normalize counts per cell and log1p-transform.

    The counts layer is preserved for HVG (seurat_v3 flavor needs raw counts)
    and for pseudo-bulk aggregation later.

    We do NOT use adata.raw to avoid doubling memory. Keep counts in a layer.
    """
    if store_raw_counts_layer:
        adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    return adata


def select_hvgs(
    adata: AnnData,
    n_top_genes: int = 2000,
    flavor: str = "seurat_v3",
    batch_key: Optional[str] = cfg.BATCH_KEY,
    protected_genes: Optional[Iterable[str]] = None,
) -> AnnData:
    """
    Select highly variable genes.

    flavor='seurat_v3' is the current best-practice default and requires
    raw counts. It correctly models the mean-variance relationship using
    a variance-stabilizing transformation.

    batch_key causes HVGs to be computed per batch then combined, which
    avoids picking genes that are variable only due to batch effects.

    protected_genes: genes forced into the HVG set regardless of rank
    (ARC in our case — it may not be variable enough to rank in top 2000).
    """
    if flavor == "seurat_v3":
        # seurat_v3 needs raw counts. Check layers first, fall back to .X
        layer = "counts" if "counts" in adata.layers else None
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_top_genes,
            flavor=flavor,
            batch_key=batch_key,
            layer=layer,
        )
    else:
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_top_genes,
            flavor=flavor,
            batch_key=batch_key,
        )

    # Force-include protected genes (e.g., ARC)
    if protected_genes is not None:
        protected = [g for g in protected_genes if g in adata.var_names]
        missing = [g for g in protected_genes if g not in adata.var_names]
        if missing:
            print(f"WARNING: protected genes not in dataset: {missing}")
        adata.var.loc[protected, "highly_variable"] = True
        print(f"Protected genes forced into HVG set: {protected}")

    print(f"HVGs selected: {adata.var['highly_variable'].sum()}")
    return adata


def run_pca(
    adata: AnnData,
    n_comps: int = 50,
    use_highly_variable: bool = True,
    random_state: int = 42,
) -> AnnData:
    """
    PCA on log-normalized data, restricted to HVGs.

    We do NOT call sc.pp.scale() before PCA. Scaling creates a dense matrix
    the size of (n_cells x n_hvgs) which on 130K cells x 2K genes is 2 GB.
    PCA without scaling still works well and saves memory.
    """
    sc.tl.pca(
        adata,
        n_comps=n_comps,
        use_highly_variable=use_highly_variable,
        random_state=random_state,
        svd_solver="arpack",
    )
    return adata


def run_harmony(
    adata: AnnData,
    batch_key: str = cfg.BATCH_KEY,
    theta: float = 2.0,
    random_state: int = 42,
) -> AnnData:
    """
    Harmony batch correction on the PCA embedding.

    Harmony iteratively corrects PCs so that cells from different batches
    mix while preserving biological structure. It operates on the embedding,
    not on raw expression, which is why it's fast and memory-efficient.

    Writes to obsm['X_pca_harmony']. Downstream steps use this embedding
    instead of 'X_pca'.
    """
    try:
        import harmonypy as hm
    except ImportError as err:
        raise ImportError(
            "harmonypy is required. Install with: pip install harmonypy"
        ) from err

    print(f"Running Harmony on '{batch_key}' (theta={theta})...")
    ho = hm.run_harmony(
        adata.obsm["X_pca"],
        adata.obs,
        vars_use=[batch_key],
        theta=theta,
        random_state=random_state,
        max_iter_harmony=20,
    )
    # Z_corr shape varies by harmonypy version.
    # Newer versions return (n_cells, n_PCs) directly — no transpose needed.
    # # Older versions return (n_PCs, n_cells) — needs .T
    Z = ho.Z_corr
    adata.obsm["X_pca_harmony"] = Z if Z.shape[0] == adata.n_obs else Z.T
    return adata


def run_neighbors_umap_leiden(
    adata: AnnData,
    use_rep: str = "X_pca_harmony",
    n_neighbors: int = 15,
    resolutions: Iterable[float] = (0.3, 0.5, 0.8, 1.0),
    random_state: int = 42,
) -> AnnData:
    """
    Build kNN graph, compute UMAP, and run Leiden at multiple resolutions.

    Multiple resolutions help decide granularity: you want enough clusters
    to resolve cell types but not so many that you fragment real populations.
    """
    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        use_rep=use_rep,
        random_state=random_state,
    )
    sc.tl.umap(adata, random_state=random_state)

    for res in resolutions:
        key = f"leiden_r{res}"
        sc.tl.leiden(
            adata,
            resolution=res,
            key_added=key,
            random_state=random_state,
            flavor="igraph",
            n_iterations=2,
            directed=False,
        )
        n_clusters = adata.obs[key].nunique()
        print(f"Leiden @ resolution {res}: {n_clusters} clusters")
    return adata


def stratified_subsample(
    adata: AnnData,
    group_key: str,
    n_per_group: int,
    random_state: int = 42,
) -> AnnData:
    """
    Subsample `n_per_group` cells from each group. Preserves group balance.

    Used for DEV_MODE execution on 8 GB Mac. Takes up to n_per_group cells
    per patient (fewer if a patient has <n_per_group cells).
    """
    rng = np.random.default_rng(random_state)
    keep_indices = []
    for group in adata.obs[group_key].unique():
        idx = np.where(adata.obs[group_key].values == group)[0]
        n_take = min(n_per_group, len(idx))
        chosen = rng.choice(idx, size=n_take, replace=False)
        keep_indices.append(chosen)
    keep_indices = np.concatenate(keep_indices)
    print(
        f"Stratified subsample on '{group_key}': "
        f"{adata.n_obs} -> {len(keep_indices)} cells "
        f"({adata.obs[group_key].nunique()} groups)"
    )
    return adata[keep_indices].copy()
