"""
Association analysis — CCA and dbRDA ordination biplots.

Links community composition to environmental/metadata variables using
constrained ordination methods.
"""
import numpy as np
import pandas as pd
from scipy.spatial.distance import braycurtis, pdist, squareform
from scipy.linalg import svd


def _center_matrix(X: np.ndarray) -> np.ndarray:
    """Column-center a matrix."""
    return X - X.mean(axis=0)


def _compute_dbRDA(dist_matrix: np.ndarray, env_matrix: np.ndarray) -> dict:
    """Distance-based Redundancy Analysis (db-RDA).

    1. PCoA on distance matrix to get sample coordinates
    2. RDA: constrained ordination of PCoA axes by environmental variables
    """
    n = dist_matrix.shape[0]

    # Step 1: PCoA (classical MDS)
    A = -0.5 * dist_matrix ** 2
    centering = np.eye(n) - np.ones((n, n)) / n
    G = centering @ A @ centering

    eigvals, eigvecs = np.linalg.eigh(G)
    idx = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]

    # Keep positive eigenvalues
    pos = eigvals > 1e-10
    eigvals_pos = eigvals[pos]
    coords = eigvecs[:, pos] * np.sqrt(eigvals_pos)

    # Step 2: RDA on PCoA coordinates
    Y = _center_matrix(coords)
    X = _center_matrix(env_matrix)

    # Fitted values: Y_hat = X(X'X)^-1 X'Y
    try:
        XtX_inv = np.linalg.pinv(X.T @ X)
        H = X @ XtX_inv @ X.T
        Y_hat = H @ Y
    except np.linalg.LinAlgError:
        Y_hat = Y

    # SVD of fitted values
    U, s, Vt = svd(Y_hat, full_matrices=False)

    # Constrained axes (take first 2)
    n_axes = min(2, len(s))
    sample_scores = U[:, :n_axes] * s[:n_axes]

    # Biplot scores for environmental variables
    # Project env variables onto constrained axes
    env_scores = X.T @ U[:, :n_axes] / np.sqrt(n - 1)

    # Proportion explained
    total_var = np.sum(eigvals_pos)
    constrained_var = np.sum(s[:n_axes] ** 2)
    prop_explained = s[:n_axes] ** 2 / total_var if total_var > 0 else s[:n_axes] * 0

    return {
        "sample_scores": sample_scores,
        "env_scores": env_scores,
        "prop_explained": prop_explained,
        "n_axes": n_axes,
    }


def _compute_CCA(count_matrix: np.ndarray, env_matrix: np.ndarray) -> dict:
    """Canonical Correspondence Analysis (CCA).

    Chi-square-based constrained ordination for community composition data.
    """
    # Chi-square transformation
    Y = count_matrix.astype(float)
    grand_total = Y.sum()
    if grand_total == 0:
        raise ValueError("Count matrix is all zeros")

    row_sums = Y.sum(axis=1, keepdims=True)
    col_sums = Y.sum(axis=0, keepdims=True)

    # Avoid division by zero
    row_sums = np.where(row_sums == 0, 1, row_sums)
    col_sums = np.where(col_sums == 0, 1, col_sums)

    row_weights = np.sqrt(row_sums.flatten() / grand_total)
    col_weights = np.sqrt(col_sums.flatten() / grand_total)

    # Chi-square standardized residuals
    expected = row_sums * col_sums / grand_total
    Q = (Y - expected) / np.sqrt(expected)

    # Weight by row masses
    Q_weighted = Q / row_weights[:, np.newaxis]

    # Center and constrain by environmental variables
    X = _center_matrix(env_matrix)
    # Weight X by row masses
    W = np.diag(row_weights ** 2)
    try:
        XtWX_inv = np.linalg.pinv(X.T @ W @ X)
        H = X @ XtWX_inv @ X.T @ W
    except np.linalg.LinAlgError:
        H = np.eye(Q_weighted.shape[0])

    Q_hat = H @ Q_weighted

    # SVD
    U, s, Vt = svd(Q_hat, full_matrices=False)
    n_axes = min(2, len(s))

    sample_scores = U[:, :n_axes] * s[:n_axes]

    # Environmental variable scores
    n = env_matrix.shape[0]
    env_scores = X.T @ W @ U[:, :n_axes] / np.sqrt(n - 1)

    # Species scores
    species_scores = Vt[:n_axes, :].T * s[:n_axes]

    # Proportion explained
    total_inertia = np.sum(Q ** 2) / grand_total
    prop_explained = s[:n_axes] ** 2 / total_inertia if total_inertia > 0 else s[:n_axes] * 0

    return {
        "sample_scores": sample_scores,
        "env_scores": env_scores,
        "species_scores": species_scores,
        "prop_explained": prop_explained,
        "n_axes": n_axes,
    }


def run_association(count_df: pd.DataFrame, meta_df: pd.DataFrame,
                    sid_col: str, env_cols: list[str],
                    group_col: str | None = None,
                    method: str = "dbRDA",
                    top_n: int = 50) -> dict:
    """Run constrained ordination analysis.

    Parameters
    ----------
    count_df : DataFrame
        Features (rows) x samples (columns).
    meta_df : DataFrame
        Metadata with sample IDs and environmental variables.
    sid_col : str
        Column name for sample IDs in metadata.
    env_cols : list[str]
        Columns to use as environmental/constraining variables.
    group_col : str or None
        Column for coloring samples in the biplot.
    method : str
        "dbRDA" or "CCA".
    top_n : int
        Number of top taxa to include.

    Returns
    -------
    dict with sample_scores, env_scores, env_labels, sample_ids, groups,
         prop_explained, species_labels (CCA only), species_scores (CCA only)
    """
    # Select top taxa
    rel = count_df.div(count_df.sum(axis=0), axis=1)
    top_taxa = rel.mean(axis=1).nlargest(top_n).index.tolist()
    subset = count_df.loc[top_taxa]

    # Align samples
    meta_indexed = meta_df.set_index(meta_df[sid_col].astype(str))
    common = sorted(set(subset.columns) & set(meta_indexed.index))
    if len(common) < 4:
        raise ValueError(f"Need at least 4 matching samples, found {len(common)}")

    count_aligned = subset[common].T.values  # samples x taxa
    env_data = meta_indexed.loc[common, env_cols].apply(pd.to_numeric, errors="coerce")

    # Drop columns/rows with NaN
    valid_cols = env_data.columns[env_data.notna().all()].tolist()
    if not valid_cols:
        raise ValueError("No numeric environmental variables without missing values")
    env_data = env_data[valid_cols].values

    groups = None
    if group_col:
        groups = meta_indexed.loc[common, group_col].tolist()

    if method == "CCA":
        result = _compute_CCA(count_aligned, env_data)
        # Include top species scores
        top_sp = min(10, len(top_taxa))
        sp_norms = np.linalg.norm(result["species_scores"], axis=1)
        top_sp_idx = np.argsort(sp_norms)[-top_sp:]
        result["species_labels"] = [top_taxa[i] for i in top_sp_idx]
        result["species_scores"] = result["species_scores"][top_sp_idx]
    else:
        # dbRDA: compute Bray-Curtis distance matrix
        dist_vec = pdist(count_aligned, metric="braycurtis")
        dist_mat = squareform(dist_vec)
        result = _compute_dbRDA(dist_mat, env_data)

    result["sample_ids"] = common
    result["groups"] = groups
    result["env_labels"] = valid_cols

    return result
