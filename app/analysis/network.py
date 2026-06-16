"""
Network analysis — SparCC-like correlation network inference.

Computes a sparse correlation network from compositional count data using
iterative log-ratio variance estimation (SparCC algorithm).
"""
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def _sparcc_simple(count_df: pd.DataFrame, n_iter: int = 20,
                   threshold: float = 0.1) -> pd.DataFrame:
    """Simplified SparCC-like correlation estimation.

    Uses iterative exclusion of highly correlated pairs from log-ratio
    variance estimation.  This is a simplified version suitable for
    interactive use — for publication-quality results, users should
    consider the full SparCC implementation.

    Parameters
    ----------
    count_df : DataFrame
        Taxa (rows) x samples (columns), raw counts.
    n_iter : int
        Number of exclusion iterations.
    threshold : float
        Correlation magnitude below which a pair is excluded from refinement.

    Returns
    -------
    DataFrame (taxa x taxa) estimated correlations.
    """
    # Add pseudocount, convert to relative abundance, then log
    data = count_df.values.astype(float) + 1
    fracs = data / data.sum(axis=0, keepdims=True)
    log_fracs = np.log(fracs)

    p = log_fracs.shape[0]  # number of taxa

    # Compute basis variance-covariance from log-ratios
    # Var(log(xi/xj)) = var_i + var_j - 2*cov_ij
    t_data = log_fracs.T  # samples x taxa
    var_log = np.var(t_data, axis=0)

    # Initialize with sample correlation
    corr, _ = spearmanr(t_data)
    if p == 2:
        corr = np.array([[1.0, corr], [corr, 1.0]])

    # Iterative refinement
    for _ in range(n_iter):
        # Estimate basis variances using median of log-ratio variances
        basis_var = np.zeros(p)
        for i in range(p):
            ratios = []
            for j in range(p):
                if i == j:
                    continue
                lr = t_data[:, i] - t_data[:, j]
                ratios.append(np.var(lr))
            basis_var[i] = np.median(ratios)

        # Re-estimate correlations
        new_corr = np.eye(p)
        for i in range(p):
            for j in range(i + 1, p):
                lr_var = np.var(t_data[:, i] - t_data[:, j])
                cov_ij = 0.5 * (basis_var[i] + basis_var[j] - lr_var)
                denom = np.sqrt(basis_var[i] * basis_var[j])
                if denom > 0:
                    r = cov_ij / denom
                    r = np.clip(r, -1, 1)
                else:
                    r = 0.0
                new_corr[i, j] = new_corr[j, i] = r

        corr = new_corr

    taxa = count_df.index.tolist()
    return pd.DataFrame(corr, index=taxa, columns=taxa)


def build_network(count_df: pd.DataFrame, top_n: int = 50,
                  corr_threshold: float = 0.3,
                  method: str = "sparcc") -> dict:
    """Build a co-occurrence network from count data.

    Parameters
    ----------
    count_df : DataFrame
        Features (rows) x samples (columns).
    top_n : int
        Number of top taxa by mean abundance.
    corr_threshold : float
        Minimum absolute correlation to include an edge.
    method : str
        "sparcc" or "spearman".

    Returns
    -------
    dict with keys: nodes (list of dicts), edges (list of dicts), corr_matrix
    """
    # Select top taxa
    rel = count_df.div(count_df.sum(axis=0), axis=1)
    top_taxa = rel.mean(axis=1).nlargest(top_n).index.tolist()
    subset = count_df.loc[top_taxa]

    if method == "sparcc":
        corr_df = _sparcc_simple(subset)
    else:
        corr_mat, _ = spearmanr(subset.values.T)
        if len(top_taxa) == 2:
            corr_mat = np.array([[1.0, corr_mat], [corr_mat, 1.0]])
        corr_df = pd.DataFrame(corr_mat, index=top_taxa, columns=top_taxa)

    # Build nodes and edges
    taxa = corr_df.index.tolist()
    mean_abund = rel.loc[top_taxa].mean(axis=1)

    nodes = []
    for t in taxa:
        nodes.append({
            "id": t,
            "abundance": float(mean_abund[t]),
        })

    edges = []
    for i in range(len(taxa)):
        for j in range(i + 1, len(taxa)):
            r = corr_df.iloc[i, j]
            if abs(r) >= corr_threshold:
                edges.append({
                    "source": taxa[i],
                    "target": taxa[j],
                    "weight": float(r),
                    "positive": bool(r > 0),
                })

    return {
        "nodes": nodes,
        "edges": edges,
        "corr_matrix": corr_df,
    }
