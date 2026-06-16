"""
Correlation analysis — taxon-taxon and taxon-metadata Spearman correlations.
"""
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def compute_taxon_correlation(count_df: pd.DataFrame, top_n: int = 30,
                              method: str = "spearman") -> dict:
    """Compute pairwise correlations between the top-N most abundant taxa.

    Parameters
    ----------
    count_df : DataFrame
        Features (rows) x samples (columns) count matrix.
    top_n : int
        Number of top taxa by mean abundance to include.
    method : str
        "spearman" or "pearson".

    Returns
    -------
    dict with keys: corr_matrix, pval_matrix, taxa_labels
    """
    # Select top taxa by mean relative abundance
    rel = count_df.div(count_df.sum(axis=0), axis=1)
    top_taxa = rel.mean(axis=1).nlargest(top_n).index.tolist()
    subset = count_df.loc[top_taxa].T  # samples x taxa

    n = len(top_taxa)
    corr_mat = np.ones((n, n))
    pval_mat = np.ones((n, n))

    for i in range(n):
        for j in range(i + 1, n):
            r, p = spearmanr(subset.iloc[:, i], subset.iloc[:, j])
            corr_mat[i, j] = corr_mat[j, i] = r
            pval_mat[i, j] = pval_mat[j, i] = p

    return {
        "corr_matrix": pd.DataFrame(corr_mat, index=top_taxa, columns=top_taxa),
        "pval_matrix": pd.DataFrame(pval_mat, index=top_taxa, columns=top_taxa),
        "taxa_labels": top_taxa,
    }


def compute_taxa_metadata_correlation(count_df: pd.DataFrame, meta_df: pd.DataFrame,
                                      sid_col: str, meta_cols: list[str],
                                      top_n: int = 30) -> dict:
    """Compute Spearman correlations between top taxa and numeric metadata variables.

    Returns
    -------
    dict with keys: corr_matrix, pval_matrix, taxa_labels, meta_labels
    """
    # Select top taxa
    rel = count_df.div(count_df.sum(axis=0), axis=1)
    top_taxa = rel.mean(axis=1).nlargest(top_n).index.tolist()
    taxa_df = count_df.loc[top_taxa].T  # samples x taxa

    # Align samples
    meta_indexed = meta_df.set_index(meta_df[sid_col].astype(str))
    common = sorted(set(taxa_df.index) & set(meta_indexed.index))
    if not common:
        return None

    taxa_aligned = taxa_df.loc[common]
    meta_aligned = meta_indexed.loc[common, meta_cols]

    # Convert to numeric, drop non-numeric columns
    meta_num = meta_aligned.apply(pd.to_numeric, errors="coerce")
    valid_cols = meta_num.columns[meta_num.notna().sum() >= 3].tolist()
    if not valid_cols:
        return None
    meta_num = meta_num[valid_cols]

    n_taxa = len(top_taxa)
    n_meta = len(valid_cols)
    corr_mat = np.zeros((n_taxa, n_meta))
    pval_mat = np.ones((n_taxa, n_meta))

    for i in range(n_taxa):
        for j in range(n_meta):
            mask = meta_num.iloc[:, j].notna()
            if mask.sum() < 3:
                continue
            r, p = spearmanr(taxa_aligned.iloc[:, i][mask], meta_num.iloc[:, j][mask])
            corr_mat[i, j] = r
            pval_mat[i, j] = p

    return {
        "corr_matrix": pd.DataFrame(corr_mat, index=top_taxa, columns=valid_cols),
        "pval_matrix": pd.DataFrame(pval_mat, index=top_taxa, columns=valid_cols),
        "taxa_labels": top_taxa,
        "meta_labels": valid_cols,
    }
