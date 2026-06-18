"""
16S Analyzer — Longitudinal analysis: trajectories, volatility, temporal DA, heatmaps.
"""
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from skbio import DistanceMatrix

from app.analysis.alpha import compute_alpha
from app.analysis.beta import compute_distance
from app.analysis.r_runner import run_r_script
from app.analysis.shared import biom_to_count_df
from app.analysis.taxonomy import aggregate_taxonomy


def _parse_time_column(meta_df: pd.DataFrame, time_col: str) -> pd.Series:
    """Convert time column to numeric. Tries numeric first, then datetime."""
    vals = meta_df[time_col].copy()
    numeric = pd.to_numeric(vals, errors="coerce")
    if numeric.notna().sum() >= len(vals) * 0.5:
        return numeric
    # Try datetime conversion
    dt = pd.to_datetime(vals, errors="coerce")
    if dt.notna().sum() >= len(vals) * 0.5:
        min_dt = dt.min()
        return (dt - min_dt).dt.total_seconds() / 86400.0  # days
    return numeric


def compute_alpha_trajectories(
    count_df: pd.DataFrame,
    meta_df: pd.DataFrame,
    sample_id_col: str,
    subject_col: str,
    time_col: str,
    metric: str,
    matched_samples: list[str],
) -> pd.DataFrame:
    """Compute alpha diversity over time per subject.

    Returns DataFrame with columns: [subject, time, metric_value, sample_id].
    """
    count_sub = count_df[[s for s in count_df.columns if s in matched_samples]]
    diversity_df = compute_alpha(count_sub, [metric])

    # Merge with metadata
    meta = meta_df.copy()
    meta[sample_id_col] = meta[sample_id_col].astype(str)
    meta = meta[meta[sample_id_col].isin(matched_samples)]
    meta["_time_numeric"] = _parse_time_column(meta, time_col)

    result_rows = []
    for _, row in meta.iterrows():
        sid = str(row[sample_id_col])
        if sid in diversity_df.index:
            result_rows.append({
                "subject": row[subject_col],
                "time": row["_time_numeric"],
                "metric_value": diversity_df.loc[sid, metric],
                "sample_id": sid,
            })

    return pd.DataFrame(result_rows).dropna(subset=["time"])


def compute_volatility(
    count_df: pd.DataFrame,
    meta_df: pd.DataFrame,
    sample_id_col: str,
    subject_col: str,
    time_col: str,
    matched_samples: list[str],
    distance_metric: str = "braycurtis",
) -> pd.DataFrame:
    """Compute within-subject temporal beta diversity (volatility).

    Returns DataFrame with columns: [subject, time_from, time_to, distance].
    """
    count_sub = count_df[[s for s in count_df.columns if s in matched_samples]]
    dm = compute_distance(count_sub, distance_metric)
    dm_df = dm.to_data_frame()

    meta = meta_df.copy()
    meta[sample_id_col] = meta[sample_id_col].astype(str)
    meta = meta[meta[sample_id_col].isin(matched_samples)]
    meta["_time_numeric"] = _parse_time_column(meta, time_col)

    result_rows = []
    for subject in meta[subject_col].unique():
        subj_meta = meta[meta[subject_col] == subject].sort_values("_time_numeric")
        samples = subj_meta[sample_id_col].tolist()
        times = subj_meta["_time_numeric"].tolist()

        for i in range(len(samples) - 1):
            s1, s2 = samples[i], samples[i + 1]
            if s1 in dm_df.index and s2 in dm_df.columns:
                result_rows.append({
                    "subject": subject,
                    "time_from": times[i],
                    "time_to": times[i + 1],
                    "distance": dm_df.loc[s1, s2],
                })

    return pd.DataFrame(result_rows)


def compute_volatility_summary(volatility_df: pd.DataFrame) -> pd.DataFrame:
    """Compute mean volatility per subject from the volatility DataFrame."""
    if volatility_df.empty:
        return pd.DataFrame(columns=["subject", "mean_volatility"])
    return (
        volatility_df.groupby("subject")["distance"]
        .mean()
        .reset_index()
        .rename(columns={"distance": "mean_volatility"})
    )


def run_longitudinal_da(
    biom_path: str,
    meta_df: pd.DataFrame,
    sample_id_col: str,
    matched_samples: list[str],
    subject_col: str,
    time_col: str,
    group_col: str | None = None,
    level: str = "ASV",
    threads: int | None = None,
    on_log=None,
) -> pd.DataFrame:
    """Run MaAsLin2 with random effects for longitudinal differential abundance.

    Returns DataFrame with feature, coef, log2fc, pvalue, qvalue.
    """
    from app.analysis.taxonomy import aggregate_counts_by_level
    from app.config import DADA2_DEFAULTS

    count_df = biom_to_count_df(biom_path)
    if level != "ASV":
        count_df = aggregate_counts_by_level(biom_path, level)

    # Subset to matched samples
    keep = [s for s in count_df.columns if s in matched_samples]
    count_df = count_df[keep]
    count_df = count_df[count_df.sum(axis=1) > 0]

    # Prepare metadata
    meta_sub = meta_df[meta_df[sample_id_col].isin(keep)].copy()
    meta_sub["_time_numeric"] = _parse_time_column(meta_sub, time_col)

    # Write temp files
    tmp_dir = Path(tempfile.mkdtemp())
    counts_path = str(tmp_dir / "count_matrix.tsv")
    count_df.index.name = "feature_id"
    count_df.to_csv(counts_path, sep="\t")

    meta_path = str(tmp_dir / "metadata.tsv")
    meta_out_cols = [sample_id_col, subject_col, time_col]
    if group_col:
        meta_out_cols.append(group_col)
    meta_out = meta_sub[meta_out_cols].copy()
    # Use numeric time
    meta_out[time_col] = meta_sub["_time_numeric"]
    meta_out.columns = ["SampleID"] + list(meta_out.columns[1:])
    meta_out.to_csv(meta_path, sep="\t", index=False)

    output_path = str(tmp_dir / "results.tsv")

    n_threads = threads or DADA2_DEFAULTS.get("threads", 1)
    args = {
        "counts": counts_path,
        "metadata": meta_path,
        "subject_col": subject_col,
        "time_col": time_col,
        "output": output_path,
        "threads": n_threads,
    }
    if group_col:
        args["group_col"] = group_col

    result = run_r_script("run_maaslin2_longitudinal.R", args, on_line=on_log)

    if not result["success"]:
        raise RuntimeError(f"MaAsLin2 longitudinal failed: {result.get('error', 'unknown')}")

    if Path(output_path).exists():
        return pd.read_csv(output_path, sep="\t")
    return pd.DataFrame(columns=["feature", "coef", "log2fc", "pvalue", "qvalue"])


def compute_temporal_heatmap(
    biom_path: str,
    meta_df: pd.DataFrame,
    sample_id_col: str,
    subject_col: str,
    time_col: str,
    matched_samples: list[str],
    level: str = "Genus",
    top_n: int = 20,
) -> pd.DataFrame:
    """Compute mean relative abundance per timepoint for top N taxa.

    Returns DataFrame with rows=taxa, columns=timepoints (sorted).
    """
    tax_df = aggregate_taxonomy(biom_path, level, top_n=top_n)
    if tax_df.empty:
        return tax_df

    meta = meta_df.copy()
    meta[sample_id_col] = meta[sample_id_col].astype(str)
    meta = meta[meta[sample_id_col].isin(matched_samples)]
    meta["_time_numeric"] = _parse_time_column(meta, time_col)

    # Map sample -> timepoint
    sample_to_time = dict(zip(meta[sample_id_col], meta["_time_numeric"]))

    # Filter tax_df to matched samples
    common = [s for s in tax_df.columns if s in sample_to_time]
    if not common:
        return pd.DataFrame()

    tax_sub = tax_df[common]

    # Group by timepoint and compute mean
    timepoints = pd.Series([sample_to_time[s] for s in common], index=common)
    result = {}
    for tp in sorted(timepoints.unique()):
        if pd.isna(tp):
            continue
        tp_samples = timepoints[timepoints == tp].index.tolist()
        result[tp] = tax_sub[tp_samples].mean(axis=1)

    return pd.DataFrame(result)
