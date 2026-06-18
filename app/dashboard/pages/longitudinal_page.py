"""
16S Analyzer — Longitudinal Analysis page.
"""
import io
import traceback

import dash_bootstrap_components as dbc
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from biom import load_table
from dash import Input, Output, State, callback_context, dcc, html, no_update

from app.analysis.longitudinal import (
    compute_alpha_trajectories,
    compute_temporal_heatmap,
    compute_volatility,
    compute_volatility_summary,
    run_longitudinal_da,
)
from app.analysis.shared import (
    biom_to_count_df,
    find_metadata_for_samples,
    get_dataset_metadata_df,
    get_group_columns,
    get_pipeline_biom_options,
    parse_uploaded_biom,
    parse_uploaded_metadata,
    validate_metadata_vs_biom,
)
from app.dashboard.app import app as dash_app

ANALYSIS_OPTIONS = [
    {"label": "Alpha Trajectories", "value": "alpha_traj"},
    {"label": "Beta Volatility", "value": "volatility"},
    {"label": "Temporal Diff. Abundance", "value": "temporal_da"},
    {"label": "Temporal Heatmap", "value": "heatmap"},
]

ALPHA_METRICS = [
    {"label": "Shannon", "value": "shannon"},
    {"label": "Simpson", "value": "simpson"},
    {"label": "Observed OTUs", "value": "observed_otus"},
    {"label": "Chao1", "value": "chao1"},
]

DISTANCE_METRICS = [
    {"label": "Bray-Curtis", "value": "braycurtis"},
    {"label": "Jaccard", "value": "jaccard"},
]

TAXONOMY_LEVELS = [
    {"label": "ASV", "value": "ASV"},
    {"label": "Phylum", "value": "Phylum"},
    {"label": "Class", "value": "Class"},
    {"label": "Order", "value": "Order"},
    {"label": "Family", "value": "Family"},
    {"label": "Genus", "value": "Genus"},
    {"label": "Species", "value": "Species"},
]

COLORS = ["#3498db", "#e74c3c", "#2ecc71", "#f39c12", "#9b59b6",
           "#1abc9c", "#e67e22", "#34495e", "#16a085", "#c0392b",
           "#8e44ad", "#2980b9", "#d35400", "#27ae60"]


def get_layout():
    pipeline_opts = get_pipeline_biom_options()
    return dbc.Container([
        html.H3("Longitudinal Analysis", className="mb-2"),
        html.P(
            "Analyze microbiome changes over time with repeated measures.",
            className="text-muted mb-4",
        ),
        dcc.Store(id="lo-biom-path"),
        dcc.Store(id="lo-meta-store"),
        dcc.Store(id="lo-sample-id-col"),
        dcc.Store(id="lo-result-csv"),

        dbc.Row([
            # Left column: inputs
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader("Input Data"),
                    dbc.CardBody([
                        dbc.Label("BIOM Table", className="fw-bold"),
                        dbc.Select(
                            id="lo-select-pipeline",
                            options=[{"label": "— none —", "value": ""}]
                            + pipeline_opts,
                            value="",
                            className="mb-2",
                        ),
                        html.Div("— or upload —", className="text-center text-muted small mb-2"),
                        dcc.Upload(
                            id="lo-upload-biom",
                            children=html.Div(["Drag & drop or ", html.A("select .biom")]),
                            style={"borderWidth": "2px", "borderStyle": "dashed",
                                   "borderRadius": "5px", "borderColor": "#555",
                                   "textAlign": "center", "padding": "8px"},
                            multiple=False,
                        ),
                        html.Div(id="lo-biom-status", className="mt-1 small mb-3"),

                        dbc.Label("Metadata (CSV/TSV)", className="fw-bold"),
                        html.Div(
                            "Must contain subject ID and time columns.",
                            className="text-muted small mb-1",
                        ),
                        dcc.Upload(
                            id="lo-upload-meta",
                            children=html.Div(["Drag & drop or ", html.A("select metadata file")]),
                            style={"borderWidth": "2px", "borderStyle": "dashed",
                                   "borderRadius": "5px", "borderColor": "#555",
                                   "textAlign": "center", "padding": "8px"},
                            multiple=False,
                        ),
                        html.Div(id="lo-meta-status", className="mt-1 small mb-3"),
                    ]),
                ], className="mb-3"),

                dbc.Card([
                    dbc.CardHeader("Longitudinal Settings"),
                    dbc.CardBody([
                        dbc.Label("Subject ID Column", className="fw-bold"),
                        dbc.Select(id="lo-subject-col", placeholder="Select...",
                                   className="mb-2"),

                        dbc.Label("Time Column", className="fw-bold"),
                        dbc.Select(id="lo-time-col", placeholder="Select...",
                                   className="mb-2"),

                        dbc.Label("Group Column (optional)", className="fw-bold"),
                        dbc.Select(id="lo-group-col", placeholder="None (no grouping)",
                                   className="mb-2"),

                        html.Hr(),

                        dbc.Label("Analysis Type", className="fw-bold"),
                        dbc.RadioItems(
                            id="lo-analysis-type",
                            options=ANALYSIS_OPTIONS,
                            value="alpha_traj",
                            className="mb-3",
                        ),

                        dbc.Label("Metric / Level", className="fw-bold"),
                        dbc.Select(
                            id="lo-metric",
                            options=ALPHA_METRICS,
                            value="shannon",
                            className="mb-3",
                        ),

                        html.Div([
                            dbc.Label("Top N taxa", className="fw-bold"),
                            dbc.Input(id="lo-top-n", type="number", value=20, min=5, max=50,
                                      className="mb-3"),
                        ], id="lo-topn-div", style={"display": "none"}),

                        dbc.Button("Run Analysis", id="lo-btn-run", color="primary",
                                   className="w-100", disabled=True),
                    ]),
                ], className="mb-3"),
            ], md=4),

            # Right column: results
            dbc.Col([
                dbc.Spinner([
                    html.Div(id="lo-error", className="mb-2"),
                    dcc.Graph(id="lo-plot", style={"display": "none"},
                             config={"toImageButtonOptions": {"format": "svg", "scale": 2}}),
                    html.Div(id="lo-stats-table"),
                    html.Div(id="lo-download-area"),
                ], color="primary"),
            ], md=8),
        ]),
    ], fluid=True)


# ── Callbacks ────────────────────────────────────────────────────────────────


@dash_app.callback(
    Output("lo-biom-path", "data"),
    Output("lo-biom-status", "children"),
    Output("lo-meta-store", "data"),
    Output("lo-sample-id-col", "data"),
    Output("lo-meta-status", "children"),
    Output("lo-subject-col", "options"),
    Output("lo-time-col", "options"),
    Output("lo-group-col", "options"),
    Output("lo-btn-run", "disabled"),
    Input("lo-upload-biom", "contents"),
    Input("lo-select-pipeline", "value"),
    Input("lo-upload-meta", "contents"),
    State("lo-upload-biom", "filename"),
    State("lo-upload-meta", "filename"),
    State("lo-biom-path", "data"),
    State("lo-meta-store", "data"),
    State("lo-sample-id-col", "data"),
    prevent_initial_call=True,
)
def on_input_change(biom_contents, pipeline_value, meta_contents,
                    biom_filename, meta_filename,
                    prev_biom_path, prev_meta_json, prev_sid_col):
    trigger = callback_context.triggered_id

    biom_path = prev_biom_path
    biom_status = no_update
    meta_json = prev_meta_json
    sid_col = prev_sid_col
    meta_status = no_update
    col_opts = no_update
    time_opts = no_update
    group_opts = no_update
    btn_disabled = no_update

    def _build_col_options(df, sid):
        cols = [c for c in df.columns if c != sid]
        return [{"label": c, "value": c} for c in cols]

    # ── Handle BIOM input ──
    if trigger == "lo-select-pipeline" and not pipeline_value:
        return None, "", None, None, "", [], [], [], True

    if trigger == "lo-select-pipeline" and pipeline_value:
        try:
            table = load_table(pipeline_value)
            n = len(table.ids(axis="sample"))
            biom_path = pipeline_value
            biom_status = dbc.Alert(f"Pipeline dataset loaded: {n} samples", color="success")

            db_meta, db_sid = get_dataset_metadata_df(pipeline_value)
            if db_meta is not None:
                meta_json = db_meta.to_json(date_format="iso", orient="split")
                sid_col = db_sid
                opts = _build_col_options(db_meta, db_sid)
                col_opts = opts
                time_opts = opts
                group_opts = [{"label": "— None —", "value": ""}] + opts
                meta_status = dbc.Alert(
                    f"Metadata auto-loaded: {len(db_meta)} samples",
                    color="success",
                )
                btn_disabled = False
            else:
                meta_json = None
                sid_col = None
                col_opts = []
                time_opts = []
                group_opts = []
                meta_status = dbc.Alert(
                    "No metadata found. Upload a metadata file with subject and time columns.",
                    color="info",
                )
                btn_disabled = True
        except Exception as e:
            biom_path = None
            biom_status = dbc.Alert(f"Error: {e}", color="danger")
            btn_disabled = True

    elif trigger == "lo-upload-biom" and biom_contents:
        path, error = parse_uploaded_biom(biom_contents, biom_filename or "")
        if error:
            return None, dbc.Alert(error, color="danger"), prev_meta_json, prev_sid_col, no_update, no_update, no_update, no_update, True
        table = load_table(path)
        sample_ids = list(table.ids(axis="sample"))
        n = len(sample_ids)
        biom_path = path
        biom_status = dbc.Alert(f"Uploaded BIOM: {n} samples", color="success")

        match_df, match_sid, match_name = find_metadata_for_samples(sample_ids)
        if match_df is not None:
            meta_json = match_df.to_json(date_format="iso", orient="split")
            sid_col = match_sid
            opts = _build_col_options(match_df, match_sid)
            col_opts = opts
            time_opts = opts
            group_opts = [{"label": "— None —", "value": ""}] + opts
            meta_status = dbc.Alert(
                f"Metadata auto-matched from \"{match_name}\": {len(match_df)} samples",
                color="success",
            )
            btn_disabled = False
        else:
            btn_disabled = meta_json is not None

    elif trigger == "lo-upload-meta" and meta_contents:
        df, s_col, error = parse_uploaded_metadata(meta_contents, meta_filename or "")
        if error:
            return prev_biom_path, no_update, None, None, dbc.Alert(error, color="danger"), [], [], [], True
        meta_json = df.to_json(date_format="iso", orient="split")
        sid_col = s_col
        opts = _build_col_options(df, s_col)
        col_opts = opts
        time_opts = opts
        group_opts = [{"label": "— None —", "value": ""}] + opts
        meta_status = dbc.Alert(f"Metadata loaded: {len(df)} samples", color="success")
        btn_disabled = biom_path is None

    return biom_path, biom_status, meta_json, sid_col, meta_status, col_opts, time_opts, group_opts, btn_disabled


@dash_app.callback(
    Output("lo-metric", "options"),
    Output("lo-metric", "value"),
    Output("lo-topn-div", "style"),
    Input("lo-analysis-type", "value"),
    prevent_initial_call=True,
)
def on_analysis_type_change(analysis_type):
    if analysis_type == "alpha_traj":
        return ALPHA_METRICS, "shannon", {"display": "none"}
    elif analysis_type == "volatility":
        return DISTANCE_METRICS, "braycurtis", {"display": "none"}
    elif analysis_type == "temporal_da":
        return TAXONOMY_LEVELS, "ASV", {"display": "none"}
    elif analysis_type == "heatmap":
        return TAXONOMY_LEVELS, "Genus", {"display": "block"}
    return ALPHA_METRICS, "shannon", {"display": "none"}


@dash_app.callback(
    Output("lo-btn-run", "disabled", allow_duplicate=True),
    Input("lo-subject-col", "value"),
    Input("lo-time-col", "value"),
    State("lo-biom-path", "data"),
    State("lo-meta-store", "data"),
    prevent_initial_call=True,
)
def on_settings_change(subject_col, time_col, biom_path, meta_json):
    if biom_path and meta_json and subject_col and time_col:
        return False
    return True


@dash_app.callback(
    Output("lo-plot", "figure"),
    Output("lo-plot", "style"),
    Output("lo-stats-table", "children"),
    Output("lo-error", "children"),
    Output("lo-download-area", "children"),
    Output("lo-result-csv", "data"),
    Input("lo-btn-run", "n_clicks"),
    State("lo-biom-path", "data"),
    State("lo-meta-store", "data"),
    State("lo-sample-id-col", "data"),
    State("lo-subject-col", "value"),
    State("lo-time-col", "value"),
    State("lo-group-col", "value"),
    State("lo-analysis-type", "value"),
    State("lo-metric", "value"),
    State("lo-top-n", "value"),
    prevent_initial_call=True,
)
def on_run(n_clicks, biom_path, meta_json, sid_col, subject_col, time_col,
           group_col, analysis_type, metric, top_n):
    if not all([biom_path, meta_json, subject_col, time_col]):
        return no_update, no_update, no_update, dbc.Alert("Please select subject and time columns.", color="warning"), no_update, no_update

    # Treat empty string group_col as None
    if not group_col:
        group_col = None

    try:
        meta_df = pd.read_json(io.StringIO(meta_json), orient="split")
        table = load_table(biom_path)
        biom_ids = list(table.ids(axis="sample"))
        match_info = validate_metadata_vs_biom(meta_df, sid_col, biom_ids)

        if not match_info["matched"]:
            return no_update, no_update, no_update, dbc.Alert("No matching sample IDs.", color="danger"), no_update, no_update

        matched = match_info["matched"]
        count_df = biom_to_count_df(biom_path)

        # Validate longitudinal requirements
        meta_sub = meta_df[meta_df[sid_col].astype(str).isin(matched)]
        n_subjects = meta_sub[subject_col].nunique()
        samples_per_subject = meta_sub.groupby(subject_col).size()

        if n_subjects < 2:
            return no_update, no_update, no_update, dbc.Alert("Need at least 2 subjects for longitudinal analysis.", color="danger"), no_update, no_update
        if (samples_per_subject < 2).all():
            return no_update, no_update, no_update, dbc.Alert("No subjects have repeated measures (≥2 timepoints).", color="danger"), no_update, no_update

        if analysis_type == "alpha_traj":
            return _run_alpha_trajectories(count_df, meta_df, sid_col, subject_col, time_col, group_col, metric, matched)
        elif analysis_type == "volatility":
            return _run_volatility(count_df, meta_df, sid_col, subject_col, time_col, group_col, metric, matched)
        elif analysis_type == "temporal_da":
            return _run_temporal_da(biom_path, meta_df, sid_col, subject_col, time_col, group_col, metric, matched)
        elif analysis_type == "heatmap":
            return _run_heatmap(biom_path, meta_df, sid_col, subject_col, time_col, matched, metric, top_n or 20)

    except Exception as e:
        return no_update, no_update, no_update, dbc.Alert(f"Error: {e}\n{traceback.format_exc()}", color="danger"), no_update, no_update


def _run_alpha_trajectories(count_df, meta_df, sid_col, subject_col, time_col, group_col, metric, matched):
    traj_df = compute_alpha_trajectories(count_df, meta_df, sid_col, subject_col, time_col, metric, matched)

    if traj_df.empty:
        return no_update, no_update, no_update, dbc.Alert("No data to plot.", color="warning"), no_update, no_update

    fig = go.Figure()

    # Get group mapping if available
    if group_col:
        meta_map = meta_df.set_index(meta_df[sid_col].astype(str))
        subject_to_group = {}
        for _, row in meta_df.iterrows():
            subject_to_group[row[subject_col]] = row[group_col]
        groups = sorted(set(subject_to_group.values()))
    else:
        subject_to_group = None
        groups = None

    subjects = sorted(traj_df["subject"].unique())

    for i, subject in enumerate(subjects):
        subj_data = traj_df[traj_df["subject"] == subject].sort_values("time")
        if group_col and subject_to_group:
            grp = subject_to_group.get(subject, "Unknown")
            color_idx = groups.index(grp) if grp in groups else 0
            legend_group = grp
        else:
            color_idx = i
            legend_group = str(subject)

        fig.add_trace(go.Scatter(
            x=subj_data["time"],
            y=subj_data["metric_value"],
            mode="lines+markers",
            name=str(subject),
            legendgroup=legend_group,
            line=dict(color=COLORS[color_idx % len(COLORS)], width=1),
            marker=dict(size=5),
            showlegend=(i < 30),  # limit legend entries
        ))

    # Add group mean lines if grouping
    if group_col and groups:
        for g_idx, grp in enumerate(groups):
            grp_subjects = [s for s, g in subject_to_group.items() if g == grp]
            grp_data = traj_df[traj_df["subject"].isin(grp_subjects)]
            mean_by_time = grp_data.groupby("time")["metric_value"].mean().reset_index()
            mean_by_time = mean_by_time.sort_values("time")
            fig.add_trace(go.Scatter(
                x=mean_by_time["time"],
                y=mean_by_time["metric_value"],
                mode="lines",
                name=f"{grp} (mean)",
                legendgroup=grp,
                line=dict(color=COLORS[g_idx % len(COLORS)], width=4, dash="dash"),
                showlegend=True,
            ))

    metric_label = metric.replace("_", " ").title()
    fig.update_layout(
        template="plotly_dark",
        title=f"{metric_label} Over Time",
        xaxis_title="Time",
        yaxis_title=metric_label,
        height=500,
        legend=dict(orientation="v", x=1.02),
    )

    csv_data = traj_df.to_csv(index=False)
    dl = _download_area(csv_data)
    stats = _trajectory_stats(traj_df, subject_col)

    return fig, {"display": "block"}, stats, "", dl, csv_data


def _run_volatility(count_df, meta_df, sid_col, subject_col, time_col, group_col, metric, matched):
    vol_df = compute_volatility(count_df, meta_df, sid_col, subject_col, time_col, matched, metric)

    if vol_df.empty:
        return no_update, no_update, no_update, dbc.Alert("No consecutive timepoints found.", color="warning"), no_update, no_update

    summary = compute_volatility_summary(vol_df)

    fig = go.Figure()

    if group_col:
        # Map subject to group
        subject_to_group = {}
        for _, row in meta_df.iterrows():
            subject_to_group[row[subject_col]] = row[group_col]
        summary["group"] = summary["subject"].map(subject_to_group)
        groups = sorted(summary["group"].dropna().unique())

        for g_idx, grp in enumerate(groups):
            grp_data = summary[summary["group"] == grp]
            fig.add_trace(go.Box(
                y=grp_data["mean_volatility"],
                name=str(grp),
                marker_color=COLORS[g_idx % len(COLORS)],
                boxpoints="all",
                jitter=0.3,
            ))
    else:
        fig.add_trace(go.Box(
            y=summary["mean_volatility"],
            name="All subjects",
            marker_color=COLORS[0],
            boxpoints="all",
            jitter=0.3,
        ))

    metric_label = metric.replace("_", " ").title()
    fig.update_layout(
        template="plotly_dark",
        title=f"Temporal Volatility ({metric_label})",
        yaxis_title=f"Mean {metric_label} Distance",
        height=450,
    )

    csv_data = summary.to_csv(index=False)
    dl = _download_area(csv_data)

    # Stats summary
    stats_rows = [{
        "Metric": metric_label,
        "N subjects": len(summary),
        "Mean volatility": f"{summary['mean_volatility'].mean():.4f}",
        "Std": f"{summary['mean_volatility'].std():.4f}",
        "Min": f"{summary['mean_volatility'].min():.4f}",
        "Max": f"{summary['mean_volatility'].max():.4f}",
    }]
    stats_table = dbc.Table.from_dataframe(
        pd.DataFrame(stats_rows), striped=True, bordered=True, hover=True,
        color="dark", size="sm",
    )

    return fig, {"display": "block"}, html.Div([html.H5("Volatility Summary", className="mt-3"), stats_table]), "", dl, csv_data


def _run_temporal_da(biom_path, meta_df, sid_col, subject_col, time_col, group_col, level, matched):
    da_df = run_longitudinal_da(
        biom_path, meta_df, sid_col, matched,
        subject_col, time_col, group_col, level,
    )

    if da_df.empty:
        return no_update, no_update, no_update, dbc.Alert("No results returned from MaAsLin2.", color="warning"), no_update, no_update

    # Build volcano-style plot
    da_df["neg_log10_q"] = -np.log10(da_df["qvalue"].clip(lower=1e-300))
    sig = da_df["qvalue"] < 0.05

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=da_df.loc[~sig, "log2fc"],
        y=da_df.loc[~sig, "neg_log10_q"],
        mode="markers",
        name="Not significant",
        marker=dict(color="gray", size=6, opacity=0.5),
    ))
    fig.add_trace(go.Scatter(
        x=da_df.loc[sig, "log2fc"],
        y=da_df.loc[sig, "neg_log10_q"],
        mode="markers",
        name="Significant (q < 0.05)",
        marker=dict(color="#e74c3c", size=8),
        text=da_df.loc[sig, "feature"],
        hovertemplate="%{text}<br>log2FC: %{x:.3f}<br>-log10(q): %{y:.3f}<extra></extra>",
    ))
    fig.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="yellow",
                  annotation_text="q = 0.05")
    fig.update_layout(
        template="plotly_dark",
        title=f"Temporal Differential Abundance ({level})",
        xaxis_title="log2 Fold Change (per time unit)",
        yaxis_title="-log10(q-value)",
        height=500,
    )

    # Summary table
    n_sig = sig.sum()
    n_up = ((da_df["log2fc"] > 0) & sig).sum()
    n_down = ((da_df["log2fc"] < 0) & sig).sum()
    stats_rows = [{
        "Total features": len(da_df),
        "Significant (q<0.05)": n_sig,
        "Increasing over time": n_up,
        "Decreasing over time": n_down,
    }]
    stats_table = dbc.Table.from_dataframe(
        pd.DataFrame(stats_rows), striped=True, bordered=True, hover=True,
        color="dark", size="sm",
    )

    # Top significant features table
    top_sig = da_df[sig].head(20)[["feature", "log2fc", "pvalue", "qvalue"]].copy()
    top_sig["log2fc"] = top_sig["log2fc"].map(lambda x: f"{x:.4f}")
    top_sig["pvalue"] = top_sig["pvalue"].map(lambda x: f"{x:.2e}")
    top_sig["qvalue"] = top_sig["qvalue"].map(lambda x: f"{x:.2e}")

    tables = [html.H5("Summary", className="mt-3"), stats_table]
    if not top_sig.empty:
        tables.append(html.H5("Top Significant Features", className="mt-3"))
        tables.append(dbc.Table.from_dataframe(
            top_sig, striped=True, bordered=True, hover=True,
            color="dark", size="sm",
        ))

    csv_data = da_df.to_csv(index=False)
    dl = _download_area(csv_data)

    return fig, {"display": "block"}, html.Div(tables), "", dl, csv_data


def _run_heatmap(biom_path, meta_df, sid_col, subject_col, time_col, matched, level, top_n):
    heatmap_df = compute_temporal_heatmap(
        biom_path, meta_df, sid_col, subject_col, time_col, matched, level, top_n,
    )

    if heatmap_df.empty:
        return no_update, no_update, no_update, dbc.Alert("No taxonomy data available.", color="warning"), no_update, no_update

    fig = go.Figure(data=go.Heatmap(
        z=heatmap_df.values,
        x=[f"T={t:.1f}" if isinstance(t, float) else str(t) for t in heatmap_df.columns],
        y=heatmap_df.index.tolist(),
        colorscale="Viridis",
        colorbar=dict(title="Rel. Abundance"),
        hoverongaps=False,
    ))
    fig.update_layout(
        template="plotly_dark",
        title=f"Mean Relative Abundance Over Time ({level}, top {top_n})",
        xaxis_title="Timepoint",
        yaxis_title="Taxon",
        height=max(400, 25 * len(heatmap_df)),
    )

    csv_data = heatmap_df.to_csv()
    dl = _download_area(csv_data)

    return fig, {"display": "block"}, "", "", dl, csv_data


def _trajectory_stats(traj_df, subject_col):
    """Build a summary of trajectory data."""
    n_subjects = traj_df["subject"].nunique()
    n_timepoints = traj_df["time"].nunique()
    rows = [{
        "Subjects": n_subjects,
        "Timepoints": n_timepoints,
        "Total observations": len(traj_df),
        "Mean value": f"{traj_df['metric_value'].mean():.4f}",
        "Std": f"{traj_df['metric_value'].std():.4f}",
    }]
    return html.Div([
        html.H5("Summary", className="mt-3"),
        dbc.Table.from_dataframe(
            pd.DataFrame(rows), striped=True, bordered=True, hover=True,
            color="dark", size="sm",
        ),
    ])


def _download_area(csv_data):
    return html.Div([
        dcc.Download(id="lo-download"),
        dbc.Button("Download Results (CSV)", id="lo-btn-download",
                   color="secondary", size="sm", className="mt-2"),
        dcc.Store(id="lo-dl-csv", data=csv_data),
    ])


@dash_app.callback(
    Output("lo-download", "data"),
    Input("lo-btn-download", "n_clicks"),
    State("lo-dl-csv", "data"),
    prevent_initial_call=True,
)
def on_download(n_clicks, csv_data):
    if csv_data:
        return dcc.send_string(csv_data, "longitudinal_results.csv")
    return no_update
