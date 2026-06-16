"""
Association Biplot page — CCA and db-RDA constrained ordination.
"""
import io
import traceback

import dash_bootstrap_components as dbc
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from biom import load_table
from dash import Input, Output, State, callback_context, dcc, html, no_update

from app.analysis.association import run_association
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


def _shorten_taxon(taxon_id: str) -> str:
    if ";" in taxon_id:
        parts = [p.strip() for p in taxon_id.split(";") if p.strip()]
        for part in reversed(parts):
            cleaned = part.split("__", 1)[-1] if "__" in part else part
            if cleaned and cleaned.lower() not in ("", "unclassified", "uncultured"):
                return cleaned
        return parts[-1] if parts else taxon_id
    return taxon_id


def get_layout():
    pipeline_opts = get_pipeline_biom_options()
    return dbc.Container([
        html.H3("Association Biplot", className="mb-2"),
        html.P(
            "Constrained ordination linking community composition to environmental variables "
            "(db-RDA or CCA).",
            className="text-muted mb-4",
        ),
        dcc.Store(id="assoc-biom-path"),
        dcc.Store(id="assoc-meta-store"),
        dcc.Store(id="assoc-sample-id-col"),

        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader("Input Data"),
                    dbc.CardBody([
                        dbc.Label("BIOM Table", className="fw-bold"),
                        dbc.Select(
                            id="assoc-select-pipeline",
                            options=[{"label": "-- none --", "value": ""}] + pipeline_opts,
                            value="",
                            className="mb-2",
                        ),
                        html.Div("-- or upload --", className="text-center text-muted small mb-2"),
                        dcc.Upload(
                            id="assoc-upload-biom",
                            children=html.Div(["Drag & drop or ", html.A("select .biom")]),
                            style={"borderWidth": "2px", "borderStyle": "dashed",
                                   "borderRadius": "5px", "borderColor": "#555",
                                   "textAlign": "center", "padding": "8px"},
                            multiple=False,
                        ),
                        html.Div(id="assoc-biom-status", className="mt-1 small mb-3"),

                        dbc.Label("Metadata (CSV/TSV)", className="fw-bold"),
                        html.Div(
                            "Required. Must contain numeric environmental variables.",
                            className="text-muted small mb-1",
                        ),
                        dcc.Upload(
                            id="assoc-upload-meta",
                            children=html.Div(["Drag & drop or ", html.A("select metadata file")]),
                            style={"borderWidth": "2px", "borderStyle": "dashed",
                                   "borderRadius": "5px", "borderColor": "#555",
                                   "textAlign": "center", "padding": "8px"},
                            multiple=False,
                        ),
                        html.Div(id="assoc-meta-status", className="mt-1 small mb-3"),

                        dbc.Label("Method", className="fw-bold"),
                        dbc.RadioItems(
                            id="assoc-method",
                            options=[
                                {"label": "db-RDA (Bray-Curtis)", "value": "dbRDA"},
                                {"label": "CCA", "value": "CCA"},
                            ],
                            value="dbRDA",
                            className="mb-3",
                        ),

                        dbc.Label("Environmental Variables", className="fw-bold"),
                        dbc.Checklist(id="assoc-env-cols", className="mb-3"),

                        dbc.Label("Group Column (for coloring)", className="fw-bold"),
                        dbc.Select(id="assoc-group-col", placeholder="Optional",
                                   className="mb-3"),

                        dbc.Button("Run Analysis", id="assoc-btn-run", color="primary",
                                   className="w-100", disabled=True),
                    ]),
                ], className="mb-3"),
            ], md=4),

            dbc.Col([
                dbc.Spinner([
                    html.Div(id="assoc-error", className="mb-2"),
                    dcc.Graph(id="assoc-biplot", style={"display": "none"},
                             config={"toImageButtonOptions": {"format": "svg", "scale": 2}}),
                ], color="primary"),
            ], md=8),
        ]),
    ], fluid=True)


# ── Input handling ───────────────────────────────────────────────────────────

@dash_app.callback(
    Output("assoc-biom-path", "data"),
    Output("assoc-biom-status", "children"),
    Output("assoc-meta-store", "data"),
    Output("assoc-sample-id-col", "data"),
    Output("assoc-meta-status", "children"),
    Output("assoc-env-cols", "options"),
    Output("assoc-env-cols", "value"),
    Output("assoc-group-col", "options"),
    Output("assoc-btn-run", "disabled"),
    Input("assoc-upload-biom", "contents"),
    Input("assoc-select-pipeline", "value"),
    Input("assoc-upload-meta", "contents"),
    State("assoc-upload-biom", "filename"),
    State("assoc-upload-meta", "filename"),
    State("assoc-biom-path", "data"),
    State("assoc-meta-store", "data"),
    State("assoc-sample-id-col", "data"),
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
    env_opts = no_update
    env_vals = no_update
    group_opts = no_update
    btn_disabled = no_update

    if trigger == "assoc-select-pipeline" and not pipeline_value:
        return None, "", None, None, "", [], [], [], True

    if trigger == "assoc-select-pipeline" and pipeline_value:
        try:
            table = load_table(pipeline_value)
            n = len(table.ids(axis="sample"))
            biom_path = pipeline_value
            biom_status = dbc.Alert(f"Dataset loaded: {n} samples", color="success")

            db_meta, db_sid = get_dataset_metadata_df(pipeline_value)
            if db_meta is not None:
                meta_json = db_meta.to_json(date_format="iso", orient="split")
                sid_col = db_sid
                env_opts, env_vals, group_opts = _meta_to_options(db_meta, db_sid)
                meta_status = dbc.Alert(f"Metadata auto-loaded: {len(db_meta)} samples", color="success")
                btn_disabled = len(env_vals) == 0
            else:
                meta_json = None
                sid_col = None
                env_opts, env_vals, group_opts = [], [], []
                meta_status = dbc.Alert("No metadata found. Upload metadata with numeric variables.", color="info")
                btn_disabled = True
        except Exception as e:
            return None, dbc.Alert(f"Error: {e}", color="danger"), None, None, "", [], [], [], True

    elif trigger == "assoc-upload-biom" and biom_contents:
        path, error = parse_uploaded_biom(biom_contents, biom_filename or "")
        if error:
            return None, dbc.Alert(error, color="danger"), prev_meta_json, prev_sid_col, no_update, no_update, no_update, no_update, True
        table = load_table(path)
        biom_path = path
        biom_status = dbc.Alert(f"Uploaded: {len(table.ids(axis='sample'))} samples", color="success")

        sample_ids = list(table.ids(axis="sample"))
        match_df, match_sid, match_name = find_metadata_for_samples(sample_ids)
        if match_df is not None:
            meta_json = match_df.to_json(date_format="iso", orient="split")
            sid_col = match_sid
            env_opts, env_vals, group_opts = _meta_to_options(match_df, match_sid)
            meta_status = dbc.Alert(f"Metadata auto-matched from \"{match_name}\"", color="success")
            btn_disabled = len(env_vals) == 0
        else:
            btn_disabled = True

    elif trigger == "assoc-upload-meta" and meta_contents:
        df, s_col, error = parse_uploaded_metadata(meta_contents, meta_filename or "")
        if error:
            return prev_biom_path, no_update, None, None, dbc.Alert(error, color="danger"), [], [], [], True
        meta_json = df.to_json(date_format="iso", orient="split")
        sid_col = s_col
        env_opts, env_vals, group_opts = _meta_to_options(df, s_col)
        meta_status = dbc.Alert(f"Metadata loaded: {len(df)} samples", color="success")
        btn_disabled = biom_path is None or len(env_vals) == 0

    return biom_path, biom_status, meta_json, sid_col, meta_status, env_opts, env_vals, group_opts, btn_disabled


def _meta_to_options(meta_df, sid_col):
    """Extract numeric columns for env variables and all columns for grouping."""
    all_cols = get_group_columns(meta_df, sid_col)
    numeric_cols = []
    for c in all_cols:
        numeric = pd.to_numeric(meta_df[c], errors="coerce")
        if numeric.notna().sum() >= 3:
            numeric_cols.append(c)

    env_opts = [{"label": c, "value": c} for c in numeric_cols]
    env_vals = numeric_cols
    group_opts = [{"label": "-- none --", "value": ""}] + [{"label": c, "value": c} for c in all_cols]
    return env_opts, env_vals, group_opts


# ── Run analysis ─────────────────────────────────────────────────────────────

@dash_app.callback(
    Output("assoc-biplot", "figure"),
    Output("assoc-biplot", "style"),
    Output("assoc-error", "children"),
    Input("assoc-btn-run", "n_clicks"),
    State("assoc-biom-path", "data"),
    State("assoc-meta-store", "data"),
    State("assoc-sample-id-col", "data"),
    State("assoc-method", "value"),
    State("assoc-env-cols", "value"),
    State("assoc-group-col", "value"),
    prevent_initial_call=True,
)
def on_run(n_clicks, biom_path, meta_json, sid_col, method, env_cols, group_col):
    if not biom_path or not meta_json:
        return no_update, no_update, dbc.Alert("Please provide BIOM and metadata.", color="warning")

    if not env_cols:
        return no_update, no_update, dbc.Alert("Select at least one environmental variable.", color="warning")

    try:
        count_df = biom_to_count_df(biom_path)
        meta_df = pd.read_json(io.StringIO(meta_json), orient="split")

        result = run_association(
            count_df, meta_df, sid_col,
            env_cols=env_cols,
            group_col=group_col if group_col else None,
            method=method,
        )

        sample_scores = result["sample_scores"]
        env_scores = result["env_scores"]
        prop = result["prop_explained"]

        axis1_label = f"{method} 1 ({prop[0]*100:.1f}%)" if len(prop) > 0 else f"{method} 1"
        axis2_label = f"{method} 2 ({prop[1]*100:.1f}%)" if len(prop) > 1 else f"{method} 2"

        fig = go.Figure()

        colors = ["#3498db", "#e74c3c", "#2ecc71", "#f39c12", "#9b59b6",
                  "#1abc9c", "#e67e22", "#34495e"]

        # Plot samples
        if result["groups"]:
            groups = result["groups"]
            unique_groups = sorted(set(groups))
            for g_idx, group in enumerate(unique_groups):
                mask = [i for i, g in enumerate(groups) if g == group]
                fig.add_trace(go.Scatter(
                    x=sample_scores[mask, 0],
                    y=sample_scores[mask, 1],
                    mode="markers",
                    name=str(group),
                    marker=dict(size=10, color=colors[g_idx % len(colors)]),
                    text=[result["sample_ids"][i] for i in mask],
                    hoverinfo="text",
                ))
        else:
            fig.add_trace(go.Scatter(
                x=sample_scores[:, 0],
                y=sample_scores[:, 1],
                mode="markers",
                name="Samples",
                marker=dict(size=10, color="#3498db"),
                text=result["sample_ids"],
                hoverinfo="text",
            ))

        # Plot environmental variable arrows
        env_labels = result["env_labels"]
        # Scale arrows to fit the sample space
        sample_range = max(
            np.ptp(sample_scores[:, 0]) if sample_scores.shape[1] > 0 else 1,
            np.ptp(sample_scores[:, 1]) if sample_scores.shape[1] > 1 else 1,
        )
        env_max = np.max(np.abs(env_scores)) if env_scores.size > 0 else 1
        scale = sample_range * 0.4 / max(env_max, 1e-10)

        for i, label in enumerate(env_labels):
            ex = env_scores[i, 0] * scale if env_scores.shape[1] > 0 else 0
            ey = env_scores[i, 1] * scale if env_scores.shape[1] > 1 else 0
            fig.add_annotation(
                x=ex, y=ey, ax=0, ay=0,
                xref="x", yref="y", axref="x", ayref="y",
                showarrow=True,
                arrowhead=2, arrowsize=1.5, arrowwidth=2,
                arrowcolor="#e74c3c",
            )
            fig.add_trace(go.Scatter(
                x=[ex * 1.1], y=[ey * 1.1],
                mode="text",
                text=[label],
                textfont=dict(color="#e74c3c", size=12),
                showlegend=False,
                hoverinfo="skip",
            ))

        # Plot species scores (CCA only)
        if method == "CCA" and "species_scores" in result and "species_labels" in result:
            sp_scores = result["species_scores"]
            sp_labels = [_shorten_taxon(s) for s in result["species_labels"]]
            fig.add_trace(go.Scatter(
                x=sp_scores[:, 0],
                y=sp_scores[:, 1] if sp_scores.shape[1] > 1 else np.zeros(len(sp_scores)),
                mode="markers+text",
                name="Species",
                marker=dict(size=6, color="#95a5a6", symbol="diamond"),
                text=sp_labels,
                textposition="top center",
                textfont=dict(size=8, color="#95a5a6"),
                hoverinfo="text",
            ))

        fig.update_layout(
            template="plotly_dark",
            height=650,
            xaxis_title=axis1_label,
            yaxis_title=axis2_label,
            legend=dict(orientation="h", yanchor="bottom", y=1.02),
        )

        return fig, {"display": "block"}, ""

    except Exception as e:
        return no_update, no_update, dbc.Alert(f"Error: {e}\n{traceback.format_exc()}", color="danger")
