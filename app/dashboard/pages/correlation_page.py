"""
Correlation Heatmap page — taxon-taxon and taxon-metadata correlations.
"""
import io
import traceback

import dash_bootstrap_components as dbc
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from biom import load_table
from dash import Input, Output, State, callback_context, dcc, html, no_update

from app.analysis.correlation import compute_taxon_correlation, compute_taxa_metadata_correlation
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


def get_layout():
    pipeline_opts = get_pipeline_biom_options()
    return dbc.Container([
        html.H3("Correlation Heatmap", className="mb-2"),
        html.P(
            "Visualize Spearman correlations between taxa, or between taxa and numeric metadata variables.",
            className="text-muted mb-4",
        ),
        dcc.Store(id="corr-biom-path"),
        dcc.Store(id="corr-meta-store"),
        dcc.Store(id="corr-sample-id-col"),

        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader("Input Data"),
                    dbc.CardBody([
                        dbc.Label("BIOM Table", className="fw-bold"),
                        dbc.Select(
                            id="corr-select-pipeline",
                            options=[{"label": "-- none --", "value": ""}] + pipeline_opts,
                            value="",
                            className="mb-2",
                        ),
                        html.Div("-- or upload --", className="text-center text-muted small mb-2"),
                        dcc.Upload(
                            id="corr-upload-biom",
                            children=html.Div(["Drag & drop or ", html.A("select .biom")]),
                            style={"borderWidth": "2px", "borderStyle": "dashed",
                                   "borderRadius": "5px", "borderColor": "#555",
                                   "textAlign": "center", "padding": "8px"},
                            multiple=False,
                        ),
                        html.Div(id="corr-biom-status", className="mt-1 small mb-3"),

                        dbc.Label("Metadata (CSV/TSV)", className="fw-bold"),
                        html.Div(
                            "Optional. Upload to enable taxa-metadata correlations.",
                            className="text-muted small mb-1",
                        ),
                        dcc.Upload(
                            id="corr-upload-meta",
                            children=html.Div(["Drag & drop or ", html.A("select metadata file")]),
                            style={"borderWidth": "2px", "borderStyle": "dashed",
                                   "borderRadius": "5px", "borderColor": "#555",
                                   "textAlign": "center", "padding": "8px"},
                            multiple=False,
                        ),
                        html.Div(id="corr-meta-status", className="mt-1 small mb-3"),

                        dbc.Label("Analysis Type", className="fw-bold"),
                        dbc.RadioItems(
                            id="corr-type",
                            options=[
                                {"label": "Taxon-Taxon", "value": "taxon"},
                                {"label": "Taxon-Metadata", "value": "metadata"},
                            ],
                            value="taxon",
                            className="mb-3",
                        ),

                        dbc.Label("Top N Taxa", className="fw-bold"),
                        dbc.Input(id="corr-top-n", type="number", value=30, min=5, max=100,
                                  className="mb-3"),

                        dbc.Button("Run Analysis", id="corr-btn-run", color="primary",
                                   className="w-100", disabled=True),
                    ]),
                ], className="mb-3"),
            ], md=4),

            dbc.Col([
                dbc.Spinner([
                    html.Div(id="corr-error", className="mb-2"),
                    dcc.Graph(id="corr-heatmap", style={"display": "none"},
                             config={"toImageButtonOptions": {"format": "svg", "scale": 2}}),
                ], color="primary"),
            ], md=8),
        ]),
    ], fluid=True)


# ── Input handling callback ──────────────────────────────────────────────────

@dash_app.callback(
    Output("corr-biom-path", "data"),
    Output("corr-biom-status", "children"),
    Output("corr-meta-store", "data"),
    Output("corr-sample-id-col", "data"),
    Output("corr-meta-status", "children"),
    Output("corr-btn-run", "disabled"),
    Input("corr-upload-biom", "contents"),
    Input("corr-select-pipeline", "value"),
    Input("corr-upload-meta", "contents"),
    State("corr-upload-biom", "filename"),
    State("corr-upload-meta", "filename"),
    State("corr-biom-path", "data"),
    State("corr-meta-store", "data"),
    State("corr-sample-id-col", "data"),
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
    btn_disabled = no_update

    if trigger == "corr-select-pipeline" and not pipeline_value:
        return None, "", None, None, "", True

    if trigger == "corr-select-pipeline" and pipeline_value:
        try:
            table = load_table(pipeline_value)
            n = len(table.ids(axis="sample"))
            biom_path = pipeline_value
            biom_status = dbc.Alert(f"Pipeline dataset loaded: {n} samples", color="success")
            btn_disabled = False

            db_meta, db_sid = get_dataset_metadata_df(pipeline_value)
            if db_meta is not None:
                meta_json = db_meta.to_json(date_format="iso", orient="split")
                sid_col = db_sid
                meta_status = dbc.Alert(f"Metadata auto-loaded: {len(db_meta)} samples", color="success")
            else:
                meta_json = None
                sid_col = None
                meta_status = dbc.Alert("No metadata found. Upload for taxa-metadata correlations.", color="info")
        except Exception as e:
            biom_path = None
            biom_status = dbc.Alert(f"Error: {e}", color="danger")
            btn_disabled = True

    elif trigger == "corr-upload-biom" and biom_contents:
        path, error = parse_uploaded_biom(biom_contents, biom_filename or "")
        if error:
            return None, dbc.Alert(error, color="danger"), prev_meta_json, prev_sid_col, no_update, True
        table = load_table(path)
        biom_path = path
        biom_status = dbc.Alert(f"Uploaded BIOM: {len(table.ids(axis='sample'))} samples", color="success")
        btn_disabled = False

        sample_ids = list(table.ids(axis="sample"))
        match_df, match_sid, match_name = find_metadata_for_samples(sample_ids)
        if match_df is not None:
            meta_json = match_df.to_json(date_format="iso", orient="split")
            sid_col = match_sid
            meta_status = dbc.Alert(f"Metadata auto-matched from \"{match_name}\"", color="success")

    elif trigger == "corr-upload-meta" and meta_contents:
        df, s_col, error = parse_uploaded_metadata(meta_contents, meta_filename or "")
        if error:
            return prev_biom_path, no_update, None, None, dbc.Alert(error, color="danger"), prev_biom_path is None
        meta_json = df.to_json(date_format="iso", orient="split")
        sid_col = s_col
        meta_status = dbc.Alert(f"Metadata loaded: {len(df)} samples", color="success")
        btn_disabled = biom_path is None

    return biom_path, biom_status, meta_json, sid_col, meta_status, btn_disabled


# ── Run analysis callback ───────────────────────────────────────────────────

@dash_app.callback(
    Output("corr-heatmap", "figure"),
    Output("corr-heatmap", "style"),
    Output("corr-error", "children"),
    Input("corr-btn-run", "n_clicks"),
    State("corr-biom-path", "data"),
    State("corr-meta-store", "data"),
    State("corr-sample-id-col", "data"),
    State("corr-type", "value"),
    State("corr-top-n", "value"),
    prevent_initial_call=True,
)
def on_run(n_clicks, biom_path, meta_json, sid_col, corr_type, top_n):
    if not biom_path:
        return no_update, no_update, dbc.Alert("Please select a BIOM table.", color="warning")

    try:
        count_df = biom_to_count_df(biom_path)
        top_n = int(top_n or 30)

        if corr_type == "metadata":
            if not meta_json:
                return no_update, no_update, dbc.Alert("Upload metadata for taxa-metadata correlations.", color="warning")

            meta_df = pd.read_json(io.StringIO(meta_json), orient="split")
            meta_cols = [c for c in meta_df.columns if c != sid_col]
            result = compute_taxa_metadata_correlation(count_df, meta_df, sid_col, meta_cols, top_n=top_n)

            if result is None:
                return no_update, no_update, dbc.Alert("No numeric metadata variables with enough data found.", color="warning")

            corr = result["corr_matrix"]
            pvals = result["pval_matrix"]

            # Build hover text with p-values
            hover = []
            for i in range(corr.shape[0]):
                row = []
                for j in range(corr.shape[1]):
                    row.append(f"r={corr.iloc[i,j]:.3f}<br>p={pvals.iloc[i,j]:.4f}")
                hover.append(row)

            # Shorten taxa labels
            labels_y = [_shorten_taxon(t) for t in corr.index]
            labels_x = corr.columns.tolist()

            fig = go.Figure(data=go.Heatmap(
                z=corr.values,
                x=labels_x,
                y=labels_y,
                hovertext=hover,
                hoverinfo="text",
                colorscale="RdBu_r",
                zmid=0,
                zmin=-1, zmax=1,
                colorbar=dict(title="Spearman r"),
            ))
            fig.update_layout(
                template="plotly_dark",
                height=max(400, 20 * len(labels_y) + 150),
                xaxis=dict(tickangle=45),
                margin=dict(l=200, b=120),
            )

        else:
            # Taxon-taxon
            result = compute_taxon_correlation(count_df, top_n=top_n)
            corr = result["corr_matrix"]
            pvals = result["pval_matrix"]

            hover = []
            for i in range(corr.shape[0]):
                row = []
                for j in range(corr.shape[1]):
                    row.append(f"r={corr.iloc[i,j]:.3f}<br>p={pvals.iloc[i,j]:.4f}")
                hover.append(row)

            labels = [_shorten_taxon(t) for t in corr.index]

            fig = go.Figure(data=go.Heatmap(
                z=corr.values,
                x=labels,
                y=labels,
                hovertext=hover,
                hoverinfo="text",
                colorscale="RdBu_r",
                zmid=0,
                zmin=-1, zmax=1,
                colorbar=dict(title="Spearman r"),
            ))
            fig.update_layout(
                template="plotly_dark",
                height=max(500, 18 * len(labels) + 150),
                width=max(500, 18 * len(labels) + 200),
                xaxis=dict(tickangle=45),
                margin=dict(l=200, b=150),
            )

        return fig, {"display": "block"}, ""

    except Exception as e:
        return no_update, no_update, dbc.Alert(f"Error: {e}\n{traceback.format_exc()}", color="danger")


def _shorten_taxon(taxon_id: str) -> str:
    """Shorten a full taxonomy string to the lowest classified rank."""
    if ";" in taxon_id:
        parts = [p.strip() for p in taxon_id.split(";") if p.strip()]
        # Return last non-empty part
        for part in reversed(parts):
            cleaned = part.split("__", 1)[-1] if "__" in part else part
            if cleaned and cleaned.lower() not in ("", "unclassified", "uncultured"):
                return cleaned
        return parts[-1] if parts else taxon_id
    return taxon_id
