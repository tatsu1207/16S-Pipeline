"""
Random Forest classification page — identify discriminant taxa between groups.
"""
import io
import traceback

import dash_bootstrap_components as dbc
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from biom import load_table
from dash import Input, Output, State, callback_context, dcc, html, no_update
from plotly.subplots import make_subplots

from app.analysis.random_forest import run_random_forest
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
        html.H3("Random Forest", className="mb-2"),
        html.P(
            "Machine learning classification to identify taxa that best discriminate between groups.",
            className="text-muted mb-4",
        ),
        dcc.Store(id="rf-biom-path"),
        dcc.Store(id="rf-meta-store"),
        dcc.Store(id="rf-sample-id-col"),

        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader("Input Data"),
                    dbc.CardBody([
                        dbc.Label("BIOM Table", className="fw-bold"),
                        dbc.Select(
                            id="rf-select-pipeline",
                            options=[{"label": "-- none --", "value": ""}] + pipeline_opts,
                            value="",
                            className="mb-2",
                        ),
                        html.Div("-- or upload --", className="text-center text-muted small mb-2"),
                        dcc.Upload(
                            id="rf-upload-biom",
                            children=html.Div(["Drag & drop or ", html.A("select .biom")]),
                            style={"borderWidth": "2px", "borderStyle": "dashed",
                                   "borderRadius": "5px", "borderColor": "#555",
                                   "textAlign": "center", "padding": "8px"},
                            multiple=False,
                        ),
                        html.Div(id="rf-biom-status", className="mt-1 small mb-3"),

                        dbc.Label("Metadata (CSV/TSV)", className="fw-bold"),
                        html.Div(
                            "Auto-loaded from pipeline dataset if available.",
                            className="text-muted small mb-1",
                        ),
                        dcc.Upload(
                            id="rf-upload-meta",
                            children=html.Div(["Drag & drop or ", html.A("select metadata file")]),
                            style={"borderWidth": "2px", "borderStyle": "dashed",
                                   "borderRadius": "5px", "borderColor": "#555",
                                   "textAlign": "center", "padding": "8px"},
                            multiple=False,
                        ),
                        html.Div(id="rf-meta-status", className="mt-1 small mb-3"),

                        dbc.Label("Group Column", className="fw-bold"),
                        dbc.Select(id="rf-group-col", placeholder="Select metadata...",
                                   className="mb-3"),

                        dbc.Label("Number of Trees", className="fw-bold"),
                        dbc.Input(id="rf-n-trees", type="number", value=500, min=100, max=5000,
                                  step=100, className="mb-3"),

                        dbc.Label("Top Features to Show", className="fw-bold"),
                        dbc.Input(id="rf-top-n", type="number", value=20, min=5, max=100,
                                  className="mb-3"),

                        dbc.Button("Run Classification", id="rf-btn-run", color="primary",
                                   className="w-100", disabled=True),
                    ]),
                ], className="mb-3"),
            ], md=4),

            dbc.Col([
                dbc.Spinner([
                    html.Div(id="rf-error", className="mb-2"),
                    html.Div(id="rf-summary", className="mb-2"),
                    dcc.Graph(id="rf-importance-plot", style={"display": "none"},
                             config={"toImageButtonOptions": {"format": "svg", "scale": 2}}),
                    dcc.Graph(id="rf-confusion-plot", style={"display": "none"},
                             config={"toImageButtonOptions": {"format": "svg", "scale": 2}}),
                    html.Div(id="rf-report-table"),
                    html.Div(id="rf-download-area"),
                ], color="primary"),
            ], md=8),
        ]),
    ], fluid=True)


# ── Input handling ───────────────────────────────────────────────────────────

@dash_app.callback(
    Output("rf-biom-path", "data"),
    Output("rf-biom-status", "children"),
    Output("rf-meta-store", "data"),
    Output("rf-sample-id-col", "data"),
    Output("rf-meta-status", "children"),
    Output("rf-group-col", "options"),
    Output("rf-btn-run", "disabled"),
    Input("rf-upload-biom", "contents"),
    Input("rf-select-pipeline", "value"),
    Input("rf-upload-meta", "contents"),
    State("rf-upload-biom", "filename"),
    State("rf-upload-meta", "filename"),
    State("rf-biom-path", "data"),
    State("rf-meta-store", "data"),
    State("rf-sample-id-col", "data"),
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
    group_opts = no_update
    btn_disabled = no_update

    if trigger == "rf-select-pipeline" and not pipeline_value:
        return None, "", None, None, "", [], True

    if trigger == "rf-select-pipeline" and pipeline_value:
        try:
            table = load_table(pipeline_value)
            n = len(table.ids(axis="sample"))
            biom_path = pipeline_value
            biom_status = dbc.Alert(f"Dataset loaded: {n} samples", color="success")

            db_meta, db_sid = get_dataset_metadata_df(pipeline_value)
            if db_meta is not None:
                meta_json = db_meta.to_json(date_format="iso", orient="split")
                sid_col = db_sid
                gcols = get_group_columns(db_meta, db_sid)
                group_opts = [{"label": c, "value": c} for c in gcols]
                meta_status = dbc.Alert(
                    f"Metadata auto-loaded: {len(db_meta)} samples, {len(gcols)} group columns",
                    color="success",
                )
                btn_disabled = len(gcols) == 0
            else:
                meta_json = None
                sid_col = None
                group_opts = []
                meta_status = dbc.Alert("No metadata found. Upload a metadata file.", color="info")
                btn_disabled = True
        except Exception as e:
            return None, dbc.Alert(f"Error: {e}", color="danger"), None, None, "", [], True

    elif trigger == "rf-upload-biom" and biom_contents:
        path, error = parse_uploaded_biom(biom_contents, biom_filename or "")
        if error:
            return None, dbc.Alert(error, color="danger"), prev_meta_json, prev_sid_col, no_update, no_update, True
        table = load_table(path)
        biom_path = path
        biom_status = dbc.Alert(f"Uploaded: {len(table.ids(axis='sample'))} samples", color="success")

        sample_ids = list(table.ids(axis="sample"))
        match_df, match_sid, match_name = find_metadata_for_samples(sample_ids)
        if match_df is not None:
            meta_json = match_df.to_json(date_format="iso", orient="split")
            sid_col = match_sid
            gcols = get_group_columns(match_df, match_sid)
            group_opts = [{"label": c, "value": c} for c in gcols]
            meta_status = dbc.Alert(
                f"Metadata auto-matched from \"{match_name}\": {len(match_df)} samples",
                color="success",
            )
            btn_disabled = len(gcols) == 0
        else:
            btn_disabled = meta_json is None

    elif trigger == "rf-upload-meta" and meta_contents:
        df, s_col, error = parse_uploaded_metadata(meta_contents, meta_filename or "")
        if error:
            return prev_biom_path, no_update, None, None, dbc.Alert(error, color="danger"), [], True
        meta_json = df.to_json(date_format="iso", orient="split")
        sid_col = s_col
        gcols = get_group_columns(df, s_col)
        group_opts = [{"label": c, "value": c} for c in gcols]
        meta_status = dbc.Alert(f"Metadata loaded: {len(df)} samples, {len(gcols)} group columns", color="success")
        btn_disabled = biom_path is None or len(gcols) == 0

    return biom_path, biom_status, meta_json, sid_col, meta_status, group_opts, btn_disabled


# ── Run classification ───────────────────────────────────────────────────────

@dash_app.callback(
    Output("rf-importance-plot", "figure"),
    Output("rf-importance-plot", "style"),
    Output("rf-confusion-plot", "figure"),
    Output("rf-confusion-plot", "style"),
    Output("rf-summary", "children"),
    Output("rf-report-table", "children"),
    Output("rf-error", "children"),
    Output("rf-download-area", "children"),
    Input("rf-btn-run", "n_clicks"),
    State("rf-biom-path", "data"),
    State("rf-meta-store", "data"),
    State("rf-sample-id-col", "data"),
    State("rf-group-col", "value"),
    State("rf-n-trees", "value"),
    State("rf-top-n", "value"),
    prevent_initial_call=True,
)
def on_run(n_clicks, biom_path, meta_json, sid_col, group_col, n_trees, top_n):
    empty = (no_update, no_update, no_update, no_update, no_update, no_update)
    if not all([biom_path, meta_json, group_col]):
        return *empty, dbc.Alert("Please fill all inputs.", color="warning"), no_update

    try:
        count_df = biom_to_count_df(biom_path)
        meta_df = pd.read_json(io.StringIO(meta_json), orient="split")
        n_trees = int(n_trees or 500)
        top_n = int(top_n or 20)

        result = run_random_forest(
            count_df, meta_df, sid_col, group_col,
            n_estimators=n_trees, top_n=top_n,
        )

        # Feature importance plot
        imp_df = result["importance_df"]
        imp_df = imp_df.iloc[::-1]  # reverse for horizontal bar
        labels = [_shorten_taxon(f) for f in imp_df["Feature"]]

        fig_imp = go.Figure(go.Bar(
            x=imp_df["Importance"].values,
            y=labels,
            orientation="h",
            marker_color="#3498db",
        ))
        fig_imp.update_layout(
            template="plotly_dark",
            height=max(400, 22 * len(labels) + 100),
            xaxis_title="Mean Decrease in Impurity",
            yaxis_title="",
            title=f"Top {top_n} Important Features",
            margin=dict(l=250),
        )

        # Confusion matrix
        cm = np.array(result["confusion_matrix"])
        class_labels = result["class_labels"]

        fig_cm = go.Figure(data=go.Heatmap(
            z=cm,
            x=class_labels,
            y=class_labels,
            colorscale="Blues",
            text=cm.astype(str),
            texttemplate="%{text}",
            hoverinfo="z",
            showscale=False,
        ))
        fig_cm.update_layout(
            template="plotly_dark",
            height=350,
            xaxis_title="Predicted",
            yaxis_title="Actual",
            title="Confusion Matrix (Cross-Validated)",
        )

        # Summary
        summary = dbc.Alert([
            html.Strong(f"Cross-validated accuracy ({result['cv_folds']}-fold): "
                        f"{result['cv_accuracy']*100:.1f}%"),
            html.Br(),
            f"Training accuracy: {result['accuracy']*100:.1f}% | "
            f"Trees: {n_trees}",
        ], color="info")

        # Classification report table
        report = result["classification_report"]
        report_rows = []
        for cls in class_labels:
            if cls in report:
                r = report[cls]
                report_rows.append({
                    "Class": cls,
                    "Precision": f"{r['precision']:.3f}",
                    "Recall": f"{r['recall']:.3f}",
                    "F1-Score": f"{r['f1-score']:.3f}",
                    "Support": int(r["support"]),
                })

        report_table = html.Div([
            html.H5("Classification Report (Cross-Validated)", className="mt-3"),
            dbc.Table.from_dataframe(
                pd.DataFrame(report_rows), striped=True, bordered=True, hover=True,
                color="dark", size="sm",
            ),
        ])

        # Download
        dl = html.Div([
            dcc.Download(id="rf-download"),
            dbc.Button("Download Importances (CSV)", id="rf-btn-download",
                       color="secondary", size="sm", className="mt-2"),
            dcc.Store(id="rf-csv-data", data=result["full_importance_df"].to_csv(index=False)),
        ])

        return (fig_imp, {"display": "block"},
                fig_cm, {"display": "block"},
                summary, report_table, "", dl)

    except Exception as e:
        return *empty, dbc.Alert(f"Error: {e}\n{traceback.format_exc()}", color="danger"), no_update


@dash_app.callback(
    Output("rf-download", "data"),
    Input("rf-btn-download", "n_clicks"),
    State("rf-csv-data", "data"),
    prevent_initial_call=True,
)
def on_download(n_clicks, csv_data):
    if csv_data:
        return dcc.send_string(csv_data, "rf_feature_importance.csv")
    return no_update
