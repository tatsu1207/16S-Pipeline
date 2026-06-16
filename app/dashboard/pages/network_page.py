"""
Network Analysis page — microbial co-occurrence network visualization.
"""
import traceback

import dash_bootstrap_components as dbc
import numpy as np
import plotly.graph_objects as go
from biom import load_table
from dash import Input, Output, State, callback_context, dcc, html, no_update

from app.analysis.network import build_network
from app.analysis.shared import (
    biom_to_count_df,
    find_metadata_for_samples,
    get_pipeline_biom_options,
    parse_uploaded_biom,
)
from app.dashboard.app import app as dash_app


def _spring_layout(nodes, edges, iterations=50, k=1.0):
    """Simple force-directed layout (Fruchterman-Reingold)."""
    n = len(nodes)
    if n == 0:
        return {}

    id_to_idx = {node["id"]: i for i, node in enumerate(nodes)}
    pos = np.random.RandomState(42).randn(n, 2)

    for it in range(iterations):
        disp = np.zeros((n, 2))
        temp = k * (1 - it / iterations)

        # Repulsive forces
        for i in range(n):
            for j in range(i + 1, n):
                delta = pos[i] - pos[j]
                dist = max(np.linalg.norm(delta), 0.01)
                force = k * k / dist
                direction = delta / dist
                disp[i] += direction * force
                disp[j] -= direction * force

        # Attractive forces (edges)
        for edge in edges:
            i = id_to_idx.get(edge["source"])
            j = id_to_idx.get(edge["target"])
            if i is None or j is None:
                continue
            delta = pos[j] - pos[i]
            dist = max(np.linalg.norm(delta), 0.01)
            force = dist * dist / k
            direction = delta / dist
            disp[i] += direction * force
            disp[j] -= direction * force

        # Apply displacement with temperature
        for i in range(n):
            mag = max(np.linalg.norm(disp[i]), 0.01)
            pos[i] += disp[i] / mag * min(mag, temp)

    return {nodes[i]["id"]: pos[i] for i in range(n)}


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
        html.H3("Network Analysis", className="mb-2"),
        html.P(
            "Infer microbial co-occurrence networks using correlation-based methods.",
            className="text-muted mb-4",
        ),
        dcc.Store(id="net-biom-path"),

        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader("Input Data"),
                    dbc.CardBody([
                        dbc.Label("BIOM Table", className="fw-bold"),
                        dbc.Select(
                            id="net-select-pipeline",
                            options=[{"label": "-- none --", "value": ""}] + pipeline_opts,
                            value="",
                            className="mb-2",
                        ),
                        html.Div("-- or upload --", className="text-center text-muted small mb-2"),
                        dcc.Upload(
                            id="net-upload-biom",
                            children=html.Div(["Drag & drop or ", html.A("select .biom")]),
                            style={"borderWidth": "2px", "borderStyle": "dashed",
                                   "borderRadius": "5px", "borderColor": "#555",
                                   "textAlign": "center", "padding": "8px"},
                            multiple=False,
                        ),
                        html.Div(id="net-biom-status", className="mt-1 small mb-3"),

                        dbc.Label("Method", className="fw-bold"),
                        dbc.RadioItems(
                            id="net-method",
                            options=[
                                {"label": "SparCC", "value": "sparcc"},
                                {"label": "Spearman", "value": "spearman"},
                            ],
                            value="sparcc",
                            className="mb-3",
                        ),

                        dbc.Label("Top N Taxa", className="fw-bold"),
                        dbc.Input(id="net-top-n", type="number", value=40, min=5, max=100,
                                  className="mb-3"),

                        dbc.Label("Correlation Threshold", className="fw-bold"),
                        dbc.Input(id="net-threshold", type="number", value=0.3, min=0.1,
                                  max=0.9, step=0.05, className="mb-3"),

                        dbc.Button("Build Network", id="net-btn-run", color="primary",
                                   className="w-100", disabled=True),
                    ]),
                ], className="mb-3"),
            ], md=3),

            dbc.Col([
                dbc.Spinner([
                    html.Div(id="net-error", className="mb-2"),
                    html.Div(id="net-summary", className="mb-2"),
                    dcc.Graph(id="net-graph", style={"display": "none"},
                             config={"toImageButtonOptions": {"format": "svg", "scale": 2}}),
                ], color="primary"),
            ], md=9),
        ]),
    ], fluid=True)


# ── Input handling ───────────────────────────────────────────────────────────

@dash_app.callback(
    Output("net-biom-path", "data"),
    Output("net-biom-status", "children"),
    Output("net-btn-run", "disabled"),
    Input("net-upload-biom", "contents"),
    Input("net-select-pipeline", "value"),
    State("net-upload-biom", "filename"),
    State("net-biom-path", "data"),
    prevent_initial_call=True,
)
def on_input_change(biom_contents, pipeline_value, biom_filename, prev_biom_path):
    trigger = callback_context.triggered_id

    if trigger == "net-select-pipeline":
        if not pipeline_value:
            return None, "", True
        try:
            table = load_table(pipeline_value)
            n = len(table.ids(axis="sample"))
            return pipeline_value, dbc.Alert(f"Dataset loaded: {n} samples", color="success"), False
        except Exception as e:
            return None, dbc.Alert(f"Error: {e}", color="danger"), True

    if trigger == "net-upload-biom" and biom_contents:
        path, error = parse_uploaded_biom(biom_contents, biom_filename or "")
        if error:
            return None, dbc.Alert(error, color="danger"), True
        table = load_table(path)
        return path, dbc.Alert(f"Uploaded: {len(table.ids(axis='sample'))} samples", color="success"), False

    return prev_biom_path, no_update, no_update


# ── Build network ────────────────────────────────────────────────────────────

@dash_app.callback(
    Output("net-graph", "figure"),
    Output("net-graph", "style"),
    Output("net-summary", "children"),
    Output("net-error", "children"),
    Input("net-btn-run", "n_clicks"),
    State("net-biom-path", "data"),
    State("net-method", "value"),
    State("net-top-n", "value"),
    State("net-threshold", "value"),
    prevent_initial_call=True,
)
def on_run(n_clicks, biom_path, method, top_n, threshold):
    if not biom_path:
        return no_update, no_update, no_update, dbc.Alert("Select a BIOM table.", color="warning")

    try:
        count_df = biom_to_count_df(biom_path)
        top_n = int(top_n or 40)
        threshold = float(threshold or 0.3)

        result = build_network(count_df, top_n=top_n, corr_threshold=threshold, method=method)
        nodes = result["nodes"]
        edges = result["edges"]

        if not nodes:
            return no_update, no_update, no_update, dbc.Alert("No taxa found.", color="warning")

        # Layout
        pos = _spring_layout(nodes, edges)

        # Scale node sizes by abundance
        abundances = [n["abundance"] for n in nodes]
        max_ab = max(abundances) if abundances else 1
        node_sizes = [max(8, 40 * (a / max_ab)) for a in abundances]

        # Draw edges
        edge_traces = []
        for edge in edges:
            x0, y0 = pos[edge["source"]]
            x1, y1 = pos[edge["target"]]
            color = "rgba(46, 204, 113, 0.4)" if edge["positive"] else "rgba(231, 76, 60, 0.4)"
            width = max(1, abs(edge["weight"]) * 3)
            edge_traces.append(go.Scatter(
                x=[x0, x1, None], y=[y0, y1, None],
                mode="lines",
                line=dict(width=width, color=color),
                hoverinfo="none",
                showlegend=False,
            ))

        # Draw nodes
        node_x = [pos[n["id"]][0] for n in nodes]
        node_y = [pos[n["id"]][1] for n in nodes]
        node_labels = [_shorten_taxon(n["id"]) for n in nodes]
        node_hover = [
            f"{_shorten_taxon(n['id'])}<br>Rel. abundance: {n['abundance']:.4f}"
            for n in nodes
        ]

        node_trace = go.Scatter(
            x=node_x, y=node_y,
            mode="markers+text",
            marker=dict(
                size=node_sizes,
                color="#3498db",
                line=dict(width=1, color="#fff"),
            ),
            text=node_labels,
            textposition="top center",
            textfont=dict(size=9),
            hovertext=node_hover,
            hoverinfo="text",
            showlegend=False,
        )

        fig = go.Figure(data=edge_traces + [node_trace])

        # Legend for edge colors
        fig.add_trace(go.Scatter(
            x=[None], y=[None], mode="lines",
            line=dict(color="rgba(46, 204, 113, 0.8)", width=2),
            name="Positive correlation",
        ))
        fig.add_trace(go.Scatter(
            x=[None], y=[None], mode="lines",
            line=dict(color="rgba(231, 76, 60, 0.8)", width=2),
            name="Negative correlation",
        ))

        fig.update_layout(
            template="plotly_dark",
            height=700,
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            legend=dict(orientation="h", yanchor="bottom", y=1.02),
        )

        n_pos = sum(1 for e in edges if e["positive"])
        n_neg = len(edges) - n_pos
        summary = dbc.Alert(
            f"Network: {len(nodes)} nodes, {len(edges)} edges "
            f"({n_pos} positive, {n_neg} negative) | "
            f"Threshold: |r| >= {threshold}",
            color="info",
        )

        return fig, {"display": "block"}, summary, ""

    except Exception as e:
        return no_update, no_update, no_update, dbc.Alert(f"Error: {e}\n{traceback.format_exc()}", color="danger")
