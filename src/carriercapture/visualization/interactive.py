"""
Interactive Dash application for CarrierCapture visualization.

Provides a web-based dashboard for:
- Loading and visualizing potential energy surfaces
- Interactive parameter adjustment and real-time fitting
- Parameter scanning with 2D heatmap visualization
- Multi-potential comparison mode
- Capture coefficient calculation
- Export functionality (JSON, HDF5, PNG)
"""

from typing import Optional, Dict, Any, List
from pathlib import Path
import json
import io
import base64
import warnings

import numpy as np
import dash
from dash import dcc, html, Input, Output, State, callback, ctx
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from carriercapture.core.potential import Potential
from carriercapture.core.config_coord import ConfigCoordinate
from carriercapture.analysis.parameter_scan import ScanResult
from carriercapture.io.readers import load_potential_from_file
from carriercapture.visualization.static import (
    plot_potential,
    plot_capture_coefficient,
    plot_configuration_coordinate,
)
from carriercapture.visualization.themes import create_dash_theme, COLORS


def create_app(debug: bool = False) -> dash.Dash:
    """
    Create and configure the Dash application.

    Parameters
    ----------
    debug : bool, default=False
        Whether to run in debug mode

    Returns
    -------
    dash.Dash
        Configured Dash app

    Examples
    --------
    >>> app = create_app()
    >>> app.run(port=8050)
    """
    app = dash.Dash(
        __name__,
        title="CarrierCapture Visualization",
        suppress_callback_exceptions=True,
    )

    theme = create_dash_theme()

    # App layout with tabs
    app.layout = html.Div(
        style=theme["container"],
        children=[
            # Header
            html.Div(
                style=theme["card"],
                children=[
                    html.H1(
                        "CarrierCapture Visualization Dashboard",
                        style=theme["header"],
                    ),
                    html.P(
                        "Interactive visualization tool for carrier capture rate calculations",
                        style=theme["text"],
                    ),
                ],
            ),
            # Tabs
            dcc.Tabs(
                id="main-tabs",
                value="potential-tab",
                children=[
                    dcc.Tab(
                        label="Potential Viewer",
                        value="potential-tab",
                        children=[create_potential_tab(theme)],
                    ),
                    dcc.Tab(
                        label="Parameter Scan",
                        value="scan-tab",
                        children=[create_scan_tab(theme)],
                    ),
                    dcc.Tab(
                        label="Comparison",
                        value="comparison-tab",
                        children=[create_comparison_tab(theme)],
                    ),
                    dcc.Tab(
                        label="Capture Calculation",
                        value="capture-tab",
                        children=[create_capture_tab(theme)],
                    ),
                ],
            ),
            # Status bar
            html.Div(
                style={**theme["card"], "marginTop": "20px"},
                children=[
                    html.Div(id="status-bar", style=theme["text"]),
                ],
            ),
            # Hidden stores for data
            dcc.Store(id="potential-store-initial"),
            dcc.Store(id="potential-store-final"),
            dcc.Store(id="scan-results-store"),
            dcc.Store(id="comparison-potentials-store", data=[]),
        ],
    )

    # Register all callbacks
    register_callbacks(app)

    return app


def create_potential_tab(theme: Dict[str, Any]) -> html.Div:
    """Create the potential viewer tab layout."""
    return html.Div(
        style={"marginTop": "20px"},
        children=[
            html.Div(
                style={"display": "flex", "gap": "20px"},
                children=[
                    # Sidebar - Controls
                    html.Div(
                        style={**theme["card"], "width": "350px", "height": "fit-content"},
                        children=[
                            html.H2("Controls", style=theme["subheader"]),

                            # File upload
                            html.Label("Load Potential Data:", style=theme["text"]),
                            dcc.Upload(
                                id="upload-potential-data",
                                children=html.Div([
                                    "Drag and Drop or ",
                                    html.A("Select File", style={"color": COLORS["primary"]}),
                                ]),
                                style={
                                    "width": "100%",
                                    "height": "60px",
                                    "lineHeight": "60px",
                                    "borderWidth": "1px",
                                    "borderStyle": "dashed",
                                    "borderRadius": "5px",
                                    "textAlign": "center",
                                    "margin": "10px 0",
                                },
                                multiple=False,
                            ),
                            html.Div(id="upload-status", style={**theme["text"], "marginTop": "10px"}),

                            html.Hr(),

                            # Harmonic potential generator
                            html.H3("Or Generate Harmonic", style=theme["subheader"]),
                            html.Label("ℏω (eV):", style=theme["text"]),
                            dcc.Input(
                                id="harmonic-hw",
                                type="number",
                                value=0.008,
                                step=0.001,
                                style={**theme["input"], "marginBottom": "10px"},
                            ),
                            html.Label("Q₀ (amu^0.5·Å):", style=theme["text"]),
                            dcc.Input(
                                id="harmonic-q0",
                                type="number",
                                value=0.0,
                                step=0.1,
                                style={**theme["input"], "marginBottom": "10px"},
                            ),
                            html.Label("E₀ (eV):", style=theme["text"]),
                            dcc.Input(
                                id="harmonic-e0",
                                type="number",
                                value=0.0,
                                step=0.1,
                                style={**theme["input"], "marginBottom": "10px"},
                            ),
                            html.Button(
                                "Generate Harmonic",
                                id="generate-harmonic-button",
                                n_clicks=0,
                                style=theme["button"],
                            ),

                            html.Hr(),

                            # Fitting controls
                            html.H3("Fitting", style=theme["subheader"]),
                            html.Label("Fit Type:", style=theme["text"]),
                            dcc.Dropdown(
                                id="fit-type",
                                options=[
                                    {"label": "Spline", "value": "spline"},
                                    {"label": "Harmonic", "value": "harmonic"},
                                    {"label": "Morse", "value": "morse"},
                                    {"label": "Polynomial", "value": "polynomial"},
                                ],
                                value="spline",
                                style={"marginBottom": "10px"},
                            ),
                            html.Div(
                                id="fit-params-container",
                                children=[
                                    html.Label("Spline Order:", style=theme["text"]),
                                    dcc.Slider(
                                        id="spline-order",
                                        min=2,
                                        max=5,
                                        step=1,
                                        value=4,
                                        marks={i: str(i) for i in range(2, 6)},
                                        tooltip={"placement": "bottom"},
                                    ),
                                    html.Label("Smoothness:", style=theme["text"]),
                                    dcc.Slider(
                                        id="smoothness",
                                        min=-4,
                                        max=-1,
                                        step=0.5,
                                        value=-3,
                                        marks={i: f"10^{i}" for i in range(-4, 0)},
                                        tooltip={"placement": "bottom"},
                                    ),
                                ],
                            ),
                            html.Button(
                                "Fit Potential",
                                id="fit-button",
                                n_clicks=0,
                                style=theme["button"],
                            ),

                            html.Hr(),

                            # Solving controls
                            html.H3("Schrödinger Solver", style=theme["subheader"]),
                            html.Label("Number of Eigenvalues:", style=theme["text"]),
                            dcc.Input(
                                id="nev",
                                type="number",
                                value=60,
                                min=1,
                                max=500,
                                style={**theme["input"], "marginBottom": "10px"},
                            ),
                            html.Button(
                                "Solve Schrödinger",
                                id="solve-button",
                                n_clicks=0,
                                style=theme["button"],
                            ),

                            html.Hr(),

                            # Display controls
                            html.H3("Display Options", style=theme["subheader"]),
                            dcc.Checklist(
                                id="display-options",
                                options=[
                                    {"label": " Show Wavefunctions", "value": "wf"},
                                    {"label": " Show Eigenvalues", "value": "ev"},
                                    {"label": " Show Data Points", "value": "data"},
                                ],
                                value=["ev", "data"],
                                style=theme["text"],
                            ),
                            html.Label("Wavefunction Scaling:", style=theme["text"]),
                            dcc.Slider(
                                id="wf-scale",
                                min=0.1,
                                max=2.0,
                                step=0.1,
                                value=1.0,
                                marks={i/10: f"{i/10:.1f}" for i in range(2, 21, 4)},
                                tooltip={"placement": "bottom"},
                            ),
                        ],
                    ),
                    # Main plot area
                    html.Div(
                        style={**theme["card"], "flex": "1"},
                        children=[
                            html.H2("Potential Energy Surface", style=theme["subheader"]),
                            dcc.Graph(
                                id="potential-plot",
                                style={"height": "600px"},
                                figure=go.Figure(
                                    layout=dict(
                                        title="Load potential data or generate harmonic to begin",
                                        xaxis_title="Q (amu^0.5·Å)",
                                        yaxis_title="Energy (eV)",
                                        template="plotly_white",
                                    )
                                ),
                            ),
                            html.Div(id="potential-info", style=theme["text"]),
                        ],
                    ),
                ],
            ),
        ],
    )


def create_scan_tab(theme: Dict[str, Any]) -> html.Div:
    """Create the parameter scan tab layout."""
    return html.Div(
        style={"marginTop": "20px"},
        children=[
            html.Div(
                style={"display": "flex", "gap": "20px"},
                children=[
                    # Controls
                    html.Div(
                        style={**theme["card"], "width": "350px", "height": "fit-content"},
                        children=[
                            html.H2("Scan Parameters", style=theme["subheader"]),

                            # Load existing scan results
                            html.Label("Load Scan Results:", style=theme["text"]),
                            dcc.Upload(
                                id="upload-scan-results",
                                children=html.Div([
                                    "Load .npz or .h5 file",
                                ]),
                                style={
                                    "width": "100%",
                                    "height": "60px",
                                    "lineHeight": "60px",
                                    "borderWidth": "1px",
                                    "borderStyle": "dashed",
                                    "borderRadius": "5px",
                                    "textAlign": "center",
                                    "margin": "10px 0",
                                },
                            ),

                            html.Hr(),

                            # ΔQ range
                            html.H3("ΔQ Range", style=theme["subheader"]),
                            html.Label("Min (amu^0.5·Å):", style=theme["text"]),
                            dcc.Input(id="scan-dq-min", type="number", value=0, step=0.5, style=theme["input"]),
                            html.Label("Max (amu^0.5·Å):", style=theme["text"]),
                            dcc.Input(id="scan-dq-max", type="number", value=25, step=0.5, style=theme["input"]),
                            html.Label("Points:", style=theme["text"]),
                            dcc.Input(id="scan-dq-points", type="number", value=25, min=3, style=theme["input"]),

                            html.Hr(),

                            # ΔE range
                            html.H3("ΔE Range", style=theme["subheader"]),
                            html.Label("Min (eV):", style=theme["text"]),
                            dcc.Input(id="scan-de-min", type="number", value=0, step=0.1, style=theme["input"]),
                            html.Label("Max (eV):", style=theme["text"]),
                            dcc.Input(id="scan-de-max", type="number", value=2.5, step=0.1, style=theme["input"]),
                            html.Label("Points:", style=theme["text"]),
                            dcc.Input(id="scan-de-points", type="number", value=10, min=3, style=theme["input"]),

                            html.Hr(),

                            # Other parameters
                            html.H3("Parameters", style=theme["subheader"]),
                            html.Label("ℏω_i (eV):", style=theme["text"]),
                            dcc.Input(id="scan-hw-i", type="number", value=0.008, step=0.001, style=theme["input"]),
                            html.Label("ℏω_f (eV):", style=theme["text"]),
                            dcc.Input(id="scan-hw-f", type="number", value=0.008, step=0.001, style=theme["input"]),
                            html.Label("Temperature (K):", style=theme["text"]),
                            dcc.Input(id="scan-temp", type="number", value=300, step=10, style=theme["input"]),

                            html.Button(
                                "Run Scan",
                                id="run-scan-button",
                                n_clicks=0,
                                style=theme["button"],
                            ),
                            html.Div(id="scan-progress", style={**theme["text"], "marginTop": "10px"}),
                        ],
                    ),
                    # Plot area
                    html.Div(
                        style={**theme["card"], "flex": "1"},
                        children=[
                            html.H2("Parameter Scan Results", style=theme["subheader"]),
                            dcc.RadioItems(
                                id="scan-plot-type",
                                options=[
                                    {"label": " Heatmap", "value": "heatmap"},
                                    {"label": " Contour", "value": "contour"},
                                ],
                                value="heatmap",
                                inline=True,
                                style={**theme["text"], "marginBottom": "10px"},
                            ),
                            dcc.Checklist(
                                id="scan-plot-options",
                                options=[
                                    {"label": " Log Scale", "value": "log"},
                                ],
                                value=["log"],
                                inline=True,
                                style=theme["text"],
                            ),
                            dcc.Graph(
                                id="scan-plot",
                                style={"height": "600px"},
                            ),
                            html.Button(
                                "Export Results",
                                id="export-scan-button",
                                n_clicks=0,
                                style=theme["button"],
                            ),
                            dcc.Download(id="download-scan-results"),
                        ],
                    ),
                ],
            ),
        ],
    )


def create_comparison_tab(theme: Dict[str, Any]) -> html.Div:
    """Create the multi-potential comparison tab layout."""
    return html.Div(
        style={"marginTop": "20px"},
        children=[
            html.Div(
                style={"display": "flex", "gap": "20px"},
                children=[
                    # Controls
                    html.Div(
                        style={**theme["card"], "width": "350px", "height": "fit-content"},
                        children=[
                            html.H2("Comparison", style=theme["subheader"]),
                            html.P("Load multiple potentials to compare", style=theme["text"]),

                            dcc.Upload(
                                id="upload-comparison-potential",
                                children=html.Div([
                                    "Add Potential",
                                ]),
                                style={
                                    "width": "100%",
                                    "height": "60px",
                                    "lineHeight": "60px",
                                    "borderWidth": "1px",
                                    "borderStyle": "dashed",
                                    "borderRadius": "5px",
                                    "textAlign": "center",
                                    "margin": "10px 0",
                                },
                                multiple=False,
                            ),

                            html.Div(id="comparison-list", style={**theme["text"], "marginTop": "10px"}),

                            html.Button(
                                "Clear All",
                                id="clear-comparison-button",
                                n_clicks=0,
                                style=theme["button"],
                            ),
                        ],
                    ),
                    # Plot area
                    html.Div(
                        style={**theme["card"], "flex": "1"},
                        children=[
                            html.H2("Potential Comparison", style=theme["subheader"]),
                            dcc.Graph(
                                id="comparison-plot",
                                style={"height": "600px"},
                            ),
                        ],
                    ),
                ],
            ),
        ],
    )


def create_capture_tab(theme: Dict[str, Any]) -> html.Div:
    """Create the capture calculation tab layout."""
    return html.Div(
        style={"marginTop": "20px"},
        children=[
            html.Div(
                style={"display": "flex", "gap": "20px"},
                children=[
                    # Controls
                    html.Div(
                        style={**theme["card"], "width": "350px", "height": "fit-content"},
                        children=[
                            html.H2("Capture Calculation", style=theme["subheader"]),

                            # Load initial state
                            html.Label("Initial State (Excited):", style=theme["text"]),
                            dcc.Upload(
                                id="upload-initial-state",
                                children=html.Div(["Load Initial"]),
                                style={
                                    "width": "100%",
                                    "height": "50px",
                                    "lineHeight": "50px",
                                    "borderWidth": "1px",
                                    "borderStyle": "dashed",
                                    "borderRadius": "5px",
                                    "textAlign": "center",
                                    "margin": "10px 0",
                                },
                            ),
                            html.Div(id="initial-status", style=theme["text"]),

                            # Load final state
                            html.Label("Final State (Relaxed):", style=theme["text"]),
                            dcc.Upload(
                                id="upload-final-state",
                                children=html.Div(["Load Final"]),
                                style={
                                    "width": "100%",
                                    "height": "50px",
                                    "lineHeight": "50px",
                                    "borderWidth": "1px",
                                    "borderStyle": "dashed",
                                    "borderRadius": "5px",
                                    "textAlign": "center",
                                    "margin": "10px 0",
                                },
                            ),
                            html.Div(id="final-status", style=theme["text"]),

                            html.Hr(),

                            # Calculation parameters
                            html.H3("Parameters", style=theme["subheader"]),
                            html.Label("W (e-ph coupling, eV):", style=theme["text"]),
                            dcc.Input(id="capture-w", type="number", value=0.068, step=0.001, style=theme["input"]),
                            html.Label("Q₀ (crossing point):", style=theme["text"]),
                            dcc.Input(id="capture-q0", type="number", value=10.0, step=0.1, style=theme["input"]),
                            html.Label("Volume (cm³):", style=theme["text"]),
                            dcc.Input(id="capture-volume", type="number", value=1e-21, step=1e-22, style=theme["input"]),
                            html.Label("Degeneracy:", style=theme["text"]),
                            dcc.Input(id="capture-degeneracy", type="number", value=1, min=1, style=theme["input"]),

                            html.Hr(),

                            # Temperature range
                            html.H3("Temperature Range", style=theme["subheader"]),
                            html.Label("Min (K):", style=theme["text"]),
                            dcc.Input(id="temp-min", type="number", value=100, style=theme["input"]),
                            html.Label("Max (K):", style=theme["text"]),
                            dcc.Input(id="temp-max", type="number", value=500, style=theme["input"]),
                            html.Label("Points:", style=theme["text"]),
                            dcc.Input(id="temp-points", type="number", value=50, min=2, style=theme["input"]),

                            html.Button(
                                "Calculate",
                                id="calculate-capture-button",
                                n_clicks=0,
                                style=theme["button"],
                            ),
                        ],
                    ),
                    # Plot area
                    html.Div(
                        style={**theme["card"], "flex": "1"},
                        children=[
                            html.H2("Capture Coefficient", style=theme["subheader"]),
                            dcc.Graph(
                                id="capture-plot",
                                style={"height": "600px"},
                            ),
                            html.Button(
                                "Export Results",
                                id="export-capture-button",
                                n_clicks=0,
                                style=theme["button"],
                            ),
                            dcc.Download(id="download-capture-results"),
                        ],
                    ),
                ],
            ),
        ],
    )


def register_callbacks(app: dash.Dash) -> None:
    """Register all Dash callbacks."""

    # Potential tab callbacks
    register_potential_callbacks(app)

    # Scan tab callbacks
    register_scan_callbacks(app)

    # Comparison tab callbacks
    register_comparison_callbacks(app)

    # Capture tab callbacks
    register_capture_callbacks(app)


def register_potential_callbacks(app: dash.Dash) -> None:
    """Register callbacks for potential viewer tab."""

    @app.callback(
        [Output("potential-store-initial", "data"),
         Output("upload-status", "children"),
         Output("potential-plot", "figure"),
         Output("status-bar", "children")],
        [Input("upload-potential-data", "contents"),
         Input("generate-harmonic-button", "n_clicks"),
         Input("fit-button", "n_clicks"),
         Input("solve-button", "n_clicks")],
        [State("upload-potential-data", "filename"),
         State("harmonic-hw", "value"),
         State("harmonic-q0", "value"),
         State("harmonic-e0", "value"),
         State("fit-type", "value"),
         State("spline-order", "value"),
         State("smoothness", "value"),
         State("nev", "value"),
         State("display-options", "value"),
         State("wf-scale", "value"),
         State("potential-store-initial", "data")],
        prevent_initial_call=True,
    )
    def handle_potential_operations(upload_contents, gen_clicks, fit_clicks, solve_clicks,
                                    filename, hw, q0, e0, fit_type, spline_order, smoothness,
                                    nev, display_options, wf_scale, current_pot_data):
        """Handle all potential operations."""
        triggered_id = ctx.triggered_id

        try:
            # Load or generate potential
            pot = None
            status_msg = ""

            if triggered_id == "upload-potential-data" and upload_contents:
                # Parse uploaded file
                content_type, content_string = upload_contents.split(',')
                decoded = base64.b64decode(content_string)

                # Try to parse as CSV/DAT
                try:
                    data = np.loadtxt(io.BytesIO(decoded))
                    Q_data = data[:, 0]
                    E_data = data[:, 1]

                    pot = Potential(name=filename or "Uploaded")
                    pot.Q_data = Q_data
                    pot.E_data = E_data

                    # Sort and set Q0, E0
                    sort_idx = np.argsort(Q_data)
                    pot.Q_data = Q_data[sort_idx]
                    pot.E_data = E_data[sort_idx]
                    min_idx = np.argmin(E_data)
                    pot.Q0 = Q_data[min_idx]
                    pot.E0 = E_data[min_idx]

                    pot.Q = np.linspace(Q_data.min(), Q_data.max(), 3001)

                    status_msg = f"✓ Loaded {filename}: {len(Q_data)} data points"

                except Exception as e:
                    return None, f"✗ Error parsing file: {e}", go.Figure(), "Error loading file"

            elif triggered_id == "generate-harmonic-button":
                pot = Potential.from_harmonic(
                    hw=hw,
                    Q0=q0,
                    E0=e0,
                    Q_range=(-20 - abs(q0), 20 + abs(q0)),
                    npoints=3001,
                )
                pot.Q_data = pot.Q[::100]  # Subsample for display
                pot.E_data = pot.E[::100]
                status_msg = f"✓ Generated harmonic potential: ℏω={hw} eV, Q₀={q0}"

            elif triggered_id == "fit-button" and current_pot_data:
                # Fit existing potential
                pot = deserialize_potential(current_pot_data)

                if fit_type == "spline":
                    s = 10 ** smoothness
                    pot.fit(fit_type="spline", order=spline_order, smoothness=s)
                    status_msg = f"✓ Fitted with spline (order={spline_order}, s={s:.2e})"
                elif fit_type == "harmonic":
                    pot.fit(fit_type="harmonic")
                    status_msg = "✓ Fitted with harmonic potential"
                elif fit_type == "morse":
                    pot.fit(fit_type="morse")
                    status_msg = "✓ Fitted with Morse potential"
                else:
                    pot.fit(fit_type="polynomial", order=spline_order)
                    status_msg = f"✓ Fitted with polynomial (order={spline_order})"

            elif triggered_id == "solve-button" and current_pot_data:
                # Solve existing potential
                pot = deserialize_potential(current_pot_data)
                pot.solve(nev=nev)
                status_msg = f"✓ Solved Schrödinger equation: {nev} eigenvalues"

            elif current_pot_data:
                # Just update display options
                pot = deserialize_potential(current_pot_data)
                status_msg = "Display updated"

            if pot is None:
                return None, "No potential loaded", go.Figure(), ""

            # Create plot
            fig = create_potential_figure(pot, display_options, wf_scale)

            # Serialize potential
            pot_data = serialize_potential(pot)

            return pot_data, status_msg, fig, status_msg

        except Exception as e:
            return None, f"✗ Error: {e}", go.Figure(), f"Error: {e}"


def register_scan_callbacks(app: dash.Dash) -> None:
    """Register callbacks for parameter scan tab."""

    @app.callback(
        [Output("scan-results-store", "data"),
         Output("scan-plot", "figure"),
         Output("scan-progress", "children")],
        [Input("upload-scan-results", "contents"),
         Input("run-scan-button", "n_clicks"),
         Input("scan-plot-type", "value"),
         Input("scan-plot-options", "value")],
        [State("upload-scan-results", "filename"),
         State("scan-dq-min", "value"),
         State("scan-dq-max", "value"),
         State("scan-dq-points", "value"),
         State("scan-de-min", "value"),
         State("scan-de-max", "value"),
         State("scan-de-points", "value"),
         State("scan-hw-i", "value"),
         State("scan-hw-f", "value"),
         State("scan-temp", "value"),
         State("scan-results-store", "data")],
        prevent_initial_call=True,
    )
    def handle_scan_operations(upload_contents, run_clicks, plot_type, plot_options,
                               filename, dq_min, dq_max, dq_points, de_min, de_max, de_points,
                               hw_i, hw_f, temp, current_results):
        """Handle parameter scan operations."""
        triggered_id = ctx.triggered_id

        try:
            results = None
            status_msg = ""

            if triggered_id == "upload-scan-results" and upload_contents:
                # Load scan results
                content_type, content_string = upload_contents.split(',')
                decoded = base64.b64decode(content_string)

                # Save to temporary file and load
                import tempfile
                with tempfile.NamedTemporaryFile(suffix=".npz", delete=False) as tmp:
                    tmp.write(decoded)
                    tmp_path = Path(tmp.name)

                results = ScanResult.load(tmp_path, format="npz")
                tmp_path.unlink()

                status_msg = f"✓ Loaded scan results: {results.dQ_grid.shape[0]}×{results.dE_grid.shape[0]} grid"

            elif triggered_id == "run-scan-button":
                # Run new scan
                from carriercapture.analysis.parameter_scan import ParameterScanner, ScanParameters

                params = ScanParameters(
                    dQ_range=(dq_min, dq_max, dq_points),
                    dE_range=(de_min, de_max, de_points),
                    hbar_omega_i=hw_i,
                    hbar_omega_f=hw_f,
                    temperature=temp,
                )

                scanner = ParameterScanner(params, verbose=False)
                results = scanner.run_harmonic_scan(n_jobs=1, show_progress=False)

                n_success = np.sum(~np.isnan(results.capture_coefficients))
                status_msg = f"✓ Scan complete: {n_success}/{results.capture_coefficients.size} successful"

            elif current_results:
                # Update plot only
                results = deserialize_scan_results(current_results)
                status_msg = "Plot updated"

            if results is None:
                return None, go.Figure(), "No scan results"

            # Create plot
            fig = create_scan_figure(results, plot_type, "log" in plot_options)

            # Serialize results
            results_data = serialize_scan_results(results)

            return results_data, fig, status_msg

        except Exception as e:
            return None, go.Figure(), f"✗ Error: {e}"


def register_comparison_callbacks(app: dash.Dash) -> None:
    """Register callbacks for comparison tab."""

    @app.callback(
        [Output("comparison-potentials-store", "data"),
         Output("comparison-list", "children"),
         Output("comparison-plot", "figure")],
        [Input("upload-comparison-potential", "contents"),
         Input("clear-comparison-button", "n_clicks")],
        [State("upload-comparison-potential", "filename"),
         State("comparison-potentials-store", "data")],
        prevent_initial_call=True,
    )
    def handle_comparison(upload_contents, clear_clicks, filename, current_potentials):
        """Handle comparison operations."""
        triggered_id = ctx.triggered_id

        potentials = current_potentials or []

        if triggered_id == "clear-comparison-button":
            potentials = []

        elif triggered_id == "upload-comparison-potential" and upload_contents:
            # Parse and add new potential
            try:
                content_type, content_string = upload_contents.split(',')
                decoded = base64.b64decode(content_string)
                data = np.loadtxt(io.BytesIO(decoded))

                pot = Potential(name=filename or f"Potential {len(potentials)+1}")
                pot.Q_data = data[:, 0]
                pot.E_data = data[:, 1]
                pot.Q = np.linspace(pot.Q_data.min(), pot.Q_data.max(), 1000)

                potentials.append(serialize_potential(pot))
            except Exception as e:
                pass

        # Create list display
        list_items = [
            html.Div(f"{i+1}. {deserialize_potential(p).name}", style={"marginBottom": "5px"})
            for i, p in enumerate(potentials)
        ]

        # Create comparison plot
        fig = create_comparison_figure([deserialize_potential(p) for p in potentials])

        return potentials, list_items, fig


def register_capture_callbacks(app: dash.Dash) -> None:
    """Register callbacks for capture calculation tab."""

    @app.callback(
        [Output("potential-store-initial", "data", allow_duplicate=True),
         Output("potential-store-final", "data"),
         Output("initial-status", "children"),
         Output("final-status", "children"),
         Output("capture-plot", "figure"),
         Output("status-bar", "children", allow_duplicate=True)],
        [Input("upload-initial-state", "contents"),
         Input("upload-final-state", "contents"),
         Input("calculate-capture-button", "n_clicks")],
        [State("upload-initial-state", "filename"),
         State("upload-final-state", "filename"),
         State("capture-w", "value"),
         State("capture-q0", "value"),
         State("capture-volume", "value"),
         State("capture-degeneracy", "value"),
         State("temp-min", "value"),
         State("temp-max", "value"),
         State("temp-points", "value"),
         State("potential-store-initial", "data"),
         State("potential-store-final", "data")],
        prevent_initial_call=True,
    )
    def handle_capture_calculation(initial_contents, final_contents, calc_clicks,
                                   initial_filename, final_filename,
                                   W, Q0, volume, degeneracy,
                                   temp_min, temp_max, temp_points,
                                   current_initial, current_final):
        """Handle capture coefficient calculation."""
        triggered_id = ctx.triggered_id

        initial_pot = None
        final_pot = None
        initial_msg = ""
        final_msg = ""
        fig = go.Figure()
        status = ""

        try:
            # Load potentials
            if triggered_id == "upload-initial-state" and initial_contents:
                content_type, content_string = initial_contents.split(',')
                decoded = base64.b64decode(content_string)
                data = np.loadtxt(io.BytesIO(decoded))

                initial_pot = Potential(name=initial_filename or "Initial")
                initial_pot.Q_data = data[:, 0]
                initial_pot.E_data = data[:, 1]
                initial_pot.Q = np.linspace(data[:, 0].min(), data[:, 0].max(), 3001)
                initial_pot.fit(fit_type="spline")
                initial_pot.solve(nev=180)

                initial_msg = f"✓ Loaded initial state"
                final_pot = deserialize_potential(current_final) if current_final else None

            elif triggered_id == "upload-final-state" and final_contents:
                content_type, content_string = final_contents.split(',')
                decoded = base64.b64decode(content_string)
                data = np.loadtxt(io.BytesIO(decoded))

                final_pot = Potential(name=final_filename or "Final")
                final_pot.Q_data = data[:, 0]
                final_pot.E_data = data[:, 1]
                final_pot.Q = np.linspace(data[:, 0].min(), data[:, 0].max(), 3001)
                final_pot.fit(fit_type="spline")
                final_pot.solve(nev=60)

                final_msg = f"✓ Loaded final state"
                initial_pot = deserialize_potential(current_initial) if current_initial else None

            elif triggered_id == "calculate-capture-button":
                # Calculate capture coefficient
                if not current_initial or not current_final:
                    status = "Error: Load both initial and final states"
                    return current_initial, current_final, initial_msg, final_msg, fig, status

                initial_pot = deserialize_potential(current_initial)
                final_pot = deserialize_potential(current_final)

                # Ensure potentials are solved
                if initial_pot.eigenvalues is None:
                    initial_pot.solve(nev=180)
                if final_pot.eigenvalues is None:
                    final_pot.solve(nev=60)

                # Create ConfigCoordinate
                cc = ConfigCoordinate(
                    pot_i=initial_pot,
                    pot_f=final_pot,
                    W=W,
                    degeneracy=degeneracy,
                )

                # Calculate overlap
                cc.calculate_overlap(Q0=Q0, cutoff=0.25, sigma=0.01)

                # Calculate capture coefficient
                temperature = np.linspace(temp_min, temp_max, temp_points)
                cc.calculate_capture_coefficient(volume=volume, temperature=temperature)

                # Create plot
                fig = plot_capture_coefficient(cc)

                status = "✓ Capture coefficient calculated"
                initial_msg = "Initial state loaded"
                final_msg = "Final state loaded"

            # Serialize potentials
            initial_data = serialize_potential(initial_pot) if initial_pot else current_initial
            final_data = serialize_potential(final_pot) if final_pot else current_final

            return initial_data, final_data, initial_msg, final_msg, fig, status

        except Exception as e:
            return current_initial, current_final, initial_msg, final_msg, fig, f"Error: {e}"


# Helper functions for serialization
def serialize_potential(pot: Potential) -> Dict[str, Any]:
    """Serialize Potential to JSON-compatible dict."""
    data = {
        "name": pot.name,
        "Q_data": pot.Q_data.tolist() if pot.Q_data is not None else None,
        "E_data": pot.E_data.tolist() if pot.E_data is not None else None,
        "Q": pot.Q.tolist() if pot.Q is not None else None,
        "E": pot.E.tolist() if pot.E is not None else None,
        "Q0": float(pot.Q0) if pot.Q0 is not None else None,
        "E0": float(pot.E0) if pot.E0 is not None else None,
        "fit_type": pot.fit_type,
        "eigenvalues": pot.eigenvalues.tolist() if pot.eigenvalues is not None else None,
        "eigenvectors": pot.eigenvectors.tolist() if pot.eigenvectors is not None else None,
    }
    return data


def deserialize_potential(data: Dict[str, Any]) -> Potential:
    """Deserialize Potential from dict."""
    pot = Potential(name=data["name"])
    pot.Q_data = np.array(data["Q_data"]) if data["Q_data"] else None
    pot.E_data = np.array(data["E_data"]) if data["E_data"] else None
    pot.Q = np.array(data["Q"]) if data["Q"] else None
    pot.E = np.array(data["E"]) if data["E"] else None
    pot.Q0 = data["Q0"]
    pot.E0 = data["E0"]
    pot.fit_type = data["fit_type"]
    pot.eigenvalues = np.array(data["eigenvalues"]) if data["eigenvalues"] else None
    pot.eigenvectors = np.array(data["eigenvectors"]) if data["eigenvectors"] else None
    return pot


def serialize_scan_results(results: ScanResult) -> Dict[str, Any]:
    """Serialize ScanResult to dict."""
    return {
        "dQ_grid": results.dQ_grid.tolist(),
        "dE_grid": results.dE_grid.tolist(),
        "capture_coefficients": results.capture_coefficients.tolist(),
        "barrier_heights": results.barrier_heights.tolist(),
    }


def deserialize_scan_results(data: Dict[str, Any]) -> ScanResult:
    """Deserialize ScanResult from dict."""
    from carriercapture.analysis.parameter_scan import ScanParameters

    # Create dummy parameters
    params = ScanParameters(
        dQ_range=(0, 1, 2),
        dE_range=(0, 1, 2),
    )

    return ScanResult(
        dQ_grid=np.array(data["dQ_grid"]),
        dE_grid=np.array(data["dE_grid"]),
        capture_coefficients=np.array(data["capture_coefficients"]),
        barrier_heights=np.array(data["barrier_heights"]),
        parameters=params,
        metadata={},
    )


def create_potential_figure(pot: Potential, display_options: List[str], wf_scale: float) -> go.Figure:
    """Create figure for potential plot."""
    fig = go.Figure()

    # Plot fitted curve
    if pot.E is not None:
        fig.add_trace(go.Scatter(
            x=pot.Q,
            y=pot.E,
            mode="lines",
            name="Potential",
            line=dict(color=COLORS["primary"], width=2),
        ))

    # Plot data points
    if "data" in display_options and pot.Q_data is not None:
        fig.add_trace(go.Scatter(
            x=pot.Q_data,
            y=pot.E_data,
            mode="markers",
            name="Data",
            marker=dict(color=COLORS["secondary"], size=8),
        ))

    # Plot eigenvalues
    if "ev" in display_options and pot.eigenvalues is not None:
        for i, ev in enumerate(pot.eigenvalues[:20]):  # First 20
            fig.add_hline(
                y=ev,
                line=dict(color=COLORS["info"], width=1, dash="dash"),
                opacity=0.5,
            )

    # Plot wavefunctions
    if "wf" in display_options and pot.eigenvectors is not None:
        for i in range(min(5, len(pot.eigenvalues))):  # First 5
            psi = pot.eigenvectors[i, :] * wf_scale
            fig.add_trace(go.Scatter(
                x=pot.Q,
                y=pot.eigenvalues[i] + psi,
                mode="lines",
                name=f"ψ_{i}",
                line=dict(width=1),
                opacity=0.7,
            ))

    fig.update_layout(
        title=f"Potential: {pot.name}",
        xaxis_title="Q (amu^0.5·Å)",
        yaxis_title="Energy (eV)",
        template="plotly_white",
        showlegend=True,
    )

    return fig


def create_scan_figure(results: ScanResult, plot_type: str, log_scale: bool) -> go.Figure:
    """Create figure for scan results."""
    Z = results.capture_coefficients.copy()

    if log_scale:
        Z = np.log10(Z + 1e-30)
        colorbar_title = "log₁₀(C) [cm³/s]"
    else:
        colorbar_title = "C [cm³/s]"

    fig = go.Figure()

    if plot_type == "heatmap":
        fig.add_trace(go.Heatmap(
            x=results.dE_grid,
            y=results.dQ_grid,
            z=Z,
            colorscale="Viridis",
            colorbar=dict(title=colorbar_title),
        ))
    else:  # contour
        fig.add_trace(go.Contour(
            x=results.dE_grid,
            y=results.dQ_grid,
            z=Z,
            colorscale="Viridis",
            colorbar=dict(title=colorbar_title),
            contours=dict(showlabels=True),
        ))

    fig.update_layout(
        title="Parameter Scan: Capture Coefficient",
        xaxis_title="ΔE (eV)",
        yaxis_title="ΔQ (amu^0.5·Å)",
        template="plotly_white",
    )

    return fig


def create_comparison_figure(potentials: List[Potential]) -> go.Figure:
    """Create figure for potential comparison."""
    fig = go.Figure()

    colors = [COLORS["primary"], COLORS["secondary"], COLORS["info"], COLORS["success"]]

    for i, pot in enumerate(potentials):
        if pot.E is not None:
            fig.add_trace(go.Scatter(
                x=pot.Q,
                y=pot.E,
                mode="lines",
                name=pot.name,
                line=dict(color=colors[i % len(colors)], width=2),
            ))

    fig.update_layout(
        title="Potential Comparison",
        xaxis_title="Q (amu^0.5·Å)",
        yaxis_title="Energy (eV)",
        template="plotly_white",
        showlegend=True,
    )

    return fig


def run_server(port: int = 8050, debug: bool = False, host: str = "127.0.0.1") -> None:
    """
    Run the Dash server.

    Parameters
    ----------
    port : int, default=8050
        Port to run server on
    debug : bool, default=False
        Whether to run in debug mode
    host : str, default="127.0.0.1"
        Host address

    Examples
    --------
    >>> run_server(port=8050)
    """
    app = create_app(debug=debug)

    print(f"\n{'=' * 60}")
    print("CarrierCapture Visualization Dashboard")
    print(f"{'=' * 60}")
    print(f"\n🌐 Server starting at: http://{host}:{port}")
    print(f"\n📊 Features:")
    print(f"  • Potential Viewer - Interactive fitting and solving")
    print(f"  • Parameter Scan - High-throughput screening visualization")
    print(f"  • Comparison - Multi-potential overlay")
    print(f"  • Capture Calculation - Full carrier capture workflow")
    print(f"\nPress Ctrl+C to stop the server\n")

    app.run(host=host, port=port, debug=debug)


__all__ = [
    "create_app",
    "run_server",
]
