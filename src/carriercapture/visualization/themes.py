"""
Theming and styling utilities for CarrierCapture visualizations.

Provides consistent color schemes, layouts, and styling for publication-quality figures.
"""

from typing import Dict, Any
import numpy as np
import plotly.graph_objects as go


# Color palettes
COLORS = {
    "primary": "#1f77b4",  # Blue
    "secondary": "#ff7f0e",  # Orange
    "success": "#2ca02c",  # Green
    "danger": "#d62728",  # Red
    "warning": "#ff9896",  # Light red
    "info": "#17becf",  # Cyan
    "dark": "#2f4f4f",  # Dark slate gray
    "light": "#f0f0f0",  # Light gray
}

# Color scheme for multiple potentials
POTENTIAL_COLORS = [
    "#1f77b4",  # Blue
    "#ff7f0e",  # Orange
    "#2ca02c",  # Green
    "#d62728",  # Red
    "#9467bd",  # Purple
    "#8c564b",  # Brown
    "#e377c2",  # Pink
    "#7f7f7f",  # Gray
    "#bcbd22",  # Olive
    "#17becf",  # Cyan
]


def get_default_layout(
    title: str = "",
    xaxis_title: str = "",
    yaxis_title: str = "",
    width: int = 900,
    height: int = 600,
    **kwargs
) -> Dict[str, Any]:
    """
    Get default layout configuration for CarrierCapture plots.

    Parameters
    ----------
    title : str
        Plot title
    xaxis_title : str
        X-axis label
    yaxis_title : str
        Y-axis label
    width : int, default=900
        Figure width in pixels
    height : int, default=600
        Figure height in pixels
    **kwargs
        Additional layout parameters to override defaults

    Returns
    -------
    dict
        Layout configuration dictionary

    Examples
    --------
    >>> layout = get_default_layout(title="My Plot", xaxis_title="X", yaxis_title="Y")
    >>> fig = go.Figure(layout=layout)
    """
    layout = dict(
        title=dict(
            text=title,
            font=dict(size=18, family="Arial, sans-serif"),
            x=0.5,
            xanchor="center",
        ),
        xaxis=dict(
            title=dict(
                text=xaxis_title,
                font=dict(size=14, family="Arial, sans-serif"),
            ),
            showgrid=True,
            gridwidth=1,
            gridcolor="lightgray",
            showline=True,
            linewidth=2,
            linecolor="black",
            mirror=True,
        ),
        yaxis=dict(
            title=dict(
                text=yaxis_title,
                font=dict(size=14, family="Arial, sans-serif"),
            ),
            showgrid=True,
            gridwidth=1,
            gridcolor="lightgray",
            showline=True,
            linewidth=2,
            linecolor="black",
            mirror=True,
        ),
        font=dict(size=12, family="Arial, sans-serif"),
        plot_bgcolor="white",
        paper_bgcolor="white",
        width=width,
        height=height,
        hovermode="closest",
        template="plotly_white",
        legend=dict(
            orientation="v",
            yanchor="top",
            y=0.99,
            xanchor="right",
            x=0.99,
            bgcolor="rgba(255, 255, 255, 0.8)",
            bordercolor="gray",
            borderwidth=1,
        ),
    )

    # Override with any custom parameters
    layout.update(kwargs)

    return layout


def apply_publication_style(fig: go.Figure) -> go.Figure:
    """
    Apply publication-quality styling to a Plotly figure.

    Parameters
    ----------
    fig : go.Figure
        Plotly figure to style

    Returns
    -------
    go.Figure
        Styled figure

    Examples
    --------
    >>> fig = go.Figure()
    >>> fig = apply_publication_style(fig)
    """
    fig.update_layout(
        font=dict(size=14, family="Arial, sans-serif", color="black"),
        plot_bgcolor="white",
        paper_bgcolor="white",
        xaxis=dict(
            showgrid=True,
            gridwidth=1,
            gridcolor="lightgray",
            showline=True,
            linewidth=2,
            linecolor="black",
            mirror=True,
            ticks="outside",
            tickwidth=2,
            tickcolor="black",
        ),
        yaxis=dict(
            showgrid=True,
            gridwidth=1,
            gridcolor="lightgray",
            showline=True,
            linewidth=2,
            linecolor="black",
            mirror=True,
            ticks="outside",
            tickwidth=2,
            tickcolor="black",
        ),
        legend=dict(
            bgcolor="rgba(255, 255, 255, 0.9)",
            bordercolor="black",
            borderwidth=1,
        ),
    )

    return fig


def get_colorscale(name: str = "viridis") -> list:
    """
    Get a colorscale for heatmaps and contour plots.

    Parameters
    ----------
    name : str, default="viridis"
        Colorscale name: "viridis", "plasma", "inferno", "magma", "cividis"

    Returns
    -------
    list
        Colorscale specification

    Examples
    --------
    >>> colorscale = get_colorscale("plasma")
    """
    colorscales = {
        "viridis": "Viridis",
        "plasma": "Plasma",
        "inferno": "Inferno",
        "magma": "Magma",
        "cividis": "Cividis",
        "blues": "Blues",
        "reds": "Reds",
        "greens": "Greens",
    }

    return colorscales.get(name.lower(), "Viridis")


def format_scientific(value: float, precision: int = 2) -> str:
    """
    Format a number in scientific notation for display.

    Parameters
    ----------
    value : float
        Number to format
    precision : int, default=2
        Number of decimal places

    Returns
    -------
    str
        Formatted string

    Examples
    --------
    >>> format_scientific(1.23e-10, precision=2)
    '1.23×10⁻¹⁰'
    """
    if value == 0:
        return "0"

    exponent = int(np.floor(np.log10(abs(value))))
    mantissa = value / (10 ** exponent)

    # Unicode superscript digits
    superscripts = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")
    exp_str = str(exponent).translate(superscripts)

    return f"{mantissa:.{precision}f}×10{exp_str}"


def create_dash_theme() -> Dict[str, Any]:
    """
    Create theme configuration for Dash applications.

    Returns
    -------
    dict
        Dash theme configuration

    Examples
    --------
    >>> theme = create_dash_theme()
    >>> app.layout = html.Div(style=theme["container"])
    """
    return {
        "container": {
            "fontFamily": "Arial, sans-serif",
            "backgroundColor": "#f8f9fa",
            "padding": "20px",
        },
        "card": {
            "backgroundColor": "white",
            "border": "1px solid #dee2e6",
            "borderRadius": "5px",
            "padding": "20px",
            "marginBottom": "20px",
            "boxShadow": "0 2px 4px rgba(0,0,0,0.1)",
        },
        "header": {
            "fontSize": "24px",
            "fontWeight": "bold",
            "color": "#212529",
            "marginBottom": "10px",
        },
        "subheader": {
            "fontSize": "18px",
            "fontWeight": "bold",
            "color": "#495057",
            "marginBottom": "8px",
        },
        "text": {
            "fontSize": "14px",
            "color": "#6c757d",
        },
        "button": {
            "backgroundColor": COLORS["primary"],
            "color": "white",
            "border": "none",
            "borderRadius": "4px",
            "padding": "10px 20px",
            "fontSize": "14px",
            "cursor": "pointer",
            "marginRight": "10px",
        },
        "input": {
            "border": "1px solid #ced4da",
            "borderRadius": "4px",
            "padding": "8px 12px",
            "fontSize": "14px",
            "width": "100%",
        },
        "slider": {
            "marginBottom": "20px",
        },
    }


__all__ = [
    "COLORS",
    "POTENTIAL_COLORS",
    "get_default_layout",
    "apply_publication_style",
    "get_colorscale",
    "format_scientific",
    "create_dash_theme",
]
