"""
Static plotting functions using Plotly.

Provides publication-quality plots for:
- Potential energy surfaces with wavefunctions
- Capture coefficient (Arrhenius plots)
- Eigenvalue spectra
- Configuration coordinate diagrams
"""

from typing import Optional, List, Union, Tuple
import numpy as np
from numpy.typing import NDArray
import plotly.graph_objects as go
from plotly.subplots import make_subplots


def plot_potential(
    potential,
    show_wavefunctions: bool = True,
    wf_scale: float = 0.02,
    wf_sampling: int = 5,
    show_data: bool = True,
    color: str = "blue",
    title: Optional[str] = None,
    fig: Optional[go.Figure] = None,
    show_eigenvalues: bool = True,
    max_wf_to_plot: Optional[int] = None,
) -> go.Figure:
    """
    Plot potential energy surface with optional wavefunctions.

    Parameters
    ----------
    potential : Potential
        Potential object with Q, E, and optionally eigenvalues/eigenvectors
    show_wavefunctions : bool, default=True
        Whether to plot wavefunction envelopes
    wf_scale : float, default=0.02
        Scaling factor for wavefunction amplitude
    wf_sampling : int, default=5
        Plot every Nth wavefunction
    show_data : bool, default=True
        Whether to show original data points
    color : str, default="blue"
        Color for the potential curve
    title : str, optional
        Plot title (defaults to potential name)
    fig : go.Figure, optional
        Existing figure to add to
    show_eigenvalues : bool, default=True
        Whether to show eigenvalue lines
    max_wf_to_plot : int, optional
        Maximum number of wavefunctions to plot

    Returns
    -------
    go.Figure
        Plotly figure object

    Examples
    --------
    >>> fig = plot_potential(pot, show_wavefunctions=True)
    >>> fig.show()
    >>> fig.write_html("potential.html")
    """
    if fig is None:
        fig = go.Figure()

    # Get title
    if title is None:
        title = potential.name if potential.name else "Potential Energy Surface"

    # Plot original data points if available
    if show_data and potential.Q_data is not None and potential.E_data is not None:
        fig.add_trace(
            go.Scatter(
                x=potential.Q_data,
                y=potential.E_data,
                mode="markers",
                marker=dict(size=6, color=color, opacity=0.6),
                name="Data",
                showlegend=True,
            )
        )

    # Plot fitted potential curve
    if potential.Q is not None and potential.E is not None:
        fig.add_trace(
            go.Scatter(
                x=potential.Q,
                y=potential.E,
                mode="lines",
                line=dict(color=color, width=3),
                name=potential.name or "Potential",
                showlegend=True,
            )
        )

    # Plot wavefunctions as filled envelopes
    if show_wavefunctions and potential.eigenvectors is not None:
        n_wf = len(potential.eigenvalues)
        if max_wf_to_plot:
            n_wf = min(n_wf, max_wf_to_plot)

        # Determine which wavefunctions to plot
        wf_indices = range(0, n_wf, wf_sampling)

        for i in wf_indices:
            psi = potential.eigenvectors[i, :]
            E_i = potential.eigenvalues[i]

            # Scale wavefunction
            E_range = potential.E.max() - potential.E.min()
            psi_scaled = psi * wf_scale * E_range

            # Upper and lower envelopes
            upper = E_i + psi_scaled
            lower = E_i - psi_scaled

            # Plot filled area
            fig.add_trace(
                go.Scatter(
                    x=potential.Q,
                    y=upper,
                    mode="lines",
                    line=dict(width=0),
                    showlegend=False,
                    hoverinfo="skip",
                )
            )
            fig.add_trace(
                go.Scatter(
                    x=potential.Q,
                    y=lower,
                    mode="lines",
                    line=dict(width=0),
                    fillcolor=f"rgba(100, 100, 100, 0.3)",
                    fill="tonexty",
                    showlegend=False,
                    hoverinfo="skip",
                )
            )

    # Plot eigenvalue lines
    if show_eigenvalues and potential.eigenvalues is not None:
        Q_min, Q_max = potential.Q.min(), potential.Q.max()

        # Plot first few eigenvalues as dashed lines
        n_to_label = min(5, len(potential.eigenvalues))
        for i in range(n_to_label):
            E_i = potential.eigenvalues[i]
            fig.add_trace(
                go.Scatter(
                    x=[Q_min, Q_max],
                    y=[E_i, E_i],
                    mode="lines",
                    line=dict(color="gray", width=1, dash="dash"),
                    name=f"n={i}",
                    showlegend=(i < 3),
                    hovertemplate=f"n={i}<br>E={E_i:.4f} eV<extra></extra>",
                )
            )

        # Plot remaining eigenvalues without labels
        for i in range(n_to_label, len(potential.eigenvalues)):
            E_i = potential.eigenvalues[i]
            fig.add_trace(
                go.Scatter(
                    x=[Q_min, Q_max],
                    y=[E_i, E_i],
                    mode="lines",
                    line=dict(color="lightgray", width=0.5, dash="dot"),
                    showlegend=False,
                    hovertemplate=f"n={i}<br>E={E_i:.4f} eV<extra></extra>",
                )
            )

    # Update layout
    fig.update_layout(
        title=title,
        xaxis_title="Q (amu<sup>0.5</sup>·Å)",
        yaxis_title="Energy (eV)",
        template="plotly_white",
        font=dict(size=14),
        hovermode="closest",
        width=900,
        height=600,
    )

    return fig


def plot_capture_coefficient(
    cc,
    color: str = "red",
    title: Optional[str] = None,
    fig: Optional[go.Figure] = None,
    show_temps: bool = True,
) -> go.Figure:
    """
    Plot capture coefficient as Arrhenius plot (log C vs 1000/T).

    Parameters
    ----------
    cc : ConfigCoordinate
        ConfigCoordinate object with capture_coefficient and temperature
    color : str, default="red"
        Color for the line
    title : str, optional
        Plot title
    fig : go.Figure, optional
        Existing figure to add to
    show_temps : bool, default=True
        Whether to show temperature labels on secondary x-axis

    Returns
    -------
    go.Figure
        Plotly figure object

    Examples
    --------
    >>> fig = plot_capture_coefficient(cc)
    >>> fig.show()
    """
    if fig is None:
        fig = go.Figure()

    if title is None:
        title = "Capture Coefficient (Arrhenius Plot)"

    # Calculate 1000/T
    inv_T = 1000.0 / cc.temperature
    log_C = np.log10(cc.capture_coefficient)

    # Main plot: log(C) vs 1000/T
    fig.add_trace(
        go.Scatter(
            x=inv_T,
            y=log_C,
            mode="lines+markers",
            line=dict(color=color, width=3),
            marker=dict(size=6, color=color),
            name="C(T)",
            hovertemplate="T=%{customdata:.0f} K<br>"
            "1000/T=%{x:.2f} K<sup>-1</sup><br>"
            "log₁₀(C)=%{y:.2f}<br>"
            "C=%{text}<extra></extra>",
            customdata=cc.temperature,
            text=[f"{c:.2e} cm³/s" for c in cc.capture_coefficient],
        )
    )

    # Update layout
    fig.update_layout(
        title=title,
        xaxis_title="1000/T (K<sup>-1</sup>)",
        yaxis_title="log₁₀(C) [cm³/s]",
        template="plotly_white",
        font=dict(size=14),
        hovermode="closest",
        width=900,
        height=600,
    )

    # Add temperature labels on top axis if requested
    if show_temps:
        # Add secondary x-axis with temperature labels
        temp_labels = [100, 200, 300, 400, 500, 600, 700, 800]
        temp_positions = [1000.0 / t for t in temp_labels if cc.temperature.min() <= t <= cc.temperature.max()]
        temp_labels_filtered = [t for t in temp_labels if cc.temperature.min() <= t <= cc.temperature.max()]

        fig.update_layout(
            xaxis2=dict(
                overlaying="x",
                side="top",
                tickmode="array",
                tickvals=temp_positions,
                ticktext=[f"{t}K" for t in temp_labels_filtered],
                title="Temperature",
            )
        )

    return fig


def plot_eigenvalue_spectrum(
    potential,
    max_levels: int = 30,
    color: str = "blue",
    title: Optional[str] = None,
) -> go.Figure:
    """
    Plot eigenvalue spectrum as energy levels.

    Parameters
    ----------
    potential : Potential
        Potential object with eigenvalues
    max_levels : int, default=30
        Maximum number of levels to plot
    color : str, default="blue"
        Color for the levels
    title : str, optional
        Plot title

    Returns
    -------
    go.Figure
        Plotly figure object

    Examples
    --------
    >>> fig = plot_eigenvalue_spectrum(pot, max_levels=20)
    >>> fig.show()
    """
    if potential.eigenvalues is None:
        raise ValueError("Potential must be solved before plotting eigenvalue spectrum")

    fig = go.Figure()

    if title is None:
        title = f"Eigenvalue Spectrum: {potential.name or 'Potential'}"

    n_levels = min(max_levels, len(potential.eigenvalues))
    eigenvalues = potential.eigenvalues[:n_levels]

    # Plot as horizontal lines
    for i, E_i in enumerate(eigenvalues):
        fig.add_trace(
            go.Scatter(
                x=[i - 0.4, i + 0.4],
                y=[E_i, E_i],
                mode="lines",
                line=dict(color=color, width=3),
                showlegend=False,
                hovertemplate=f"n={i}<br>E={E_i:.6f} eV<extra></extra>",
            )
        )

        # Add text label for first few levels
        if i < 5:
            fig.add_annotation(
                x=i + 0.5,
                y=E_i,
                text=f"n={i}",
                showarrow=False,
                xanchor="left",
                font=dict(size=10),
            )

    # Calculate spacing statistics
    if len(eigenvalues) > 1:
        spacing = np.diff(eigenvalues)
        avg_spacing = spacing.mean()
        min_spacing = spacing.min()
        max_spacing = spacing.max()

        # Add statistics as annotation
        stats_text = (
            f"Levels: {n_levels}<br>"
            f"Avg spacing: {avg_spacing:.6f} eV<br>"
            f"Min/Max: {min_spacing:.6f} / {max_spacing:.6f} eV"
        )

        fig.add_annotation(
            x=0.02,
            y=0.98,
            text=stats_text,
            showarrow=False,
            xref="paper",
            yref="paper",
            xanchor="left",
            yanchor="top",
            bgcolor="rgba(255, 255, 255, 0.8)",
            bordercolor="gray",
            borderwidth=1,
            font=dict(size=10),
        )

    fig.update_layout(
        title=title,
        xaxis_title="Quantum Number n",
        yaxis_title="Energy (eV)",
        template="plotly_white",
        font=dict(size=14),
        width=900,
        height=600,
        xaxis=dict(dtick=5),
    )

    return fig


def plot_configuration_coordinate(
    pot_initial,
    pot_final,
    Q0: Optional[float] = None,
    title: Optional[str] = None,
    show_wavefunctions: bool = False,
) -> go.Figure:
    """
    Plot configuration coordinate diagram with both potentials.

    Parameters
    ----------
    pot_initial : Potential
        Initial state potential
    pot_final : Potential
        Final state potential
    Q0 : float, optional
        Q coordinate for overlap calculation (shown as vertical line)
    title : str, optional
        Plot title
    show_wavefunctions : bool, default=False
        Whether to show wavefunctions

    Returns
    -------
    go.Figure
        Plotly figure object

    Examples
    --------
    >>> fig = plot_configuration_coordinate(pot_i, pot_f, Q0=10.0)
    >>> fig.show()
    """
    fig = go.Figure()

    if title is None:
        title = "Configuration Coordinate Diagram"

    # Plot initial potential
    plot_potential(
        pot_initial,
        show_wavefunctions=show_wavefunctions,
        color="blue",
        fig=fig,
        show_data=False,
    )

    # Plot final potential
    plot_potential(
        pot_final,
        show_wavefunctions=show_wavefunctions,
        color="red",
        fig=fig,
        show_data=False,
    )

    # Add Q0 line if provided
    if Q0 is not None:
        E_min = min(pot_initial.E.min(), pot_final.E.min())
        E_max = max(pot_initial.E.max(), pot_final.E.max())

        fig.add_trace(
            go.Scatter(
                x=[Q0, Q0],
                y=[E_min, E_max],
                mode="lines",
                line=dict(color="green", width=2, dash="dash"),
                name=f"Q₀={Q0:.2f}",
                showlegend=True,
            )
        )

    fig.update_layout(
        title=title,
        xaxis_title="Q (amu<sup>0.5</sup>·Å)",
        yaxis_title="Energy (eV)",
        template="plotly_white",
        font=dict(size=14),
        width=900,
        height=600,
    )

    return fig


def plot_overlap_matrix(
    cc,
    title: Optional[str] = None,
    log_scale: bool = True,
) -> go.Figure:
    """
    Plot wavefunction overlap matrix as heatmap.

    Parameters
    ----------
    cc : ConfigCoordinate
        ConfigCoordinate object with overlap_matrix
    title : str, optional
        Plot title
    log_scale : bool, default=True
        Whether to use log scale for colors

    Returns
    -------
    go.Figure
        Plotly figure object

    Examples
    --------
    >>> fig = plot_overlap_matrix(cc)
    >>> fig.show()
    """
    if cc.overlap_matrix is None:
        raise ValueError("ConfigCoordinate must have overlap_matrix calculated")

    fig = go.Figure()

    if title is None:
        title = "Wavefunction Overlap Matrix"

    # Use absolute values and optionally log scale
    overlap_abs = np.abs(cc.overlap_matrix)

    if log_scale:
        # Avoid log(0) by adding small epsilon
        overlap_plot = np.log10(overlap_abs + 1e-10)
        colorbar_title = "log₁₀|⟨i|Q̂|j⟩|"
    else:
        overlap_plot = overlap_abs
        colorbar_title = "|⟨i|Q̂|j⟩|"

    fig.add_trace(
        go.Heatmap(
            z=overlap_plot,
            colorscale="Viridis",
            colorbar=dict(title=colorbar_title),
            hovertemplate="i=%{y}<br>j=%{x}<br>overlap=%{customdata:.3e}<extra></extra>",
            customdata=overlap_abs,
        )
    )

    fig.update_layout(
        title=title,
        xaxis_title="Final State j",
        yaxis_title="Initial State i",
        template="plotly_white",
        font=dict(size=14),
        width=800,
        height=700,
    )

    return fig


__all__ = [
    "plot_potential",
    "plot_capture_coefficient",
    "plot_eigenvalue_spectrum",
    "plot_configuration_coordinate",
    "plot_overlap_matrix",
]
