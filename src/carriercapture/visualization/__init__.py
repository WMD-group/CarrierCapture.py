"""Visualization utilities for CarrierCapture."""

from .static import (
    plot_potential,
    plot_capture_coefficient,
    plot_eigenvalue_spectrum,
    plot_configuration_coordinate,
    plot_overlap_matrix,
)

from .themes import (
    COLORS,
    POTENTIAL_COLORS,
    get_default_layout,
    apply_publication_style,
    get_colorscale,
    format_scientific,
    create_dash_theme,
)

from .interactive import (
    create_app,
    run_server,
)

__all__ = [
    # Static plots
    "plot_potential",
    "plot_capture_coefficient",
    "plot_eigenvalue_spectrum",
    "plot_configuration_coordinate",
    "plot_overlap_matrix",
    # Themes
    "COLORS",
    "POTENTIAL_COLORS",
    "get_default_layout",
    "apply_publication_style",
    "get_colorscale",
    "format_scientific",
    "create_dash_theme",
    # Interactive dashboard
    "create_app",
    "run_server",
]
