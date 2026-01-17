"""
Tests for visualization functions.

Tests the static plotting functions and basic dashboard creation.
"""

import pytest
import numpy as np
from pathlib import Path

# Try to import plotly - skip tests if not available
pytest.importorskip("plotly")

from carriercapture.core.potential import Potential
from carriercapture.core.config_coord import ConfigCoordinate
from carriercapture.visualization.static import (
    plot_potential,
    plot_capture_coefficient,
    plot_eigenvalue_spectrum,
    plot_configuration_coordinate,
    plot_overlap_matrix,
)
from carriercapture.visualization.themes import (
    COLORS,
    get_default_layout,
    apply_publication_style,
    format_scientific,
)


class TestStaticPlots:
    """Test static plotting functions."""

    @pytest.fixture
    def harmonic_potential(self):
        """Create a simple harmonic potential for testing."""
        Q = np.linspace(-5, 5, 100)
        E = 0.5 * Q**2
        pot = Potential(Q_data=Q, E_data=E, name="Harmonic")
        pot.fit(fit_type="spline")
        pot.solve(nev=10)
        return pot

    def test_plot_potential_basic(self, harmonic_potential):
        """Test basic potential plotting."""
        fig = plot_potential(harmonic_potential, show_wavefunctions=False)

        assert fig is not None
        assert len(fig.data) >= 1  # At least the potential curve
        assert "Harmonic" in fig.layout.title.text or "Harmonic" in str(fig.data[0].name)

    def test_plot_potential_with_wavefunctions(self, harmonic_potential):
        """Test potential plotting with wavefunctions."""
        fig = plot_potential(harmonic_potential, show_wavefunctions=True, wf_sampling=2)

        assert fig is not None
        # Should have more traces (wavefunctions add multiple traces)
        assert len(fig.data) > 2

    def test_plot_potential_with_eigenvalues(self, harmonic_potential):
        """Test potential plotting with eigenvalue lines."""
        fig = plot_potential(harmonic_potential, show_eigenvalues=True)

        assert fig is not None
        # Check that eigenvalues are shown (as horizontal lines)
        assert len(fig.data) > 1

    def test_plot_eigenvalue_spectrum(self, harmonic_potential):
        """Test eigenvalue spectrum plotting."""
        fig = plot_eigenvalue_spectrum(harmonic_potential, max_levels=10)

        assert fig is not None
        assert len(fig.data) == 10  # One line per eigenvalue
        assert "Eigenvalue Spectrum" in fig.layout.title.text

    def test_plot_eigenvalue_spectrum_requires_solved(self):
        """Test that plotting spectrum requires solved potential."""
        pot = Potential(Q_data=np.linspace(-5, 5, 50), E_data=np.zeros(50))

        with pytest.raises(ValueError, match="must be solved"):
            plot_eigenvalue_spectrum(pot)

    def test_plot_capture_coefficient(self, harmonic_potential):
        """Test capture coefficient plotting."""
        # Create a second potential
        Q = np.linspace(-5, 5, 100)
        E = 0.5 * (Q - 2)**2
        pot_f = Potential(Q_data=Q, E_data=E, name="Final")
        pot_f.fit(fit_type="spline")
        pot_f.solve(nev=10)

        # Create ConfigCoordinate
        cc = ConfigCoordinate(
            pot_i=harmonic_potential,
            pot_f=pot_f,
            W=0.1,
            degeneracy=1
        )

        # Calculate (minimal for testing)
        cc.calculate_overlap(Q0=1.0, cutoff=0.5, sigma=0.05)

        temperature = np.linspace(100, 500, 10)
        cc.calculate_capture_coefficient(volume=1e-21, temperature=temperature)

        # Plot
        fig = plot_capture_coefficient(cc)

        assert fig is not None
        assert len(fig.data) >= 1
        assert "1000/T" in fig.layout.xaxis.title.text
        assert "log" in fig.layout.yaxis.title.text.lower()

    def test_plot_configuration_coordinate(self, harmonic_potential):
        """Test configuration coordinate diagram."""
        # Create a second potential
        Q = np.linspace(-5, 5, 100)
        E = 0.5 * (Q - 3)**2 + 0.5
        pot_f = Potential(Q_data=Q, E_data=E, name="Final")
        pot_f.fit(fit_type="spline")

        fig = plot_configuration_coordinate(
            harmonic_potential,
            pot_f,
            Q0=1.5,
            show_wavefunctions=False
        )

        assert fig is not None
        # Should have both potentials plus Q0 line
        assert len(fig.data) >= 3

    def test_plot_overlap_matrix(self, harmonic_potential):
        """Test overlap matrix plotting."""
        # Create a second potential
        Q = np.linspace(-5, 5, 100)
        E = 0.5 * (Q - 2)**2
        pot_f = Potential(Q_data=Q, E_data=E, name="Final")
        pot_f.fit(fit_type="spline")
        pot_f.solve(nev=10)

        # Create ConfigCoordinate and calculate overlap
        cc = ConfigCoordinate(
            pot_i=harmonic_potential,
            pot_f=pot_f,
            W=0.1,
            degeneracy=1
        )
        cc.calculate_overlap(Q0=1.0, cutoff=0.5, sigma=0.05)

        # Plot
        fig = plot_overlap_matrix(cc)

        assert fig is not None
        assert len(fig.data) == 1  # Heatmap
        assert fig.data[0].type == "heatmap"


class TestThemes:
    """Test theming utilities."""

    def test_colors_defined(self):
        """Test that color palette is defined."""
        assert "primary" in COLORS
        assert "secondary" in COLORS
        assert isinstance(COLORS["primary"], str)

    def test_get_default_layout(self):
        """Test default layout generation."""
        layout = get_default_layout(
            title="Test Plot",
            xaxis_title="X Axis",
            yaxis_title="Y Axis"
        )

        assert isinstance(layout, dict)
        assert "title" in layout
        assert "xaxis" in layout
        assert "yaxis" in layout

    def test_apply_publication_style(self):
        """Test applying publication style to figure."""
        import plotly.graph_objects as go

        fig = go.Figure()
        fig = apply_publication_style(fig)

        assert fig.layout.plot_bgcolor == "white"
        assert fig.layout.paper_bgcolor == "white"

    def test_format_scientific(self):
        """Test scientific notation formatting."""
        result = format_scientific(1.23e-10, precision=2)

        assert isinstance(result, str)
        assert "×10" in result or "10" in result


class TestVisualizationCLI:
    """Test CLI integration."""

    def test_viz_command_exists(self):
        """Test that viz commands are registered."""
        from carriercapture.cli.main import cli

        # Check that viz command is registered
        assert "viz" in cli.commands
        assert "plot" in cli.commands

    def test_viz_command_help(self):
        """Test that viz command has help text."""
        from carriercapture.cli.commands.viz import viz_cmd

        assert viz_cmd.help is not None
        assert "dashboard" in viz_cmd.help.lower()


class TestInteractiveApp:
    """Test Dash app creation."""

    @pytest.mark.skip(reason="Requires Dash to be fully configured")
    def test_create_app(self):
        """Test that Dash app can be created."""
        from carriercapture.visualization.interactive import create_app

        app = create_app()
        assert app is not None
        assert hasattr(app, "layout")


class TestExports:
    """Test that visualization module exports expected functions."""

    def test_static_exports(self):
        """Test that static module exports plotting functions."""
        from carriercapture.visualization import static

        expected = [
            "plot_potential",
            "plot_capture_coefficient",
            "plot_eigenvalue_spectrum",
            "plot_configuration_coordinate",
            "plot_overlap_matrix",
        ]

        for func in expected:
            assert hasattr(static, func), f"Missing export: {func}"

    def test_themes_exports(self):
        """Test that themes module exports styling functions."""
        from carriercapture.visualization import themes

        expected = [
            "COLORS",
            "get_default_layout",
            "apply_publication_style",
        ]

        for item in expected:
            assert hasattr(themes, item), f"Missing export: {item}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
