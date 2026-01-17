"""
Tests for interactive Dash dashboard.

Tests the create_app function and basic dashboard structure.
"""

import pytest

from carriercapture.visualization.interactive import create_app, run_server


class TestDashApp:
    """Test Dash application creation."""

    def test_create_app(self):
        """Test that app can be created."""
        app = create_app(debug=False)

        assert app is not None
        assert app.title == "CarrierCapture Visualization"

    def test_app_layout(self):
        """Test that app has expected layout components."""
        app = create_app(debug=False)

        # Check that layout exists
        assert app.layout is not None

        # Check that main tabs exist
        # The layout should contain dcc.Tabs with 4 tabs
        layout_str = str(app.layout)
        assert "Potential Viewer" in layout_str
        assert "Parameter Scan" in layout_str
        assert "Comparison" in layout_str
        assert "Capture Calculation" in layout_str

    def test_app_callbacks_registered(self):
        """Test that callbacks are registered."""
        app = create_app(debug=False)

        # Check that callbacks exist
        assert len(app.callback_map) > 0


class TestDashComponents:
    """Test individual dashboard components."""

    def test_serialization_functions(self):
        """Test potential serialization/deserialization."""
        from carriercapture.visualization.interactive import (
            serialize_potential,
            deserialize_potential,
        )
        from carriercapture.core.potential import Potential
        import numpy as np

        # Create test potential
        pot = Potential.from_harmonic(hw=0.01, Q0=0.0, E0=0.0)
        pot.solve(nev=10)

        # Serialize
        data = serialize_potential(pot)

        # Check serialized data
        assert data["name"] is not None
        assert data["Q"] is not None
        assert data["E"] is not None
        assert data["eigenvalues"] is not None

        # Deserialize
        pot_restored = deserialize_potential(data)

        # Check restored potential
        assert pot_restored.name == pot.name
        assert np.allclose(pot_restored.Q, pot.Q)
        assert np.allclose(pot_restored.E, pot.E)
        assert np.allclose(pot_restored.eigenvalues, pot.eigenvalues)

    def test_scan_results_serialization(self):
        """Test scan results serialization."""
        from carriercapture.visualization.interactive import (
            serialize_scan_results,
            deserialize_scan_results,
        )
        from carriercapture.analysis.parameter_scan import (
            ScanParameters,
            ScanResult,
        )
        import numpy as np

        # Create test scan results
        params = ScanParameters(
            dQ_range=(0, 5, 3),
            dE_range=(0, 1, 3),
        )

        results = ScanResult(
            dQ_grid=np.array([0, 2.5, 5]),
            dE_grid=np.array([0, 0.5, 1]),
            capture_coefficients=np.random.rand(3, 3) * 1e-10,
            barrier_heights=np.random.rand(3, 3),
            parameters=params,
            metadata={},
        )

        # Serialize
        data = serialize_scan_results(results)

        # Check serialized data
        assert len(data["dQ_grid"]) == 3
        assert len(data["dE_grid"]) == 3

        # Deserialize
        results_restored = deserialize_scan_results(data)

        # Check restored results
        assert np.allclose(results_restored.dQ_grid, results.dQ_grid)
        assert np.allclose(results_restored.dE_grid, results.dE_grid)
        assert np.allclose(results_restored.capture_coefficients, results.capture_coefficients)

    def test_create_potential_figure(self):
        """Test potential figure creation."""
        from carriercapture.visualization.interactive import create_potential_figure
        from carriercapture.core.potential import Potential

        pot = Potential.from_harmonic(hw=0.01, Q0=0.0, E0=0.0)
        pot.solve(nev=10)

        fig = create_potential_figure(pot, ["ev", "data"], wf_scale=1.0)

        assert fig is not None
        assert len(fig.data) > 0  # Should have traces

    def test_create_scan_figure(self):
        """Test scan figure creation."""
        from carriercapture.visualization.interactive import create_scan_figure
        from carriercapture.analysis.parameter_scan import (
            ScanParameters,
            ScanResult,
        )
        import numpy as np

        params = ScanParameters(
            dQ_range=(0, 5, 3),
            dE_range=(0, 1, 3),
        )

        results = ScanResult(
            dQ_grid=np.array([0, 2.5, 5]),
            dE_grid=np.array([0, 0.5, 1]),
            capture_coefficients=np.random.rand(3, 3) * 1e-10,
            barrier_heights=np.random.rand(3, 3),
            parameters=params,
            metadata={},
        )

        # Test heatmap
        fig = create_scan_figure(results, "heatmap", log_scale=True)
        assert fig is not None
        assert len(fig.data) > 0

        # Test contour
        fig = create_scan_figure(results, "contour", log_scale=False)
        assert fig is not None
        assert len(fig.data) > 0

    def test_create_comparison_figure(self):
        """Test comparison figure creation."""
        from carriercapture.visualization.interactive import create_comparison_figure
        from carriercapture.core.potential import Potential

        pot1 = Potential.from_harmonic(hw=0.01, Q0=0.0, E0=0.0)
        pot2 = Potential.from_harmonic(hw=0.02, Q0=5.0, E0=0.5)

        fig = create_comparison_figure([pot1, pot2])

        assert fig is not None
        assert len(fig.data) == 2  # Two potentials


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
