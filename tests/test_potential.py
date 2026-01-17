"""
Tests for Potential class.

Validates potential fitting, solving, and serialization.
"""

import pytest
import numpy as np
import tempfile
from pathlib import Path
from carriercapture.core import Potential


class TestPotentialCreation:
    """Test Potential object creation."""

    def test_default_initialization(self):
        """Test default Potential initialization."""
        pot = Potential()
        assert pot.name == ""
        assert pot.Q0 == 0.0
        assert pot.E0 == 0.0
        assert pot.nev == 30
        assert pot.Q_data is None
        assert pot.fit_func is None

    def test_initialization_with_params(self):
        """Test Potential initialization with parameters."""
        pot = Potential(name="test", Q0=5.0, E0=1.5, nev=60)
        assert pot.name == "test"
        assert pot.Q0 == 5.0
        assert pot.E0 == 1.5
        assert pot.nev == 60


class TestHarmonicPotential:
    """Test harmonic potential creation and solving."""

    def test_from_harmonic(self):
        """Test creating harmonic potential."""
        pot = Potential.from_harmonic(hw=0.02, Q0=0.0, E0=0.0)

        assert pot.fit_type == "harmonic"
        assert pot.Q is not None
        assert pot.E is not None
        assert pot.fit_func is not None

    def test_harmonic_eigenvalues(self):
        """Test that harmonic potential gives correct eigenvalues."""
        hw = 0.02
        pot = Potential.from_harmonic(hw=hw, npoints=5000)
        pot.solve(nev=10)

        # E_n = ℏω(n + 1/2)
        expected = hw * (np.arange(10) + 0.5)
        np.testing.assert_allclose(pot.eigenvalues, expected, rtol=1e-3)

    def test_harmonic_ground_state(self):
        """Test ground state energy."""
        hw = 0.03
        pot = Potential.from_harmonic(hw=hw)
        pot.solve(nev=5)

        expected_ground = hw / 2
        np.testing.assert_allclose(pot.eigenvalues[0], expected_ground, rtol=1e-4)

    def test_harmonic_evaluation(self):
        """Test evaluating harmonic potential."""
        pot = Potential.from_harmonic(hw=0.02, Q0=0.0, E0=0.0)

        # At Q=0, should be E0
        np.testing.assert_allclose(pot(0.0), 0.0, rtol=1e-6)

        # Should be symmetric
        np.testing.assert_allclose(pot(5.0), pot(-5.0), rtol=1e-6)


class TestFileIO:
    """Test loading potential from file."""

    def test_from_file_comma_separated(self):
        """Test loading comma-separated data file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write("# Q, E\n")
            f.write("0.0, 1.0\n")
            f.write("1.0, 1.5\n")
            f.write("2.0, 2.5\n")
            f.write("3.0, 4.0\n")
            filepath = f.name

        try:
            pot = Potential.from_file(filepath, name="test")
            assert pot.name == "test"
            assert len(pot.Q_data) == 4
            assert len(pot.E_data) == 4
            assert pot.Q0 == pytest.approx(0.0)
            assert pot.E0 == pytest.approx(1.0)
        finally:
            Path(filepath).unlink()

    def test_from_file_space_separated(self):
        """Test loading space-separated data file."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".dat", delete=False
        ) as f:
            f.write("# Q E\n")
            f.write("0.0 1.0\n")
            f.write("1.0 1.5\n")
            f.write("2.0 2.5\n")
            filepath = f.name

        try:
            pot = Potential.from_file(filepath)
            assert len(pot.Q_data) == 3
            assert pot.E0 == pytest.approx(1.0)
        finally:
            Path(filepath).unlink()


class TestFitting:
    """Test potential fitting."""

    def test_fit_spline(self):
        """Test spline fitting."""
        # Create some data
        Q_data = np.linspace(-10, 10, 25)
        E_data = 0.01 * Q_data**2 + 0.5

        pot = Potential()
        pot.Q_data = Q_data
        pot.E_data = E_data
        pot.Q = np.linspace(-10, 10, 1000)

        pot.fit(fit_type="spline", order=2, smoothness=0.0)

        assert pot.fit_func is not None
        assert pot.fit_type == "spline"
        assert len(pot.E) == 1000

        # Check that fit passes through data points (interpolating spline)
        for Q, E in zip(Q_data, E_data):
            np.testing.assert_allclose(pot(Q), E, rtol=1e-3)

    def test_fit_harmonic_to_data(self):
        """Test fitting harmonic potential to parabolic data."""
        Q_data = np.linspace(-10, 10, 25)
        hw = 0.02
        from carriercapture._constants import AMU, HBAR_C

        a = (AMU / 2) * (hw / (HBAR_C * 1e10)) ** 2
        E_data = a * Q_data**2  # a already includes factor of 1/2

        pot = Potential()
        pot.Q_data = Q_data
        pot.E_data = E_data
        pot.Q0 = 0.0
        pot.E0 = 0.0
        pot.Q = np.linspace(-10, 10, 1000)

        pot.fit(fit_type="harmonic", hw=hw)

        assert pot.fit_func is not None
        assert pot.fit_type == "harmonic"

        # Check fit quality
        E_fit = pot(Q_data)
        np.testing.assert_allclose(E_fit, E_data, rtol=1e-6)


class TestSolving:
    """Test Schrödinger equation solving."""

    def test_solve_updates_eigenvalues(self):
        """Test that solve() updates eigenvalues and eigenvectors."""
        pot = Potential.from_harmonic(hw=0.02)
        assert pot.eigenvalues is None
        assert pot.eigenvectors is None

        pot.solve(nev=30)

        assert pot.eigenvalues is not None
        assert pot.eigenvectors is not None
        assert len(pot.eigenvalues) == 30
        assert pot.eigenvectors.shape == (30, len(pot.Q))

    def test_solve_without_fit_raises(self):
        """Test that solving without fitting raises error."""
        pot = Potential()
        with pytest.raises(ValueError, match="Must fit potential"):
            pot.solve()


class TestSerialization:
    """Test serialization and deserialization."""

    def test_to_dict_from_dict(self):
        """Test round-trip serialization."""
        pot = Potential.from_harmonic(hw=0.02, Q0=1.0, E0=0.5, nev=20)
        pot.solve(nev=10)

        # Serialize
        data = pot.to_dict()

        assert data["name"] == pot.name
        assert data["Q0"] == 1.0
        assert data["E0"] == 0.5
        assert data["nev"] == 10  # Updated by solve
        assert "eigenvalues" in data

        # Deserialize
        pot2 = Potential.from_dict(data)

        assert pot2.name == pot.name
        assert pot2.Q0 == pot.Q0
        assert pot2.E0 == pot.E0
        assert pot2.nev == pot.nev
        np.testing.assert_allclose(pot2.Q, pot.Q)
        np.testing.assert_allclose(pot2.eigenvalues, pot.eigenvalues)

    def test_copy(self):
        """Test deep copy."""
        pot = Potential.from_harmonic(hw=0.02)
        pot.solve(nev=10)

        pot2 = pot.copy()

        # Modify original
        pot.name = "modified"
        pot.eigenvalues[0] = 999.0

        # Copy should be unchanged
        assert pot2.name != "modified"
        assert pot2.eigenvalues[0] != 999.0


class TestEvaluation:
    """Test potential evaluation."""

    def test_call_single_value(self):
        """Test evaluating potential at single Q value."""
        pot = Potential.from_harmonic(hw=0.02, Q0=0.0, E0=0.0)
        E = pot(5.0)
        assert isinstance(E, (float, np.floating))

    def test_call_array(self):
        """Test evaluating potential at array of Q values."""
        pot = Potential.from_harmonic(hw=0.02)
        Q_test = np.array([0.0, 1.0, 2.0, 3.0])
        E = pot(Q_test)
        assert isinstance(E, np.ndarray)
        assert len(E) == 4

    def test_call_without_fit_raises(self):
        """Test that calling unfitted potential raises error."""
        pot = Potential()
        with pytest.raises(ValueError, match="Must fit potential"):
            pot(5.0)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
