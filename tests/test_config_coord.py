"""
Tests for ConfigCoordinate class.

Validates overlap calculations and capture coefficient calculations.
"""

import pytest
import numpy as np
from carriercapture.core import Potential, ConfigCoordinate


class TestConfigCoordinateCreation:
    """Test ConfigCoordinate initialization."""

    def test_initialization(self):
        """Test basic initialization."""
        pot_i = Potential.from_harmonic(hw=0.03)
        pot_f = Potential.from_harmonic(hw=0.02)

        cc = ConfigCoordinate(pot_i, pot_f, name="test", W=0.2, degeneracy=2)

        assert cc.name == "test"
        assert cc.W == 0.2
        assert cc.degeneracy == 2
        assert cc.pot_i is pot_i
        assert cc.pot_f is pot_f
        assert cc.overlap_matrix is None
        assert cc.capture_coefficient is None


class TestOverlapCalculation:
    """Test overlap matrix calculation."""

    def test_calculate_overlap_basic(self):
        """Test basic overlap calculation."""
        # Create two harmonic potentials with different frequencies
        pot_i = Potential.from_harmonic(hw=0.03, Q0=0.0, E0=1.0, npoints=3000)
        pot_f = Potential.from_harmonic(hw=0.02, Q0=5.0, E0=0.0, npoints=3000)

        pot_i.solve(nev=30)
        pot_f.solve(nev=20)

        cc = ConfigCoordinate(pot_i, pot_f, W=0.2)
        cc.calculate_overlap(Q0=5.0, cutoff=0.25, sigma=0.025)

        assert cc.overlap_matrix is not None
        assert cc.delta_matrix is not None
        assert cc.overlap_matrix.shape == (30, 20)
        assert cc.delta_matrix.shape == (30, 20)

    def test_overlap_matrix_properties(self):
        """Test properties of overlap matrix."""
        pot_i = Potential.from_harmonic(hw=0.03, Q0=0.0, E0=0.0, npoints=3000)
        pot_f = Potential.from_harmonic(hw=0.03, Q0=0.0, E0=0.0, npoints=3000)

        pot_i.solve(nev=20)
        pot_f.solve(nev=20)

        cc = ConfigCoordinate(pot_i, pot_f)
        cc.calculate_overlap(Q0=0.0, cutoff=0.5, sigma=0.025)

        # For identical potentials with Q0=0, diagonal should be zero
        # (because ⟨ψ_i|(Q-0)|ψ_i⟩ = 0 by symmetry)
        diagonal = np.diag(cc.overlap_matrix)
        np.testing.assert_allclose(diagonal, 0.0, atol=1e-10)

    def test_energy_cutoff_filtering(self):
        """Test that energy cutoff filters overlaps correctly."""
        pot_i = Potential.from_harmonic(hw=0.05, Q0=0.0, npoints=2000)
        pot_f = Potential.from_harmonic(hw=0.02, Q0=5.0, npoints=2000)

        pot_i.solve(nev=20)
        pot_f.solve(nev=20)

        cc = ConfigCoordinate(pot_i, pot_f)

        # Very small cutoff should give mostly zeros
        cc.calculate_overlap(Q0=5.0, cutoff=0.01, sigma=0.01)
        n_nonzero_small = np.count_nonzero(cc.overlap_matrix)

        # Larger cutoff should give more non-zeros
        cc.calculate_overlap(Q0=5.0, cutoff=0.5, sigma=0.01)
        n_nonzero_large = np.count_nonzero(cc.overlap_matrix)

        assert n_nonzero_large > n_nonzero_small

    def test_overlap_without_solved_potentials_raises(self):
        """Test that calculating overlap without solving raises error."""
        pot_i = Potential.from_harmonic(hw=0.03)
        pot_f = Potential.from_harmonic(hw=0.02)

        cc = ConfigCoordinate(pot_i, pot_f)

        with pytest.raises(ValueError, match="must be solved"):
            cc.calculate_overlap(Q0=0.0)

    def test_overlap_incompatible_grids_raises(self):
        """Test that incompatible Q grids raise error."""
        pot_i = Potential.from_harmonic(hw=0.03, Q_range=(-10, 10), npoints=1000)
        pot_f = Potential.from_harmonic(hw=0.02, Q_range=(-20, 20), npoints=2000)

        pot_i.solve(nev=10)
        pot_f.solve(nev=10)

        cc = ConfigCoordinate(pot_i, pot_f)

        # Will raise ValueError about incompatible grids (or broadcasting error)
        with pytest.raises(ValueError):
            cc.calculate_overlap(Q0=0.0)


class TestCaptureCoefficient:
    """Test capture coefficient calculation."""

    def test_calculate_capture_coefficient_basic(self):
        """Test basic capture coefficient calculation."""
        pot_i = Potential.from_harmonic(hw=0.03, Q0=0.0, E0=1.0, npoints=3000)
        pot_f = Potential.from_harmonic(hw=0.02, Q0=10.0, E0=0.0, npoints=3000)

        pot_i.solve(nev=60)
        pot_f.solve(nev=40)

        cc = ConfigCoordinate(pot_i, pot_f, W=0.2, degeneracy=1)
        cc.calculate_overlap(Q0=10.0, cutoff=0.25, sigma=0.025)

        temperature = np.linspace(100, 500, 50)
        volume = 1e-21  # cm³

        cc.calculate_capture_coefficient(volume=volume, temperature=temperature)

        assert cc.capture_coefficient is not None
        assert cc.temperature is not None
        assert cc.partial_capture_coefficient is not None
        assert len(cc.capture_coefficient) == 50
        assert cc.capture_coefficient.shape == (50,)
        assert cc.partial_capture_coefficient.shape == (60, 40, 50)

    def test_capture_coefficient_temperature_dependence(self):
        """Test that capture coefficient varies with temperature."""
        pot_i = Potential.from_harmonic(hw=0.03, Q0=0.0, E0=1.0, npoints=3000)
        pot_f = Potential.from_harmonic(hw=0.02, Q0=10.0, E0=0.0, npoints=3000)

        pot_i.solve(nev=60)
        pot_f.solve(nev=40)

        cc = ConfigCoordinate(pot_i, pot_f, W=0.2)
        cc.calculate_overlap(Q0=10.0)

        temperature = np.array([100.0, 200.0, 300.0, 400.0, 500.0])
        cc.calculate_capture_coefficient(volume=1e-21, temperature=temperature)

        # Capture coefficient should vary with temperature
        # Check that not all values are the same
        assert len(np.unique(cc.capture_coefficient)) > 1
        # Check that there's a reasonable spread
        assert np.std(cc.capture_coefficient) > 0

    def test_capture_coefficient_physical_units(self):
        """Test that capture coefficient has correct units (cm³/s)."""
        pot_i = Potential.from_harmonic(hw=0.03, Q0=0.0, E0=1.0, npoints=2000)
        pot_f = Potential.from_harmonic(hw=0.02, Q0=10.0, E0=0.0, npoints=2000)

        pot_i.solve(nev=30)
        pot_f.solve(nev=20)

        cc = ConfigCoordinate(pot_i, pot_f, W=0.2)
        cc.calculate_overlap(Q0=10.0)

        temperature = np.linspace(100, 500, 20)
        cc.calculate_capture_coefficient(volume=1e-21, temperature=temperature)

        # Should be positive and reasonable magnitude (1e-15 to 1e-5 cm³/s typical)
        assert np.all(cc.capture_coefficient > 0)
        assert np.all(cc.capture_coefficient < 1.0)  # Upper bound check

    def test_partition_function_convergence_check(self):
        """Test that insufficient eigenvalues raises convergence error."""
        pot_i = Potential.from_harmonic(hw=0.03, Q0=0.0, E0=1.0, npoints=2000)
        pot_f = Potential.from_harmonic(hw=0.02, Q0=10.0, E0=0.0, npoints=2000)

        # Only solve for very few eigenvalues
        pot_i.solve(nev=3)
        pot_f.solve(nev=20)

        cc = ConfigCoordinate(pot_i, pot_f, W=0.2)
        cc.calculate_overlap(Q0=10.0)

        # At high temperature with few eigenvalues, partition function won't converge
        temperature = np.array([800.0])

        with pytest.raises(ValueError, match="Partition function not converged"):
            cc.calculate_capture_coefficient(volume=1e-21, temperature=temperature)

    def test_capture_without_overlap_raises(self):
        """Test that calculating capture without overlap raises error."""
        pot_i = Potential.from_harmonic(hw=0.03)
        pot_f = Potential.from_harmonic(hw=0.02)

        pot_i.solve(nev=10)
        pot_f.solve(nev=10)

        cc = ConfigCoordinate(pot_i, pot_f)

        with pytest.raises(ValueError, match="Must calculate overlaps"):
            cc.calculate_capture_coefficient(
                volume=1e-21, temperature=np.array([300.0])
            )


class TestFullWorkflow:
    """Test complete workflow from potentials to capture coefficient."""

    def test_simple_harmonic_case(self):
        """Test full workflow with simple harmonic potentials."""
        # Create two displaced harmonic oscillators
        hw_i = 0.03  # eV
        hw_f = 0.02  # eV
        dQ = 10.0  # amu^0.5·Å
        dE = 1.0  # eV

        pot_i = Potential.from_harmonic(hw=hw_i, Q0=0.0, E0=dE, npoints=4000)
        pot_f = Potential.from_harmonic(hw=hw_f, Q0=dQ, E0=0.0, npoints=4000)

        # Solve Schrödinger equation
        pot_i.solve(nev=100)
        pot_f.solve(nev=60)

        # Create configuration coordinate
        cc = ConfigCoordinate(pot_i, pot_f, W=0.205, degeneracy=1)

        # Calculate overlap
        cc.calculate_overlap(Q0=dQ, cutoff=0.25, sigma=0.0075)

        # Calculate capture coefficient
        temperature = np.linspace(100, 500, 50)
        volume = 1.28463e-21  # cm³
        cc.calculate_capture_coefficient(volume=volume, temperature=temperature)

        # Basic sanity checks
        assert cc.capture_coefficient is not None
        assert np.all(cc.capture_coefficient > 0)
        assert np.all(np.isfinite(cc.capture_coefficient))

        # Check that capture coefficient varies with temperature
        # Values should not all be the same
        assert len(np.unique(cc.capture_coefficient)) > 1
        # Check reasonable variation (std should be > 1% of mean)
        if np.mean(cc.capture_coefficient) > 0:
            rel_std = np.std(cc.capture_coefficient) / np.mean(cc.capture_coefficient)
            assert rel_std > 0.01


class TestSerialization:
    """Test serialization and deserialization."""

    def test_to_dict_from_dict(self):
        """Test round-trip serialization."""
        pot_i = Potential.from_harmonic(hw=0.03, Q0=0.0, E0=1.0, npoints=1000)
        pot_f = Potential.from_harmonic(hw=0.02, Q0=5.0, E0=0.0, npoints=1000)

        pot_i.solve(nev=20)
        pot_f.solve(nev=15)

        cc = ConfigCoordinate(pot_i, pot_f, name="test", W=0.2, degeneracy=2)
        cc.calculate_overlap(Q0=5.0)
        cc.calculate_capture_coefficient(
            volume=1e-21, temperature=np.linspace(100, 300, 10)
        )

        # Serialize
        data = cc.to_dict()

        assert data["name"] == "test"
        assert data["W"] == 0.2
        assert data["degeneracy"] == 2
        assert "overlap_matrix" in data
        assert "capture_coefficient" in data

        # Deserialize
        cc2 = ConfigCoordinate.from_dict(data)

        assert cc2.name == cc.name
        assert cc2.W == cc.W
        assert cc2.degeneracy == cc.degeneracy
        np.testing.assert_allclose(cc2.overlap_matrix, cc.overlap_matrix)
        np.testing.assert_allclose(cc2.capture_coefficient, cc.capture_coefficient)
        np.testing.assert_allclose(
            cc2.pot_i.eigenvalues, cc.pot_i.eigenvalues, rtol=1e-6
        )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
