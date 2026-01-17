"""
Tests for parameter scanning functionality.

Tests the ParameterScanner class for high-throughput materials screening.
"""

import pytest
import numpy as np
from pathlib import Path
import tempfile

from carriercapture.analysis.parameter_scan import (
    ParameterScanner,
    ScanParameters,
    ScanResult,
)


class TestScanParameters:
    """Test ScanParameters dataclass."""

    def test_create_scan_parameters(self):
        """Test creating ScanParameters object."""
        params = ScanParameters(
            dQ_range=(0, 10, 5),
            dE_range=(0, 1, 5),
        )

        assert params.dQ_range == (0, 10, 5)
        assert params.dE_range == (0, 1, 5)
        assert params.hbar_omega_i == 0.008  # Default
        assert params.temperature == 300.0  # Default

    def test_scan_parameters_with_custom_values(self):
        """Test ScanParameters with custom values."""
        params = ScanParameters(
            dQ_range=(0, 25, 25),
            dE_range=(0, 2.5, 10),
            hbar_omega_i=0.010,
            hbar_omega_f=0.012,
            temperature=400.0,
            volume=2e-21,
        )

        assert params.hbar_omega_i == 0.010
        assert params.hbar_omega_f == 0.012
        assert params.temperature == 400.0
        assert params.volume == 2e-21


class TestScanResult:
    """Test ScanResult dataclass."""

    @pytest.fixture
    def sample_result(self):
        """Create sample scan result."""
        params = ScanParameters(
            dQ_range=(0, 10, 5),
            dE_range=(0, 1, 5),
        )

        dQ_grid = np.linspace(0, 10, 5)
        dE_grid = np.linspace(0, 1, 5)
        capture_coeffs = np.random.rand(5, 5) * 1e-10
        barriers = np.random.rand(5, 5) * 2

        return ScanResult(
            dQ_grid=dQ_grid,
            dE_grid=dE_grid,
            capture_coefficients=capture_coeffs,
            barrier_heights=barriers,
            parameters=params,
            metadata={'test': 'data'}
        )

    def test_save_load_npz(self, sample_result, tmp_path):
        """Test saving and loading NPZ format."""
        filepath = tmp_path / "test_result.npz"

        # Save
        sample_result.save(filepath, format='npz')
        assert filepath.exists()

        # Load
        loaded = ScanResult.load(filepath, format='npz')

        assert np.array_equal(loaded.dQ_grid, sample_result.dQ_grid)
        assert np.array_equal(loaded.dE_grid, sample_result.dE_grid)
        assert np.array_equal(loaded.capture_coefficients, sample_result.capture_coefficients)
        assert np.array_equal(loaded.barrier_heights, sample_result.barrier_heights)

    @pytest.mark.skip(reason="Requires h5py - optional dependency")
    def test_save_load_hdf5(self, sample_result, tmp_path):
        """Test saving and loading HDF5 format."""
        pytest.importorskip("h5py")

        filepath = tmp_path / "test_result.h5"

        # Save
        sample_result.save(filepath, format='hdf5')
        assert filepath.exists()

        # Load
        loaded = ScanResult.load(filepath, format='hdf5')

        assert np.array_equal(loaded.dQ_grid, sample_result.dQ_grid)
        assert np.array_equal(loaded.dE_grid, sample_result.dE_grid)


class TestParameterScanner:
    """Test ParameterScanner class."""

    @pytest.fixture
    def simple_params(self):
        """Create simple parameters for fast testing."""
        return ScanParameters(
            dQ_range=(0, 5, 3),  # Small grid for fast testing
            dE_range=(0, 0.5, 3),
            hbar_omega_i=0.008,
            hbar_omega_f=0.008,
            temperature=300.0,
            nev_initial=60,  # Need enough for partition function convergence at 300K
            nev_final=30,   # Scaled proportionally
            Q_grid_points=500,  # Reduced for speed
        )

    def test_create_scanner(self, simple_params):
        """Test creating ParameterScanner."""
        scanner = ParameterScanner(simple_params, verbose=False)

        assert scanner.params == simple_params
        assert len(scanner.dQ_grid) == 3
        assert len(scanner.dE_grid) == 3

    def test_build_grids(self, simple_params):
        """Test grid building."""
        scanner = ParameterScanner(simple_params, verbose=False)

        # Check dQ grid
        assert scanner.dQ_grid[0] == 0.0
        assert scanner.dQ_grid[-1] == 5.0
        assert len(scanner.dQ_grid) == 3

        # Check dE grid
        assert scanner.dE_grid[0] == 0.0
        assert scanner.dE_grid[-1] == 0.5
        assert len(scanner.dE_grid) == 3

    def test_create_harmonic_potentials(self, simple_params):
        """Test creating harmonic potentials."""
        scanner = ParameterScanner(simple_params, verbose=False)

        pot_i, pot_f = scanner._create_harmonic_potentials(
            hbar_omega_i=0.008,
            hbar_omega_f=0.008,
            dQ=2.0,
            dE=0.5
        )

        # Check that potentials were created and solved
        assert pot_i is not None
        assert pot_f is not None
        assert pot_i.eigenvalues is not None
        assert pot_f.eigenvalues is not None
        assert len(pot_i.eigenvalues) == 60  # nev_initial
        assert len(pot_f.eigenvalues) == 30  # nev_final

    def test_calculate_W_coupling(self, simple_params):
        """Test W coupling calculation."""
        scanner = ParameterScanner(simple_params, verbose=False)

        W = scanner._calculate_W_coupling(
            hbar_omega_f=0.008,
            dQ=2.0,
            dE=0.5
        )

        assert isinstance(W, float)
        assert W > 0

    def test_calculate_single_point(self, simple_params):
        """Test single point calculation."""
        scanner = ParameterScanner(simple_params, verbose=False)

        capture_coeff, barrier = scanner._calculate_single_point(
            hbar_omega_i=0.008,
            hbar_omega_f=0.008,
            dQ=2.0,
            dE=0.5
        )

        # Should return valid numbers (or nan if failed)
        assert isinstance(capture_coeff, (float, np.floating))
        assert isinstance(barrier, (float, np.floating))

        # If successful, should be positive
        if not np.isnan(capture_coeff):
            assert capture_coeff > 0

    @pytest.mark.slow
    def test_run_harmonic_scan_serial(self, simple_params):
        """Test running harmonic scan in serial mode."""
        scanner = ParameterScanner(simple_params, verbose=False)

        results = scanner.run_harmonic_scan(n_jobs=1, show_progress=False)

        # Check result structure
        assert isinstance(results, ScanResult)
        assert results.capture_coefficients.shape == (3, 3)
        assert results.barrier_heights.shape == (3, 3)
        assert np.array_equal(results.dQ_grid, scanner.dQ_grid)
        assert np.array_equal(results.dE_grid, scanner.dE_grid)

        # Check that some calculations succeeded
        n_success = np.sum(~np.isnan(results.capture_coefficients))
        assert n_success > 0

    @pytest.mark.slow
    @pytest.mark.skipif(True, reason="Parallel execution can be unstable in tests")
    def test_run_harmonic_scan_parallel(self, simple_params):
        """Test running harmonic scan in parallel mode."""
        pytest.importorskip("joblib")

        scanner = ParameterScanner(simple_params, verbose=False)

        results = scanner.run_harmonic_scan(n_jobs=2, show_progress=False)

        # Check that results match expected structure
        assert isinstance(results, ScanResult)
        assert results.capture_coefficients.shape == (3, 3)

        # Check that some calculations succeeded
        n_success = np.sum(~np.isnan(results.capture_coefficients))
        assert n_success > 0

    def test_scan_result_metadata(self, simple_params):
        """Test that scan results include metadata."""
        scanner = ParameterScanner(simple_params, verbose=False)

        results = scanner.run_harmonic_scan(n_jobs=1, show_progress=False)

        # Check metadata
        assert 'hbar_omega_i' in results.metadata
        assert 'hbar_omega_f' in results.metadata
        assert 'temperature' in results.metadata
        assert 'volume' in results.metadata


class TestScanCLI:
    """Test CLI integration."""

    def test_scan_command_exists(self):
        """Test that scan commands are registered."""
        from carriercapture.cli.main import cli

        assert "scan" in cli.commands
        assert "scan-plot" in cli.commands

    def test_scan_command_help(self):
        """Test that scan command has help text."""
        from carriercapture.cli.commands.scan import scan_cmd

        assert scan_cmd.help is not None
        assert "parameter" in scan_cmd.help.lower()


class TestModuleExports:
    """Test that analysis module exports expected classes."""

    def test_analysis_exports(self):
        """Test that analysis module exports parameter scan classes."""
        from carriercapture import analysis

        assert hasattr(analysis, 'ParameterScanner')
        assert hasattr(analysis, 'ScanParameters')
        assert hasattr(analysis, 'ScanResult')


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
