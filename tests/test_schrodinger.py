"""
Tests for Schrödinger equation solver.

Validates the 1D TISE solver against analytical solutions and Julia results.
"""

import pytest
import numpy as np
from carriercapture.core.schrodinger import (
    solve_schrodinger_1d,
    build_hamiltonian_1d,
    normalize_wavefunctions,
)


class TestHarmonicOscillator:
    """Test suite for harmonic oscillator (analytical solution)."""

    def test_eigenvalues_match_analytical(self):
        """
        Test that harmonic oscillator eigenvalues match E_n = ℏω(n + 1/2).

        This replicates the Julia test from test_tise.jl which uses:
        - NQ = 10000 points
        - Q range: -4 to 4
        - ℏω = 1 eV
        - nev = 10
        - Expected: E_n = 0.5 + n with tolerance 1e-3
        """
        # Setup (matching Julia test)
        NQ = 10000
        Qi, Qf = -4, 4
        Q = np.linspace(Qi, Qf, NQ)
        hw = 1.0  # eV
        nev = 10

        # Harmonic potential: V(Q) = a * Q²
        # where a = (amu/2) * (ℏω / (ℏc * 1e10))²
        # The factor of 1/2 is included in a, matching Julia's harmonic function
        from carriercapture._constants import AMU, HBAR_C

        a = (AMU / 2) * (hw / (HBAR_C * 1e10)) ** 2

        def harmonic_potential(Q):
            return a * Q**2  # NOT 0.5 * a * Q²!

        # Solve
        eigenvalues, eigenvectors = solve_schrodinger_1d(
            harmonic_potential, Q, nev=nev
        )

        # Analytical solution: E_n = ℏω(n + 1/2) = hw * (n + 0.5)
        expected = hw * (np.arange(nev) + 0.5)

        # Check eigenvalues match within tolerance
        np.testing.assert_allclose(eigenvalues, expected, atol=1e-3)

    def test_ground_state_energy(self):
        """Test ground state energy E_0 = ℏω/2."""
        Q = np.linspace(-10, 10, 5000)
        hw = 0.02  # 20 meV

        from carriercapture._constants import AMU, HBAR_C

        a = (AMU / 2) * (hw / (HBAR_C * 1e10)) ** 2

        def harmonic_potential(Q):
            return a * Q**2

        eigenvalues, _ = solve_schrodinger_1d(harmonic_potential, Q, nev=5)

        expected_ground = hw / 2
        np.testing.assert_allclose(eigenvalues[0], expected_ground, rtol=1e-4)

    def test_energy_spacing(self):
        """Test that energy levels are equally spaced by ℏω."""
        Q = np.linspace(-10, 10, 5000)
        hw = 0.03

        from carriercapture._constants import AMU, HBAR_C

        a = (AMU / 2) * (hw / (HBAR_C * 1e10)) ** 2

        def harmonic_potential(Q):
            return a * Q**2

        eigenvalues, _ = solve_schrodinger_1d(harmonic_potential, Q, nev=10)

        # Spacing should be ℏω
        spacings = np.diff(eigenvalues)
        np.testing.assert_allclose(spacings, hw, rtol=1e-3)


class TestWavefunctionProperties:
    """Test wavefunction normalization and orthogonality."""

    def test_normalization(self):
        """Test that wavefunctions are normalized: ∫|ψ|² dQ = 1."""
        Q = np.linspace(-10, 10, 5000)
        hw = 0.02

        from carriercapture._constants import AMU, HBAR_C

        a = (AMU / 2) * (hw / (HBAR_C * 1e10)) ** 2

        def harmonic_potential(Q):
            return a * Q**2

        _, eigenvectors = solve_schrodinger_1d(harmonic_potential, Q, nev=10)

        dQ = Q[1] - Q[0]

        for i in range(len(eigenvectors)):
            norm_sq = np.trapezoid(eigenvectors[i, :] ** 2, dx=dQ)
            np.testing.assert_allclose(norm_sq, 1.0, rtol=1e-6)

    def test_orthogonality(self):
        """Test that wavefunctions are orthogonal: ⟨ψ_i|ψ_j⟩ = δ_ij."""
        Q = np.linspace(-10, 10, 5000)
        hw = 0.02

        from carriercapture._constants import AMU, HBAR_C

        a = (AMU / 2) * (hw / (HBAR_C * 1e10)) ** 2

        def harmonic_potential(Q):
            return a * Q**2

        _, eigenvectors = solve_schrodinger_1d(harmonic_potential, Q, nev=10)

        dQ = Q[1] - Q[0]

        # Check first few pairs
        for i in range(5):
            for j in range(i + 1, 5):
                overlap = np.trapezoid(
                    eigenvectors[i, :] * eigenvectors[j, :], dx=dQ
                )
                np.testing.assert_allclose(overlap, 0.0, atol=1e-6)


class TestNormalizationFunction:
    """Test the normalize_wavefunctions utility function."""

    def test_normalize_random_wavefunction(self):
        """Test normalization of random wavefunctions."""
        np.random.seed(42)
        wavefunctions = np.random.rand(10, 1000)
        h = 0.01

        normalized = normalize_wavefunctions(wavefunctions, h)

        # Check normalization
        norms_sq = h * np.sum(normalized**2, axis=1)
        np.testing.assert_allclose(norms_sq, 1.0, rtol=1e-10)


class TestHamiltonianConstruction:
    """Test Hamiltonian matrix construction."""

    def test_hamiltonian_is_sparse(self):
        """Test that Hamiltonian is sparse tridiagonal."""
        Q = np.linspace(-10, 10, 1000)

        def potential(Q):
            return 0.1 * Q**2

        H = build_hamiltonian_1d(potential, Q)

        # Check shape
        assert H.shape == (1000, 1000)

        # Check sparsity (tridiagonal should have 3*N - 2 nonzeros)
        expected_nnz = 3 * 1000 - 2
        assert H.nnz == expected_nnz

    def test_hamiltonian_is_hermitian(self):
        """Test that Hamiltonian is Hermitian (symmetric)."""
        Q = np.linspace(-10, 10, 100)

        def potential(Q):
            return 0.1 * Q**2

        H = build_hamiltonian_1d(potential, Q)

        # For real symmetric matrices, H = H^T
        H_dense = H.toarray()
        np.testing.assert_allclose(H_dense, H_dense.T, rtol=1e-10)


class TestErrorHandling:
    """Test error handling and edge cases."""

    def test_non_uniform_grid_raises(self):
        """Test that non-uniform grid raises ValueError."""
        Q = np.array([0, 1, 2, 4, 8])  # Non-uniform spacing

        def potential(Q):
            return Q**2

        with pytest.raises(ValueError, match="uniformly spaced"):
            build_hamiltonian_1d(potential, Q)

    def test_convergence_failure(self):
        """Test that ARPACK convergence failure is handled."""
        Q = np.linspace(-10, 10, 5000)

        def potential(Q):
            return 0.1 * Q**2

        # Request too many eigenvalues with low maxiter
        with pytest.raises(RuntimeError, match="did not converge"):
            solve_schrodinger_1d(potential, Q, nev=100, maxiter=10)
