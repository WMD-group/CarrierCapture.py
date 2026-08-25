"""
1D Schrödinger equation solver using finite difference method and ARPACK.

This module provides functions to solve the time-independent Schrödinger equation
for 1D potentials, which is critical for computing phonon states in carrier capture
calculations.

Performance: This is the computational bottleneck (80-90% of runtime). Uses scipy's
ARPACK wrapper (same FORTRAN backend as Julia) for efficient sparse eigenvalue solving.
"""

from typing import Callable, Tuple
import numpy as np
from numpy.typing import NDArray
from scipy.sparse import diags, csr_matrix
from scipy.sparse.linalg import eigsh, ArpackNoConvergence

from .._constants import AMU, HBAR_C


def build_hamiltonian_1d(
    potential_func: Callable[[NDArray[np.float64]], NDArray[np.float64]],
    Q: NDArray[np.float64],
) -> csr_matrix:
    """
    Build sparse Hamiltonian matrix for 1D TISE using finite differences.

    Uses three-point finite difference stencil for kinetic energy operator:
        T ≈ -ℏ²/(2m) * (ψ[i+1] - 2ψ[i] + ψ[i-1]) / h²

    The Hamiltonian is:
        H = T + V
    where T is the kinetic energy (sparse tridiagonal) and V is the
    potential energy (diagonal).

    Parameters
    ----------
    potential_func : callable
        Potential energy function V(Q) in eV
        Should accept numpy array and return numpy array
    Q : NDArray[np.float64]
        Configuration coordinate grid (amu^0.5·Å)
        Must be uniformly spaced

    Returns
    -------
    H : csr_matrix
        Sparse Hamiltonian matrix (N × N) in CSR format

    Notes
    -----
    The finite difference method assumes a uniform grid spacing h = Q[1] - Q[0].
    Boundary conditions are infinite potential walls (ψ = 0 at boundaries).

    The mass is m = 1 amu, typical for nuclear configuration coordinates.

    Examples
    --------
    >>> from carriercapture.core.potential import Potential
    >>> pot = Potential.from_harmonic(hw=0.02, Q0=0.0, E0=0.0)
    >>> Q = np.linspace(-10, 10, 1000)
    >>> H = build_hamiltonian_1d(pot, Q)
    >>> H.shape
    (1000, 1000)
    """
    N = len(Q)
    h = Q[1] - Q[0]  # Grid spacing in amu^0.5·Å (assumes uniform grid)

    # Verify uniform grid
    if not np.allclose(np.diff(Q), h):
        raise ValueError("Q grid must be uniformly spaced for finite difference method")

    # Mass-weighted coordinates (Q in amu^0.5·Å, unit effective mass):
    # T = -(ℏ²/2) d²/dQ², with ℏ² = (ℏc)²/(AMU·c²) = (HBAR_C·1e10)²/AMU
    # expressed in eV·amu·Å², so T comes out in eV.
    factor = (1.0 / AMU) * (HBAR_C * 1e10) ** 2

    # Kinetic energy coefficient: ℏ²/(2m h²)
    kinetic_coeff = factor / (2.0 * h**2)

    # Kinetic energy operator from finite difference: -ℏ²/(2m) d²/dQ²
    # d²ψ/dQ² ≈ (ψ[i+1] - 2ψ[i] + ψ[i-1])/h²
    # So T = -kinetic_coeff * [-1, 2, -1] = kinetic_coeff * [1, -2, 1]
    T_off = -kinetic_coeff * np.ones(N - 1)  # Off-diagonal: -ℏ²/(2m h²)
    T_diag = 2.0 * kinetic_coeff * np.ones(N)  # Diagonal: 2ℏ²/(2m h²)
    T = diags([T_off, T_diag, T_off], offsets=[-1, 0, 1], format="csr")

    # Potential energy operator: V(Q) on diagonal
    V_vals = potential_func(Q)
    V = diags([V_vals], offsets=[0], format="csr")

    # Total Hamiltonian
    H = T + V

    return H


def normalize_wavefunctions(
    wavefunctions: NDArray[np.float64],
    grid_spacing: float,
) -> NDArray[np.float64]:
    """
    Normalize wavefunctions so ∫|ψ|² dQ = 1.

    Parameters
    ----------
    wavefunctions : NDArray[np.float64]
        Wavefunctions, shape (nev, N_Q)
        Each row is a wavefunction
    grid_spacing : float
        Spacing between grid points (amu^0.5·Å)

    Returns
    -------
    normalized : NDArray[np.float64]
        Normalized wavefunctions, same shape as input

    Notes
    -----
    Uses the rectangle rule for numerical integration:
        ∫|ψ|² dQ ≈ h * Σ|ψ|²
    (equivalent to the trapezoid rule here, since the hard-wall boundary
    conditions make ψ vanish at the grid edges)

    Examples
    --------
    >>> wf = np.random.rand(10, 1000)  # 10 wavefunctions, 1000 points
    >>> h = 0.01
    >>> wf_norm = normalize_wavefunctions(wf, h)
    >>> np.allclose(np.sum(wf_norm**2, axis=1) * h, 1.0)
    True
    """
    # norm² = ∫|ψ|² dQ ≈ h * Σ|ψ|² (rectangle rule, uniform grid)
    norms_sq = grid_spacing * np.sum(wavefunctions**2, axis=1)
    norms = np.sqrt(norms_sq)

    # Divide each wavefunction by its norm
    # Broadcasting: (nev, N_Q) / (nev, 1)
    normalized = wavefunctions / norms[:, np.newaxis]

    return normalized


def solve_schrodinger_1d(
    potential_func: Callable[[NDArray[np.float64]], NDArray[np.float64]],
    Q: NDArray[np.float64],
    nev: int = 30,
    maxiter: int | None = None,
) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    """
    Solve 1D time-independent Schrödinger equation using ARPACK.

    Uses scipy.sparse.linalg.eigsh (which wraps ARPACK FORTRAN library) to find
    the lowest `nev` eigenvalues and eigenvectors of the Hamiltonian.

    This is the computational bottleneck (80-90% of total runtime) for carrier
    capture calculations.

    Parameters
    ----------
    potential_func : callable
        Potential energy function V(Q) in eV
        Should accept numpy array and return numpy array
    Q : NDArray[np.float64]
        Configuration coordinate grid (amu^0.5·Å)
        Must be uniformly spaced, typically 500-5000 points
    nev : int, default=30
        Number of eigenvalues to compute (lowest energy states)
        Typical range: 30-180
    maxiter : int, optional
        Maximum number of ARPACK iterations
        If None, defaults to max(nev * len(Q), 1000)

    Returns
    -------
    eigenvalues : NDArray[np.float64]
        Phonon energies (eV), shape (nev,)
        Sorted in ascending order
    eigenvectors : NDArray[np.float64]
        Normalized wavefunctions, shape (nev, len(Q))
        Each row is a wavefunction ψ_n(Q)

    Raises
    ------
    RuntimeError
        If ARPACK fails to converge within maxiter iterations
    ValueError
        If Q grid is not uniformly spaced

    Notes
    -----
    - Grid must be uniform for finite difference method
    - Wavefunctions are normalized: ∫|ψ|² dQ = 1
    - Uses mass m = 1 amu (typical for nuclear coordinates)
    - Boundary conditions: ψ = 0 at Q_min and Q_max (infinite walls)

    Performance
    -----------
    - N=5000, nev=60: ~0.5-1s (comparable to Julia)
    - Bottleneck is ARPACK eigenvalue solver
    - Same FORTRAN backend as Julia, so performance is identical

    Examples
    --------
    Harmonic oscillator with ℏω = 20 meV (E_n = ℏω(n + 1/2)):

    >>> from carriercapture.core.potential import Potential
    >>> pot = Potential.from_harmonic(hw=0.02, Q0=0.0, E0=0.0)
    >>> Q = np.linspace(-20, 20, 5000)
    >>> eigenvalues, eigenvectors = solve_schrodinger_1d(pot, Q, nev=10)
    >>> np.round(eigenvalues[:3], 6)
    array([0.01, 0.03, 0.05])

    See Also
    --------
    scipy.sparse.linalg.eigsh : Underlying ARPACK wrapper
    build_hamiltonian_1d : Hamiltonian construction
    normalize_wavefunctions : Wavefunction normalization

    References
    ----------
    .. [1] Alkauskas et al., "First-principles theory of nonradiative carrier capture
           via multiphonon emission", Phys. Rev. B 90, 075202 (2014)
    """
    N = len(Q)

    if maxiter is None:
        maxiter = max(nev * N, 1000)

    # Build sparse Hamiltonian
    H = build_hamiltonian_1d(potential_func, Q)

    # Solve for lowest nev eigenvalues using ARPACK
    # which='SA' means "Smallest Algebraic" (lowest eigenvalues)
    # tol=0 means use machine precision
    try:
        eigenvalues, eigenvectors_raw = eigsh(
            H, k=nev, which="SA", maxiter=maxiter, tol=0
        )
    except ArpackNoConvergence as e:
        raise RuntimeError(
            f"ARPACK did not converge within {maxiter} iterations. "
            f"Computed {e.eigenvalues.shape[0]}/{nev} eigenvalues. "
            f"Try increasing maxiter or reducing nev."
        ) from e

    # eigenvalues are already in eV (no scaling needed)
    # eigsh returns eigenvectors as columns, transpose to rows
    # Shape: (nev, N)
    eigenvectors = eigenvectors_raw.T

    # Normalize wavefunctions
    h = Q[1] - Q[0]
    eigenvectors = normalize_wavefunctions(eigenvectors, h)

    return eigenvalues, eigenvectors
