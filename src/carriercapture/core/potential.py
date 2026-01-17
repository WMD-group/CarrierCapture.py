"""
Potential energy surface fitting and quantum state solving.

This module provides the Potential class for representing 1D potential energy surfaces,
fitting functions to data, and solving the Schrödinger equation for phonon states.
"""

from typing import Optional, Callable, Dict, Any, Literal
from pathlib import Path
import numpy as np
from numpy.typing import NDArray
import copy

from .schrodinger import solve_schrodinger_1d
from .._constants import AMU, HBAR_C, K_B

# Type alias for fit types
FitType = Literal["spline", "bspline", "harmonic", "polyfunc", "morse", "morse_poly"]


class Potential:
    """
    Represents a 1D potential energy surface with fitting and quantum solving capabilities.

    This class stores potential energy surface data, fits analytical or interpolated
    functions to the data, and solves the 1D Schrödinger equation to obtain phonon
    eigenvalues and wavefunctions.

    Attributes
    ----------
    name : str
        Identifier for the potential
    Q_data : NDArray[np.float64] | None
        Configuration coordinates (sample points), amu^0.5·Å
    E_data : NDArray[np.float64] | None
        Energy values at sample points, eV
    Q0 : float
        Minimum configuration coordinate, amu^0.5·Å
    E0 : float
        Minimum energy, eV
    Q : NDArray[np.float64] | None
        Interpolated Q grid for evaluation
    E : NDArray[np.float64] | None
        Energy on Q grid, eV
    fit_type : FitType | None
        Type of fitting function
    fit_func : Callable | None
        Fitted function (Q -> E)
    fit_params : Dict[str, Any]
        Hyperparameters for fitting
    eigenvalues : NDArray[np.float64] | None
        Phonon eigenvalues (eV)
    eigenvectors : NDArray[np.float64] | None
        Phonon wavefunctions, shape (nev, len(Q))
    nev : int
        Number of eigenvalues to compute
    temperature : float
        Temperature for thermal weighting (K)
    use_thermal_weights : bool
        Enable Boltzmann weighting in fitting

    Examples
    --------
    Create from file and fit:

    >>> pot = Potential.from_file("data.dat", name="DX-center")
    >>> pot.fit(fit_type="spline", order=4, smoothness=0.001)
    >>> pot.solve(nev=60)
    >>> pot.eigenvalues[:5]  # First 5 eigenvalues
    array([0.523, 0.569, 0.615, 0.661, 0.707])

    Create harmonic potential:

    >>> pot = Potential.from_harmonic(hw=0.02, Q0=0.0, E0=0.0)
    >>> pot.solve(nev=10)
    >>> pot.eigenvalues[:3]  # E_n = ℏω(n + 1/2)
    array([0.01, 0.03, 0.05])

    Evaluate potential:

    >>> pot(5.0)  # Energy at Q=5.0
    1.234
    """

    def __init__(
        self,
        name: str = "",
        Q0: float = 0.0,
        E0: float = 0.0,
        nev: int = 30,
        temperature: float = 293.15,
        Q_data: Optional[NDArray[np.float64]] = None,
        E_data: Optional[NDArray[np.float64]] = None,
    ):
        """
        Initialize a Potential.

        Parameters
        ----------
        name : str, default=""
            Identifier for the potential
        Q0 : float, default=0.0
            Minimum configuration coordinate (amu^0.5·Å)
        E0 : float, default=0.0
            Minimum energy (eV)
        nev : int, default=30
            Number of eigenvalues to compute
        temperature : float, default=293.15
            Temperature for thermal weighting (K)
        Q_data : NDArray[np.float64], optional
            Configuration coordinate data points (amu^0.5·Å)
        E_data : NDArray[np.float64], optional
            Potential energy data points (eV)
        """
        self.name = name
        self.Q_data = Q_data
        self.E_data = E_data
        self.Q0 = Q0
        self.E0 = E0
        self.Q: Optional[NDArray[np.float64]] = None
        self.E: Optional[NDArray[np.float64]] = None
        self.fit_type: Optional[FitType] = None
        self.fit_func: Optional[Callable] = None
        self.fit_params: Dict[str, Any] = {}
        self.eigenvalues: Optional[NDArray[np.float64]] = None
        self.eigenvectors: Optional[NDArray[np.float64]] = None
        self.nev = nev
        self.temperature = temperature
        self.use_thermal_weights = False

    @classmethod
    def from_file(
        cls,
        filepath: str | Path,
        name: str = "",
        resolution: int = 3001,
        **kwargs,
    ) -> "Potential":
        """
        Load potential from two-column data file.

        File format: Two columns (Q, E) separated by whitespace or comma.
        Lines starting with '#' are treated as comments.

        Parameters
        ----------
        filepath : str or Path
            Path to data file
        name : str, default=""
            Name for the potential
        resolution : int, default=3001
            Number of points for interpolation grid
        **kwargs
            Additional arguments passed to __init__

        Returns
        -------
        pot : Potential
            Potential loaded from file

        Examples
        --------
        >>> pot = Potential.from_file("potential.csv", name="DX-center")
        >>> pot.Q_data.shape
        (25,)
        """
        filepath = Path(filepath)

        # Load data (handle both comma and whitespace separated)
        try:
            data = np.loadtxt(filepath, comments="#", delimiter=",")
        except ValueError:
            data = np.loadtxt(filepath, comments="#")

        if data.shape[1] != 2:
            raise ValueError(f"Expected 2 columns (Q, E), got {data.shape[1]}")

        Q_data = data[:, 0]
        E_data = data[:, 1]

        # Sort by Q
        sort_idx = np.argsort(Q_data)
        Q_data = Q_data[sort_idx]
        E_data = E_data[sort_idx]

        # Create potential
        pot = cls(name=name or filepath.stem, **kwargs)
        pot.Q_data = Q_data
        pot.E_data = E_data

        # Find minimum
        min_idx = np.argmin(E_data)
        pot.Q0 = Q_data[min_idx]
        pot.E0 = E_data[min_idx]

        # Create interpolation grid
        pot.Q = np.linspace(Q_data.min(), Q_data.max(), resolution)

        return pot

    @classmethod
    def from_harmonic(
        cls,
        hw: float,
        Q0: float = 0.0,
        E0: float = 0.0,
        Q_range: tuple[float, float] = (-20.0, 20.0),
        npoints: int = 5000,
        **kwargs,
    ) -> "Potential":
        """
        Create harmonic potential: E = E0 + (1/2) * k * (Q - Q0)².

        Where k = (amu/2) * (ℏω/ℏc)²

        Parameters
        ----------
        hw : float
            Phonon energy ℏω (eV)
        Q0 : float, default=0.0
            Equilibrium position (amu^0.5·Å)
        E0 : float, default=0.0
            Minimum energy (eV)
        Q_range : tuple[float, float], default=(-20.0, 20.0)
            Range for Q grid (amu^0.5·Å)
        npoints : int, default=5000
            Number of grid points
        **kwargs
            Additional arguments passed to __init__

        Returns
        -------
        pot : Potential
            Harmonic potential

        Examples
        --------
        >>> pot = Potential.from_harmonic(hw=0.02)  # 20 meV phonon
        >>> pot.solve(nev=10)
        >>> pot.eigenvalues[0]  # Ground state E_0 = ℏω/2
        0.01
        """
        pot = cls(name=f"Harmonic (ℏω={hw:.4f} eV)", Q0=Q0, E0=E0, **kwargs)

        # Spring constant: k = (m/2) * ω² = (m/2) * (ℏω/ℏ)²
        # In our units: k has units of eV/(amu^0.5·Å)²
        # k = (AMU / 2) * (hw / HBAR)**2, but need to convert units properly
        # Actually, for harmonic: E = (1/2) * k * Q² where k = m * ω²
        # But we want it in terms of hw, so: k = m * (hw/ℏ)²
        # Let's use the simpler form: E = (1/2) * a * (Q - Q0)²
        # where a = (amu/2) * (ℏω/ℏc)² with proper unit conversion

        # Conversion: a has units eV/(amu^0.5·Å)²
        # a = (AMU / 2) * (hw / (HBAR_C * 1e10))**2
        a = (AMU / 2) * (hw / (HBAR_C * 1e10)) ** 2

        # Create grid
        pot.Q = np.linspace(Q_range[0], Q_range[1], npoints)

        # Create harmonic function
        # Note: a already includes the factor of 1/2
        def harmonic_func(Q):
            return E0 + a * (Q - Q0) ** 2

        pot.E = harmonic_func(pot.Q)
        pot.fit_func = harmonic_func
        pot.fit_type = "harmonic"
        pot.fit_params = {"hw": hw, "a": a}

        return pot

    def fit(
        self,
        fit_type: FitType = "spline",
        **fit_kwargs,
    ) -> None:
        """
        Fit function to data points.

        Parameters
        ----------
        fit_type : FitType, default="spline"
            Fitting method: "spline", "bspline", "harmonic", "polyfunc",
            "morse", "morse_poly"
        **fit_kwargs
            Fitting hyperparameters specific to each method:

            For "spline":
                - order : int, default=2
                    Spline degree (1-5)
                - smoothness : float, default=0.0
                    Smoothing parameter (0 = interpolating)
                - weights : array-like, optional
                    Weights for each data point

            For "harmonic":
                - hw : float
                    Phonon energy ℏω (eV)

        Raises
        ------
        ValueError
            If Q_data and E_data are not set
        NotImplementedError
            For fit types not yet implemented

        Examples
        --------
        >>> pot.fit(fit_type="spline", order=4, smoothness=0.001)
        >>> pot.fit(fit_type="harmonic", hw=0.02)
        """
        if self.Q_data is None or self.E_data is None:
            raise ValueError("Must set Q_data and E_data before fitting")

        if self.Q is None:
            # Create default interpolation grid
            self.Q = np.linspace(self.Q_data.min(), self.Q_data.max(), 3001)

        self.fit_type = fit_type
        self.fit_params = fit_kwargs

        if fit_type == "spline":
            self._fit_spline(**fit_kwargs)
        elif fit_type == "harmonic":
            self._fit_harmonic(**fit_kwargs)
        elif fit_type in ("polyfunc", "polynomial"):
            self._fit_polynomial(**fit_kwargs)
        elif fit_type == "morse":
            self._fit_morse(**fit_kwargs)
        elif fit_type == "morse_poly":
            self._fit_morse_poly(**fit_kwargs)
        else:
            raise NotImplementedError(f"Fit type '{fit_type}' not yet implemented")

        # Evaluate on grid
        self.E = self.fit_func(self.Q)

    def _fit_spline(
        self,
        order: int = 2,
        smoothness: float = 0.0,
        weights: Optional[NDArray[np.float64]] = None,
    ) -> None:
        """Fit smoothing spline using scipy's UnivariateSpline (FITPACK)."""
        from scipy.interpolate import UnivariateSpline

        # Sort data (required by FITPACK)
        sort_idx = np.argsort(self.Q_data)
        Q_sorted = self.Q_data[sort_idx]
        E_sorted = self.E_data[sort_idx]

        if weights is not None:
            weights = np.array(weights)[sort_idx]

        # Create spline
        # s=smoothness: if 0, interpolating spline
        # k=order: polynomial degree
        # ext='extrapolate': allow evaluation outside data range
        spl = UnivariateSpline(
            Q_sorted, E_sorted, w=weights, k=order, s=smoothness, ext="extrapolate"
        )

        self.fit_func = spl

    def _fit_harmonic(self, hw: float) -> None:
        """Fit harmonic potential to data."""
        # Use the minimum from data
        # Note: a already includes the factor of 1/2
        a = (AMU / 2) * (hw / (HBAR_C * 1e10)) ** 2

        def harmonic_func(Q):
            return self.E0 + a * (Q - self.Q0) ** 2

        self.fit_func = harmonic_func

    def _fit_polynomial(
        self,
        poly_order: Optional[int] = None,
        degree: Optional[int] = None,
        p0: Optional[list[float]] = None,
    ) -> None:
        """Fit polynomial potential to data using scipy.optimize.curve_fit."""
        from scipy.optimize import curve_fit

        # Handle both poly_order and degree parameter names
        if degree is not None:
            poly_order = degree
        elif poly_order is None:
            poly_order = 4  # Default

        # Polynomial: E = Σ coeffs[i] * (Q - Q0)^i for i=0..poly_order
        # Note: coeffs[0] = E0 (constant term)
        def poly_func(Q, *coeffs):
            result = np.zeros_like(Q)
            for i, coeff in enumerate(coeffs):
                result += coeff * (Q - self.Q0) ** i
            return result

        # Initial guess (degree + 1 coefficients)
        if p0 is None:
            p0 = np.zeros(poly_order + 1)
            p0[0] = self.E_data.min()  # Guess E0 from data

        # Fit using nonlinear least squares
        try:
            popt, _ = curve_fit(poly_func, self.Q_data, self.E_data, p0=p0)
        except RuntimeError as e:
            raise ValueError(f"Polynomial fitting failed: {e}") from e

        self.fit_func = lambda Q: poly_func(Q, *popt)
        self.fit_params["degree"] = poly_order
        self.fit_params["coeffs"] = popt.tolist()

        # Update E0 from fitted constant term
        self.E0 = float(popt[0])

    def _fit_morse(self, p0: Optional[list[float]] = None) -> None:
        """
        Fit Morse potential to data.

        Morse potential: E = E0 + D * (1 - exp(-a*(Q - Q0)))²

        Estimates Q0 and E0 from data, fits D and a.

        Parameters
        ----------
        p0 : list[float], optional
            Initial parameters [D, a] where:
            - D: Morse depth (eV)
            - a: Morse width parameter (1/amu^0.5·Å)
        """
        from scipy.optimize import curve_fit

        # Estimate Q0 and E0 from data minimum
        min_idx = np.argmin(self.E_data)
        Q0_fit = float(self.Q_data[min_idx])
        E0_fit = float(self.E_data[min_idx])

        # Filter out extreme outliers that would dominate the fit
        # Morse potentials can have huge values far from equilibrium
        E_range = self.E_data.max() - E0_fit
        reasonable_threshold = E0_fit + min(E_range, 100.0)  # Cap at 100 eV above minimum
        reasonable_mask = self.E_data <= reasonable_threshold

        Q_fit_data = self.Q_data[reasonable_mask]
        E_fit_data = self.E_data[reasonable_mask]

        def morse_func(Q, D, a):
            return E0_fit + D * (1 - np.exp(-a * (Q - Q0_fit))) ** 2

        # Initial guess
        if p0 is None:
            # Estimate D: The asymptotic energy (far from Q0) approaches E0 + D
            # Use data points far from minimum
            far_indices = np.abs(Q_fit_data - Q0_fit) > 0.5 * (Q_fit_data.max() - Q_fit_data.min())
            if np.any(far_indices):
                # Take reasonable subset of far points (not too extreme)
                E_far = E_fit_data[far_indices]
                E_far_reasonable = E_far[E_far < E0_fit + 1000]  # Exclude extreme outliers
                if len(E_far_reasonable) > 0:
                    D_guess = np.median(E_far_reasonable) - E0_fit
                else:
                    D_guess = (E_fit_data.max() - E0_fit)
            else:
                D_guess = (E_fit_data.max() - E0_fit)

            # Estimate a from curvature near minimum
            # For Morse: d²E/dQ² |_Q0 = 2*D*a²
            # Estimate from finite differences
            near_min = np.abs(Q_fit_data - Q0_fit) < 0.2 * (Q_fit_data.max() - Q_fit_data.min())
            if np.sum(near_min) > 3:
                Q_near = Q_fit_data[near_min]
                E_near = E_fit_data[near_min]
                # Fit parabola to get curvature
                coeffs = np.polyfit(Q_near - Q0_fit, E_near - E0_fit, 2)
                curvature = 2 * coeffs[0]  # Second derivative
                if D_guess > 0 and curvature > 0:
                    a_guess = np.sqrt(curvature / (2 * D_guess))
                else:
                    a_guess = 1.0
            else:
                a_guess = 1.0

            p0 = [D_guess, a_guess]

        try:
            popt, _ = curve_fit(
                morse_func,
                Q_fit_data,
                E_fit_data,
                p0=p0,
                maxfev=5000,
            )
        except RuntimeError as e:
            raise ValueError(f"Morse fitting failed: {e}") from e

        D, a = popt
        self.fit_func = lambda Q: morse_func(Q, D, a)
        self.fit_params["D"] = float(D)
        self.fit_params["a"] = float(a)
        self.fit_params["Q0"] = Q0_fit
        self.fit_params["E0"] = E0_fit

        # Update instance attributes
        self.Q0 = Q0_fit
        self.E0 = E0_fit

    def _fit_morse_poly(
        self,
        poly_order: int = 4,
        p0: Optional[list[float]] = None,
    ) -> None:
        """
        Fit Morse potential with polynomial corrections.

        E = E0 + Morse(Q) + Polynomial(Q - Q0)

        Parameters
        ----------
        poly_order : int, default=4
            Order of polynomial correction
        p0 : list[float], optional
            Initial parameters [A, a, r0, poly_coeffs...]
            where A is Morse depth, a is Morse width, r0 is offset,
            and poly_coeffs are polynomial coefficients
        """
        from scipy.optimize import curve_fit

        def morse_poly_func(Q, A, a, r0, *poly_coeffs):
            # Morse part
            morse = A * (1 - np.exp(-a * (Q - self.Q0 - r0))) ** 2

            # Polynomial part
            poly = np.zeros_like(Q)
            for i, coeff in enumerate(poly_coeffs):
                poly += coeff * (Q - self.Q0) ** (i + 1)

            return self.E0 + morse + poly

        # Initial guess
        if p0 is None:
            A_guess = (self.E_data.max() - self.E0) * 2
            a_guess = 1.0
            r0_guess = 0.0
            poly_guess = np.zeros(poly_order)
            p0 = [A_guess, a_guess, r0_guess] + poly_guess.tolist()

        try:
            popt, _ = curve_fit(morse_poly_func, self.Q_data, self.E_data, p0=p0)
        except RuntimeError as e:
            raise ValueError(f"Morse-polynomial fitting failed: {e}") from e

        A, a, r0 = popt[:3]
        poly_coeffs = popt[3:]

        self.fit_func = lambda Q: morse_poly_func(Q, A, a, r0, *poly_coeffs)
        self.fit_params["A"] = float(A)
        self.fit_params["a"] = float(a)
        self.fit_params["r0"] = float(r0)
        self.fit_params["poly_coeffs"] = poly_coeffs.tolist()

    def filter_thermally_accessible(
        self, thermal_energy: Optional[float] = None
    ) -> None:
        """
        Filter data points beyond thermal energy barrier.

        Keeps only the "island" of connected data points below the energy threshold
        that contains the minimum energy point. This is useful for removing
        high-energy data points that are thermally inaccessible.

        Parameters
        ----------
        thermal_energy : float, optional
            Thermal energy threshold (eV). If None, uses k_B * self.temperature.

        Raises
        ------
        ValueError
            If Q_data or E_data is not set
        ValueError
            If multiple minima are found

        Examples
        --------
        >>> pot.filter_thermally_accessible(thermal_energy=0.5)  # Keep points within 0.5 eV of minimum
        """
        if self.Q_data is None or self.E_data is None:
            raise ValueError("Must set Q_data and E_data before filtering")

        if thermal_energy is None:
            thermal_energy = K_B * self.temperature

        # Find minimum
        min_indices = np.where(self.E_data == self.E_data.min())[0]
        if len(min_indices) > 1:
            raise ValueError("Multiple minima found in data")

        min_idx = min_indices[0]

        # Points below threshold
        below_thresh = self.E_data <= (self.E0 + thermal_energy)

        # Find connected "island" containing minimum
        island = []
        for i, is_below in enumerate(below_thresh):
            if not is_below:
                # Check if we've found the island containing minimum
                if min_idx in island:
                    break
                # Reset island if we haven't reached minimum yet
                island = []
            else:
                island.append(i)

        # Filter data
        island = np.array(island)
        self.Q_data = self.Q_data[island]
        self.E_data = self.E_data[island]

    def solve(
        self,
        nev: Optional[int] = None,
        maxiter: Optional[int] = None,
    ) -> None:
        """
        Solve 1D Schrödinger equation for this potential.

        Updates self.eigenvalues and self.eigenvectors in place.

        Parameters
        ----------
        nev : int, optional
            Number of eigenvalues to compute. If None, uses self.nev.
        maxiter : int, optional
            Maximum ARPACK iterations. If None, uses default.

        Raises
        ------
        ValueError
            If potential function is not fitted
        RuntimeError
            If ARPACK fails to converge

        Examples
        --------
        >>> pot.fit(fit_type="spline")
        >>> pot.solve(nev=60)
        >>> len(pot.eigenvalues)
        60
        """
        if self.fit_func is None:
            raise ValueError("Must fit potential before solving")

        if nev is not None:
            self.nev = nev

        # Solve Schrödinger equation
        self.eigenvalues, self.eigenvectors = solve_schrodinger_1d(
            self.fit_func, self.Q, nev=self.nev, maxiter=maxiter
        )

    def __call__(self, Q: float | NDArray[np.float64]) -> float | NDArray[np.float64]:
        """
        Evaluate potential at Q.

        Parameters
        ----------
        Q : float or NDArray[np.float64]
            Configuration coordinate(s)

        Returns
        -------
        E : float or NDArray[np.float64]
            Energy at Q

        Raises
        ------
        ValueError
            If potential is not fitted
        """
        if self.fit_func is None:
            raise ValueError("Must fit potential before evaluation")
        return self.fit_func(Q)

    def copy(self) -> "Potential":
        """
        Create a deep copy of this potential.

        Returns
        -------
        pot_copy : Potential
            Independent copy
        """
        return copy.deepcopy(self)

    def to_dict(self) -> dict:
        """
        Serialize potential to dictionary.

        Returns
        -------
        data : dict
            Dictionary representation
        """
        data = {
            "name": self.name,
            "Q0": self.Q0,
            "E0": self.E0,
            "nev": self.nev,
            "temperature": self.temperature,
            "use_thermal_weights": self.use_thermal_weights,
            "fit_type": self.fit_type,
            "fit_params": self.fit_params,
        }

        # Add arrays if present
        if self.Q_data is not None:
            data["Q_data"] = self.Q_data.tolist()
        if self.E_data is not None:
            data["E_data"] = self.E_data.tolist()
        if self.Q is not None:
            data["Q"] = self.Q.tolist()
        if self.E is not None:
            data["E"] = self.E.tolist()
        if self.eigenvalues is not None:
            data["eigenvalues"] = self.eigenvalues.tolist()
        if self.eigenvectors is not None:
            data["eigenvectors"] = self.eigenvectors.tolist()

        return data

    @classmethod
    def from_dict(cls, data: dict) -> "Potential":
        """
        Deserialize potential from dictionary.

        Parameters
        ----------
        data : dict
            Dictionary representation

        Returns
        -------
        pot : Potential
            Reconstructed potential
        """
        pot = cls(
            name=data.get("name", ""),
            Q0=data.get("Q0", 0.0),
            E0=data.get("E0", 0.0),
            nev=data.get("nev", 30),
            temperature=data.get("temperature", 293.15),
        )

        pot.use_thermal_weights = data.get("use_thermal_weights", False)
        pot.fit_type = data.get("fit_type")
        pot.fit_params = data.get("fit_params", {})

        # Restore arrays
        if "Q_data" in data:
            pot.Q_data = np.array(data["Q_data"])
        if "E_data" in data:
            pot.E_data = np.array(data["E_data"])
        if "Q" in data:
            pot.Q = np.array(data["Q"])
        if "E" in data:
            pot.E = np.array(data["E"])
        if "eigenvalues" in data:
            pot.eigenvalues = np.array(data["eigenvalues"])
        if "eigenvectors" in data:
            pot.eigenvectors = np.array(data["eigenvectors"])

        # Re-fit to restore function
        if pot.fit_type and pot.Q_data is not None:
            pot.fit(pot.fit_type, **pot.fit_params)

        return pot


def find_crossing(pot1: Potential, pot2: Potential) -> tuple[float, float]:
    """
    Find crossing point between two potential energy surfaces.

    Finds the configuration coordinate Q where pot1(Q) = pot2(Q).

    Parameters
    ----------
    pot1 : Potential
        First potential
    pot2 : Potential
        Second potential

    Returns
    -------
    Q_cross : float
        Configuration coordinate at crossing (amu^0.5·Å)
    E_cross : float
        Energy at crossing (eV)

    Raises
    ------
    ValueError
        If potentials are not fitted
    RuntimeError
        If no crossing point is found

    Examples
    --------
    >>> pot1 = Potential.from_harmonic(hw=0.03, Q0=0, E0=1.0)
    >>> pot2 = Potential.from_harmonic(hw=0.02, Q0=10, E0=0.0)
    >>> Q_cross, E_cross = find_crossing(pot1, pot2)
    """
    from scipy.optimize import brentq

    if pot1.fit_func is None or pot2.fit_func is None:
        raise ValueError("Both potentials must be fitted")

    # Define difference function
    def diff_func(Q):
        return pot1(Q) - pot2(Q)

    # Find search range from Q grids
    if pot1.Q is not None and pot2.Q is not None:
        Q_min = max(pot1.Q.min(), pot2.Q.min())
        Q_max = min(pot1.Q.max(), pot2.Q.max())
    else:
        # Fallback: use range around Q0 values
        Q_min = min(pot1.Q0, pot2.Q0) - 20
        Q_max = max(pot1.Q0, pot2.Q0) + 20

    # Start search at midpoint
    Q_mid = (Q_min + Q_max) / 2

    # Check if there's a sign change
    try:
        Q_cross = brentq(diff_func, Q_min, Q_max)
    except ValueError:
        # No sign change found, try from midpoint
        try:
            if diff_func(Q_mid) * diff_func(Q_max) < 0:
                Q_cross = brentq(diff_func, Q_mid, Q_max)
            else:
                Q_cross = brentq(diff_func, Q_min, Q_mid)
        except ValueError as e:
            raise RuntimeError(f"No crossing point found between potentials: {e}") from e

    E_cross = pot1(Q_cross)

    return Q_cross, E_cross


def fit_morse(
    Q_data: NDArray[np.float64],
    E_data: NDArray[np.float64],
) -> tuple[dict, Callable]:
    """
    Fit Morse potential to data.

    Parameters
    ----------
    Q_data : NDArray[np.float64]
        Configuration coordinates (amu^0.5·Å)
    E_data : NDArray[np.float64]
        Potential energies (eV)

    Returns
    -------
    params : dict
        Fitted parameters: D, a, Q0, E0
    func : Callable
        Fitted Morse function

    Examples
    --------
    >>> params, func = fit_morse(Q_data, E_data)
    >>> E_fit = func(Q_data)
    """
    pot = Potential(Q_data=Q_data, E_data=E_data)
    pot.fit(fit_type="morse")
    return pot.fit_params, pot.fit_func


def fit_polynomial(
    Q_data: NDArray[np.float64],
    E_data: NDArray[np.float64],
    degree: int = 4,
) -> tuple[dict, Callable]:
    """
    Fit polynomial potential to data.

    Parameters
    ----------
    Q_data : NDArray[np.float64]
        Configuration coordinates (amu^0.5·Å)
    E_data : NDArray[np.float64]
        Potential energies (eV)
    degree : int, default=4
        Polynomial degree

    Returns
    -------
    params : dict
        Fitted parameters: degree, coeffs
    func : Callable
        Fitted polynomial function

    Examples
    --------
    >>> params, func = fit_polynomial(Q_data, E_data, degree=4)
    >>> E_fit = func(Q_data)
    """
    pot = Potential(Q_data=Q_data, E_data=E_data)
    pot.fit(fit_type="polynomial", degree=degree)
    return pot.fit_params, pot.fit_func


def fit_morse_poly(
    Q_data: NDArray[np.float64],
    E_data: NDArray[np.float64],
    poly_degree: int = 2,
) -> tuple[dict, Callable]:
    """
    Fit Morse + polynomial hybrid potential to data.

    Parameters
    ----------
    Q_data : NDArray[np.float64]
        Configuration coordinates (amu^0.5·Å)
    E_data : NDArray[np.float64]
        Potential energies (eV)
    poly_degree : int, default=2
        Polynomial correction degree

    Returns
    -------
    params : dict
        Fitted parameters: D, a, Q0, E0, poly_coeffs
    func : Callable
        Fitted hybrid function

    Examples
    --------
    >>> params, func = fit_morse_poly(Q_data, E_data, poly_degree=2)
    >>> E_fit = func(Q_data)
    """
    pot = Potential(Q_data=Q_data, E_data=E_data)
    pot.fit(fit_type="morse_poly", poly_degree=poly_degree)
    return pot.fit_params, pot.fit_func


def filter_thermally_accessible(
    Q_data: NDArray[np.float64],
    E_data: NDArray[np.float64],
    temperature: float = 300.0,
    n_kBT: float = 3.0,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """
    Filter data to thermally accessible region.

    Keeps only points within E_min + n*k_B*T of the minimum energy.

    Parameters
    ----------
    Q_data : NDArray[np.float64]
        Configuration coordinates (amu^0.5·Å)
    E_data : NDArray[np.float64]
        Potential energies (eV)
    temperature : float, default=300.0
        Temperature (K)
    n_kBT : float, default=3.0
        Number of k_B*T above minimum to keep

    Returns
    -------
    Q_filtered : NDArray[np.float64]
        Filtered configuration coordinates
    E_filtered : NDArray[np.float64]
        Filtered potential energies

    Examples
    --------
    >>> Q_filt, E_filt = filter_thermally_accessible(Q_data, E_data, temperature=300, n_kBT=3)
    """
    pot = Potential(Q_data=Q_data, E_data=E_data)
    pot.filter_thermally_accessible(temperature=temperature, n_kBT=n_kBT)
    return pot.Q_data, pot.E_data
