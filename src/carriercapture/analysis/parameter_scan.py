"""
Parameter scanning for high-throughput materials screening.

Provides tools for scanning parameter spaces (ΔQ, ΔE, ℏω) to screen
materials for carrier capture properties. Supports parallel execution
and progress reporting.
"""

from typing import Dict, Optional, Tuple, Any, Union
from pathlib import Path
import json
import numpy as np
from numpy.typing import NDArray
import warnings
from dataclasses import dataclass, field, asdict

from carriercapture.core.potential import Potential
from carriercapture.core.config_coord import ConfigCoordinate


@dataclass
class ScanParameters:
    """
    Container for parameter scan configuration.

    Attributes
    ----------
    dQ_range : tuple
        (min, max, n_points) for ΔQ grid (amu^0.5·Å)
    dE_range : tuple
        (min, max, n_points) for ΔE grid (eV)
    hbar_omega_i : float or tuple
        ℏω for initial state (eV). If tuple: (min, max, n_points)
    hbar_omega_f : float or tuple
        ℏω for final state (eV). If tuple: (min, max, n_points)
    W : float
        Electron-phonon coupling matrix element (eV/(amu^0.5·Å)).
        Must be calculated for the defect transition of interest;
        since C ∝ W², it scales the whole map uniformly.
    temperature : float or NDArray
        Temperature(s) for calculation (K)
    volume : float
        Supercell volume (cm³)
    degeneracy : int
        Degeneracy factor
    sigma : float
        Gaussian delta width (eV)
    cutoff : float
        Energy cutoff for overlaps (eV)
    """
    dQ_range: Tuple[float, float, int]
    dE_range: Tuple[float, float, int]
    hbar_omega_i: Union[float, Tuple[float, float, int]] = 0.008  # 8 meV default
    hbar_omega_f: Union[float, Tuple[float, float, int]] = 0.008
    W: Optional[float] = None
    temperature: Union[float, NDArray[np.float64]] = 300.0
    volume: float = 1e-21
    degeneracy: int = 1
    sigma: float = 0.01
    cutoff: float = 0.25
    nev_initial: int = 180
    nev_final: int = 60
    Q_grid_points: int = 5000
    Q_range_padding: float = 20.0


@dataclass
class ScanResult:
    """
    Container for parameter scan results.

    Attributes
    ----------
    dQ_grid : NDArray
        ΔQ values tested
    dE_grid : NDArray
        ΔE values tested
    capture_coefficients : NDArray
        Capture coefficient array, shape (n_dQ, n_dE) or (n_dQ, n_dE, n_temp)
    barrier_heights : NDArray
        Classical barrier heights, shape (n_dQ, n_dE)
    parameters : ScanParameters
        Parameters used for scan
    metadata : dict
        Additional metadata
    """
    dQ_grid: NDArray[np.float64]
    dE_grid: NDArray[np.float64]
    capture_coefficients: NDArray[np.float64]
    barrier_heights: NDArray[np.float64]
    parameters: ScanParameters
    metadata: Dict[str, Any] = field(default_factory=dict)

    def save(self, filepath: Union[str, Path], format: str = "npz") -> None:
        """
        Save scan results to file.

        Parameters
        ----------
        filepath : str or Path
            Output file path
        format : str, default="npz"
            File format: "npz" or "hdf5"
        """
        filepath = Path(filepath)

        params_dict = asdict(self.parameters)
        if isinstance(params_dict.get("temperature"), np.ndarray):
            params_dict["temperature"] = params_dict["temperature"].tolist()
        params_json = json.dumps(params_dict)

        if format == "npz":
            np.savez_compressed(
                filepath,
                dQ_grid=self.dQ_grid,
                dE_grid=self.dE_grid,
                capture_coefficients=self.capture_coefficients,
                barrier_heights=self.barrier_heights,
                parameters_json=params_json,
                **{f"param_{k}": v for k, v in self.metadata.items()}
            )
        elif format == "hdf5":
            try:
                import h5py
                with h5py.File(filepath, 'w') as f:
                    f.create_dataset('dQ_grid', data=self.dQ_grid)
                    f.create_dataset('dE_grid', data=self.dE_grid)
                    f.create_dataset('capture_coefficients', data=self.capture_coefficients)
                    f.create_dataset('barrier_heights', data=self.barrier_heights)
                    f.attrs['parameters_json'] = params_json
                    for k, v in self.metadata.items():
                        f.attrs[k] = v
            except ImportError:
                raise ImportError("h5py not installed. Install with: pip install h5py")
        else:
            raise ValueError(f"Unknown format: {format}. Use 'npz' or 'hdf5'")

    @classmethod
    def load(cls, filepath: Union[str, Path], format: str = "npz") -> "ScanResult":
        """
        Load scan results from file.

        Parameters
        ----------
        filepath : str or Path
            Input file path
        format : str, default="npz"
            File format: "npz" or "hdf5"

        Returns
        -------
        ScanResult
            Loaded scan results
        """
        filepath = Path(filepath)

        if format == "npz":
            data = np.load(filepath, allow_pickle=True)
            metadata = {k.replace('param_', ''): v for k, v in data.items()
                       if k.startswith('param_')}
            params_json = str(data['parameters_json']) if 'parameters_json' in data else None

            return cls(
                dQ_grid=data['dQ_grid'],
                dE_grid=data['dE_grid'],
                capture_coefficients=data['capture_coefficients'],
                barrier_heights=data['barrier_heights'],
                parameters=cls._params_from_json(params_json),
                metadata=dict(metadata)
            )
        elif format == "hdf5":
            try:
                import h5py
                with h5py.File(filepath, 'r') as f:
                    metadata = dict(f.attrs)
                    params_json = metadata.pop('parameters_json', None)
                    return cls(
                        dQ_grid=f['dQ_grid'][:],
                        dE_grid=f['dE_grid'][:],
                        capture_coefficients=f['capture_coefficients'][:],
                        barrier_heights=f['barrier_heights'][:],
                        parameters=cls._params_from_json(params_json),
                        metadata=metadata
                    )
            except ImportError:
                raise ImportError("h5py not installed. Install with: pip install h5py")
        else:
            raise ValueError(f"Unknown format: {format}. Use 'npz' or 'hdf5'")

    @staticmethod
    def _params_from_json(params_json: Optional[str]) -> ScanParameters:
        """Reconstruct ScanParameters from the JSON stored by save()."""
        if not params_json:
            # File predates parameter persistence
            return ScanParameters(dQ_range=(0, 0, 0), dE_range=(0, 0, 0))
        d = json.loads(params_json)
        for key in ("dQ_range", "dE_range", "hbar_omega_i", "hbar_omega_f"):
            if isinstance(d.get(key), list):
                d[key] = tuple(d[key])
        if isinstance(d.get("temperature"), list):
            d["temperature"] = np.array(d["temperature"])
        return ScanParameters(**d)


class ParameterScanner:
    """
    High-throughput parameter scanner for materials screening.

    Scans parameter space (ΔQ, ΔE, ℏω) to compute capture coefficients
    for many material configurations. Supports parallel execution.

    Parameters
    ----------
    params : ScanParameters
        Scan configuration
    verbose : bool, default=True
        Print progress information

    Examples
    --------
    >>> params = ScanParameters(
    ...     dQ_range=(0, 25, 25),
    ...     dE_range=(0, 2.5, 10),
    ...     W=0.05,
    ... )
    >>> scanner = ParameterScanner(params)
    >>> results = scanner.run_harmonic_scan(n_jobs=4)
    >>> results.save("scan_results.npz")
    """

    def __init__(self, params: ScanParameters, verbose: bool = True):
        if params.W is None:
            raise ValueError(
                "ScanParameters.W (electron-phonon coupling, eV/(amu^0.5·Å)) "
                "must be set: it must be calculated for the defect transition "
                "of interest, and C ∝ W² scales the whole scan."
            )
        self.params = params
        self.verbose = verbose

        # Build parameter grids
        self._build_grids()

    def _build_grids(self) -> None:
        """Build parameter grids from ranges."""
        # ΔQ grid
        self.dQ_grid = np.linspace(
            self.params.dQ_range[0],
            self.params.dQ_range[1],
            self.params.dQ_range[2]
        )

        # ΔE grid
        self.dE_grid = np.linspace(
            self.params.dE_range[0],
            self.params.dE_range[1],
            self.params.dE_range[2]
        )

        # ℏω grids (if scanning)
        if isinstance(self.params.hbar_omega_i, tuple):
            self.hbar_omega_i_grid = np.linspace(*self.params.hbar_omega_i)
        else:
            self.hbar_omega_i_grid = np.array([self.params.hbar_omega_i])

        if isinstance(self.params.hbar_omega_f, tuple):
            self.hbar_omega_f_grid = np.linspace(*self.params.hbar_omega_f)
        else:
            self.hbar_omega_f_grid = np.array([self.params.hbar_omega_f])

    def _create_harmonic_potentials(
        self,
        hbar_omega_i: float,
        hbar_omega_f: float,
        dQ: float,
        dE: float
    ) -> Tuple[Potential, Potential]:
        """
        Create harmonic potentials for given parameters.

        Parameters
        ----------
        hbar_omega_i : float
            ℏω for initial state (eV)
        hbar_omega_f : float
            ℏω for final state (eV)
        dQ : float
            Horizontal shift (amu^0.5·Å)
        dE : float
            Vertical shift (eV)

        Returns
        -------
        pot_i, pot_f : Potential
            Initial and final potentials
        """
        # Create Q grid
        Q = np.linspace(
            -self.params.Q_range_padding - dQ,
            self.params.Q_range_padding + dQ,
            self.params.Q_grid_points
        )

        # Initial state (excited): centered at Q=0, shifted up by dE
        pot_i = Potential.from_harmonic(
            hw=hbar_omega_i,
            Q0=0.0,
            E0=dE,
            Q_range=(Q[0], Q[-1]),
            npoints=len(Q)
        )

        # Final state (relaxed): centered at Q=dQ, no vertical shift
        pot_f = Potential.from_harmonic(
            hw=hbar_omega_f,
            Q0=dQ,
            E0=0.0,
            Q_range=(Q[0], Q[-1]),
            npoints=len(Q)
        )

        # Solve both potentials
        pot_i.solve(nev=self.params.nev_initial)
        pot_f.solve(nev=self.params.nev_final)

        return pot_i, pot_f

    def _calculate_single_point(
        self,
        hbar_omega_i: float,
        hbar_omega_f: float,
        dQ: float,
        dE: float,
    ) -> Tuple[float, float]:
        """
        Calculate capture coefficient for single parameter point.

        Parameters
        ----------
        hbar_omega_i : float
            ℏω for initial state
        hbar_omega_f : float
            ℏω for final state
        dQ : float
            Horizontal shift
        dE : float
            Vertical shift

        Returns
        -------
        capture_coeff : float
            Capture coefficient (cm³/s)
        barrier_height : float
            Classical barrier height (eV)
        """
        try:
            # Create potentials
            pot_i, pot_f = self._create_harmonic_potentials(
                hbar_omega_i, hbar_omega_f, dQ, dE
            )

            # Create ConfigCoordinate
            cc = ConfigCoordinate(
                pot_i=pot_i,
                pot_f=pot_f,
                W=self.params.W,
                degeneracy=self.params.degeneracy
            )

            # Calculate Q0 (midpoint between minima)
            Q0 = dQ / 2

            # Calculate overlap
            cc.calculate_overlap(
                Q0=Q0,
                cutoff=self.params.cutoff,
                sigma=self.params.sigma
            )

            # Calculate capture coefficient
            if isinstance(self.params.temperature, (int, float)):
                temperature = np.array([self.params.temperature])
            else:
                temperature = self.params.temperature

            cc.calculate_capture_coefficient(
                volume=self.params.volume,
                temperature=temperature
            )

            # Get capture coefficient (first temperature point)
            capture_coeff = cc.capture_coefficient[0] if len(cc.capture_coefficient) > 0 else np.nan

            # Calculate classical barrier height
            try:
                from carriercapture.core.potential import find_crossing
                crossing_Q, crossing_E = find_crossing(pot_f, pot_i)
                barrier_height = crossing_E - dE
            except (ValueError, RuntimeError) as e:
                warnings.warn(f"Could not find crossing at dQ={dQ:.2f}, dE={dE:.2f}: {e}")
                barrier_height = np.nan

            return capture_coeff, barrier_height

        except Exception as e:
            if self.verbose:
                warnings.warn(f"Error at dQ={dQ:.2f}, dE={dE:.2f}: {e}")
            return np.nan, np.nan

    def run_harmonic_scan(
        self,
        n_jobs: int = 1,
        show_progress: bool = True,
    ) -> ScanResult:
        """
        Run parameter scan with harmonic potentials.

        Parameters
        ----------
        n_jobs : int, default=1
            Number of parallel jobs. Use -1 for all cores.
        show_progress : bool, default=True
            Show progress bar

        Returns
        -------
        ScanResult
            Scan results with capture coefficients and barrier heights

        Examples
        --------
        >>> scanner = ParameterScanner(params)
        >>> results = scanner.run_harmonic_scan(n_jobs=4)
        >>> print(f"Computed {results.capture_coefficients.size} points")
        """
        n_dQ = len(self.dQ_grid)
        n_dE = len(self.dE_grid)

        # Initialize result arrays
        capture_coeffs = np.zeros((n_dQ, n_dE))
        barrier_heights = np.zeros((n_dQ, n_dE))

        if self.verbose:
            print(f"Starting parameter scan:")
            print(f"  ΔQ: {n_dQ} points from {self.dQ_grid[0]:.2f} to {self.dQ_grid[-1]:.2f}")
            print(f"  ΔE: {n_dE} points from {self.dE_grid[0]:.2f} to {self.dE_grid[-1]:.2f}")
            print(f"  Total: {n_dQ * n_dE} calculations")
            print(f"  Using {n_jobs} job(s)")

        # Get ℏω values (use first if scanning)
        hbar_omega_i = self.hbar_omega_i_grid[0]
        hbar_omega_f = self.hbar_omega_f_grid[0]

        # Prepare parameter combinations
        params_list = [
            (hbar_omega_i, hbar_omega_f, dQ, dE, i, j)
            for i, dQ in enumerate(self.dQ_grid)
            for j, dE in enumerate(self.dE_grid)
        ]

        if n_jobs == 1:
            # Serial execution
            if show_progress:
                try:
                    from rich.progress import track
                    iterator = track(params_list, description="Scanning...")
                except ImportError:
                    iterator = params_list
                    if self.verbose:
                        print("Note: Install 'rich' for progress bars")
            else:
                iterator = params_list

            for hbar_omega_i, hbar_omega_f, dQ, dE, i, j in iterator:
                capture_coeff, barrier_height = self._calculate_single_point(
                    hbar_omega_i, hbar_omega_f, dQ, dE
                )
                capture_coeffs[i, j] = capture_coeff
                barrier_heights[i, j] = barrier_height

        else:
            # Parallel execution
            try:
                from joblib import Parallel, delayed
            except ImportError:
                raise ImportError("joblib not installed. Install with: pip install joblib")

            def _calc_with_index(args):
                return self._calculate_single_point(*args[:4]) + args[4:]

            verbose = 10 if (show_progress and self.verbose) else 0
            results = Parallel(n_jobs=n_jobs, verbose=verbose)(
                delayed(_calc_with_index)(p) for p in params_list
            )

            for capture_coeff, barrier_height, i, j in results:
                capture_coeffs[i, j] = capture_coeff
                barrier_heights[i, j] = barrier_height

        if self.verbose:
            n_success = np.sum(~np.isnan(capture_coeffs))
            print(f"✓ Scan complete: {n_success}/{capture_coeffs.size} successful")

        # Create result object
        result = ScanResult(
            dQ_grid=self.dQ_grid,
            dE_grid=self.dE_grid,
            capture_coefficients=capture_coeffs,
            barrier_heights=barrier_heights,
            parameters=self.params,
            metadata={
                'hbar_omega_i': hbar_omega_i,
                'hbar_omega_f': hbar_omega_f,
                'temperature': self.params.temperature,
                'volume': self.params.volume,
            }
        )

        return result


__all__ = [
    "ScanParameters",
    "ScanResult",
    "ParameterScanner",
]
