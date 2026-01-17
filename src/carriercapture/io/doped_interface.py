"""
Interface for doped package integration.

Provides functions to load defect data from the doped package and convert
it into CarrierCapture Potential objects for carrier capture rate calculations.

The doped package (https://github.com/SMTG-Bham/doped) automates defect
calculations and provides tools for configuration coordinate diagram generation.
This module enables seamless integration between doped and CarrierCapture.
"""

from pathlib import Path
from typing import Union, Tuple, Dict, Any, List, Optional
import numpy as np
from numpy.typing import NDArray
import warnings

try:
    from monty.serialization import loadfn
    from pymatgen.core import Structure
    from pymatgen.io.vasp import Vasprun, Outcar

    MONTY_AVAILABLE = True
except ImportError:
    MONTY_AVAILABLE = False
    warnings.warn(
        "monty and/or pymatgen not available. Install with: pip install carriercapture[doped]",
        ImportWarning
    )

try:
    from doped.core import DefectEntry, Defect
    from doped.utils.configurations import orient_s2_like_s1, get_path_structures

    DOPED_AVAILABLE = True
except ImportError:
    DOPED_AVAILABLE = False
    # Only warn if monty is available (otherwise would get duplicate warning)
    if MONTY_AVAILABLE:
        warnings.warn(
            "doped package not available. Install with: pip install carriercapture[doped]",
            ImportWarning
        )


def _check_doped_available() -> None:
    """Check if doped and dependencies are available."""
    if not DOPED_AVAILABLE or not MONTY_AVAILABLE:
        raise ImportError(
            "doped integration requires doped package and dependencies. "
            "Install with: pip install carriercapture[doped]"
        )


def load_defect_entry(file_path: Union[str, Path]) -> Any:
    """
    Load DefectEntry from doped JSON.GZ file.

    Parameters
    ----------
    file_path : str or Path
        Path to DefectEntry JSON.GZ file saved by doped

    Returns
    -------
    DefectEntry
        Loaded defect entry object

    Raises
    ------
    ImportError
        If doped package is not installed
    FileNotFoundError
        If file does not exist

    Examples
    --------
    >>> defect = load_defect_entry("v_O_defect.json.gz")
    >>> print(defect.name)
    'v_O_0'
    """
    _check_doped_available()

    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"DefectEntry file not found: {file_path}")

    defect_entry = loadfn(file_path)

    # Validate that we loaded a DefectEntry
    if not isinstance(defect_entry, DefectEntry):
        raise ValueError(
            f"File does not contain a DefectEntry object. "
            f"Got type: {type(defect_entry).__name__}"
        )

    return defect_entry


def get_available_charge_states(defect_entry: Any) -> List[int]:
    """
    Get list of available charge states from DefectEntry.

    Parameters
    ----------
    defect_entry : DefectEntry
        Loaded defect entry

    Returns
    -------
    List[int]
        Available charge states

    Examples
    --------
    >>> charges = get_available_charge_states(defect)
    >>> print(charges)
    [-2, -1, 0, +1, +2]
    """
    _check_doped_available()

    # DefectEntry has charge_state attribute for the single charge state it represents
    # For multiple charge states, user needs to load multiple DefectEntry files
    # or access from DefectThermodynamics object

    if hasattr(defect_entry, 'charge_state'):
        return [defect_entry.charge_state]
    elif hasattr(defect_entry, 'charge'):
        return [defect_entry.charge]
    else:
        raise AttributeError(
            "DefectEntry does not have 'charge_state' or 'charge' attribute"
        )


def validate_charge_states(
    defect_entry: Any,
    charge_initial: int,
    charge_final: int
) -> None:
    """
    Validate that requested charge states are available in DefectEntry.

    Parameters
    ----------
    defect_entry : DefectEntry
        Loaded defect entry
    charge_initial : int
        Initial charge state
    charge_final : int
        Final charge state

    Raises
    ------
    ValueError
        If requested charge states are not available
    """
    _check_doped_available()

    available_charges = get_available_charge_states(defect_entry)

    for charge, label in [(charge_initial, "initial"), (charge_final, "final")]:
        if charge not in available_charges:
            raise ValueError(
                f"{label.capitalize()} charge state {charge:+d} not available in DefectEntry. "
                f"Available charges: {[f'{c:+d}' for c in available_charges]}"
            )


def _get_dQ_from_structures(
    struct1: Structure,
    struct2: Structure
) -> float:
    """
    Calculate mass-weighted displacement between two structures.

    This reproduces doped's internal _get_dQ function for consistency.

    Parameters
    ----------
    struct1 : Structure
        Initial structure
    struct2 : Structure
        Final structure

    Returns
    -------
    float
        Mass-weighted displacement (amu^0.5·Å)
    """
    _check_doped_available()

    try:
        return np.sqrt(
            sum(
                (a.distance(b) ** 2) * a.specie.atomic_mass
                for a, b in zip(struct1, struct2, strict=False)
            )
        )
    except Exception as e:
        raise ValueError(
            f"Failed to calculate mass-weighted displacement. "
            f"Structures may not have matching compositions or atom ordering. "
            f"Error: {e}"
        )


def suggest_Q0(
    struct_initial: Structure,
    struct_final: Structure,
    align: bool = True
) -> float:
    """
    Suggest Q0 value based on structure displacement.

    Q0 is typically chosen as a point along the configuration coordinate
    where the wavefunction overlap is significant. A common choice is
    Q0 ~ 0.5 * ΔQ (midpoint between initial and final structures).

    Parameters
    ----------
    struct_initial : Structure
        Initial defect structure
    struct_final : Structure
        Final defect structure
    align : bool, default=True
        Whether to align structures before calculating displacement

    Returns
    -------
    float
        Suggested Q0 value (amu^0.5·Å)

    Examples
    --------
    >>> Q0 = suggest_Q0(struct_i, struct_f)
    >>> print(f"Suggested Q0: {Q0:.2f} amu^0.5·Å")
    """
    _check_doped_available()

    if align:
        # Align struct_final to match struct_initial
        struct_final_aligned = orient_s2_like_s1(struct_initial, struct_final)
    else:
        struct_final_aligned = struct_final

    dQ = _get_dQ_from_structures(struct_initial, struct_final_aligned)

    # Common choice: Q0 at midpoint
    Q0 = 0.5 * dQ

    return Q0


def load_path_calculations(
    path_dir: Union[str, Path],
    ref_structure: Optional[Structure] = None,
    energy_key: str = "energy",
    verbose: bool = False
) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    """
    Load completed VASP calculations along configuration coordinate path.

    Expects directory structure:
        path_dir/
            image_000/  (VASP calculation)
            image_001/
            ...
            image_NNN/

    Parameters
    ----------
    path_dir : str or Path
        Directory containing image_XXX/ subdirectories with VASP calculations
    ref_structure : Structure, optional
        Reference structure for alignment (typically initial structure).
        If None, uses structure from image_000 as reference.
    energy_key : str, default="energy"
        How to extract energy: "energy" (electronic), "energy_per_atom", or "free_energy"
    verbose : bool, default=False
        Print progress information

    Returns
    -------
    Q : NDArray[np.float64]
        Configuration coordinates (amu^0.5·Å), shape (n_images,)
    E : NDArray[np.float64]
        Potential energies (eV), shape (n_images,)

    Raises
    ------
    FileNotFoundError
        If path directory or image subdirectories not found
    ValueError
        If VASP calculations failed or are incomplete

    Examples
    --------
    >>> Q, E = load_path_calculations("cc_path/")
    >>> print(f"Loaded {len(Q)} images")
    >>> print(f"ΔQ = {Q[-1]:.2f} amu^0.5·Å")
    """
    _check_doped_available()

    path_dir = Path(path_dir)
    if not path_dir.exists():
        raise FileNotFoundError(f"Path directory not found: {path_dir}")

    # Find all image directories
    image_dirs = sorted(path_dir.glob("image_*"))
    if not image_dirs:
        raise FileNotFoundError(
            f"No image_XXX directories found in {path_dir}. "
            f"Expected directory structure: image_000/, image_001/, ..."
        )

    if verbose:
        print(f"Found {len(image_dirs)} image directories")

    # Load structures and energies
    structures = []
    energies = []

    for img_dir in image_dirs:
        # Try to load vasprun.xml
        vasprun_path = img_dir / "vasprun.xml"
        if not vasprun_path.exists():
            raise FileNotFoundError(
                f"vasprun.xml not found in {img_dir}. "
                f"Ensure VASP calculations are complete."
            )

        try:
            vasprun = Vasprun(vasprun_path, parse_potcar_file=False)

            # Check if calculation converged
            if not vasprun.converged:
                warnings.warn(
                    f"VASP calculation in {img_dir} did not converge! "
                    f"Results may be unreliable."
                )

            # Extract structure and energy
            structure = vasprun.final_structure

            if energy_key == "energy":
                energy = vasprun.final_energy
            elif energy_key == "energy_per_atom":
                energy = vasprun.final_energy / len(structure)
            elif energy_key == "free_energy":
                # Free energy from OUTCAR
                outcar_path = img_dir / "OUTCAR"
                if outcar_path.exists():
                    outcar = Outcar(outcar_path)
                    energy = outcar.final_energy
                else:
                    warnings.warn(f"OUTCAR not found in {img_dir}, using electronic energy")
                    energy = vasprun.final_energy
            else:
                raise ValueError(
                    f"Unknown energy_key: {energy_key}. "
                    f"Must be 'energy', 'energy_per_atom', or 'free_energy'"
                )

            structures.append(structure)
            energies.append(energy)

        except Exception as e:
            raise ValueError(
                f"Failed to load VASP calculation from {img_dir}: {e}"
            )

    if verbose:
        print(f"Loaded {len(structures)} structures with energies")

    # Use first structure as reference if not provided
    if ref_structure is None:
        ref_structure = structures[0]
        if verbose:
            print("Using first structure (image_000) as reference")

    # Calculate Q values (mass-weighted displacements from reference)
    Q_values = []
    for struct in structures:
        # Align to reference
        struct_aligned = orient_s2_like_s1(ref_structure, struct, verbose=False)
        dQ = _get_dQ_from_structures(ref_structure, struct_aligned)
        Q_values.append(dQ)

    Q = np.array(Q_values)
    E = np.array(energies)

    # Shift energies so minimum is at 0
    E = E - E.min()

    if verbose:
        print(f"Q range: {Q.min():.2f} to {Q.max():.2f} amu^0.5·Å")
        print(f"E range: {E.min():.4f} to {E.max():.4f} eV")

    return Q, E


def extract_cc_data_from_structures(
    struct_initial: Structure,
    struct_final: Structure,
    energy_initial: float,
    energy_final: float,
    n_images: int = 10,
    align: bool = True,
    verbose: bool = False
) -> Tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """
    Extract Q-E data for configuration coordinate diagram from two structures.

    Generates linear interpolation path between initial and final structures,
    with energies linearly interpolated. This is useful for quick estimates or
    when full NEB/path calculations are not available.

    Parameters
    ----------
    struct_initial : Structure
        Initial defect structure (e.g., charge state 0)
    struct_final : Structure
        Final defect structure (e.g., charge state +1)
    energy_initial : float
        Total energy of initial structure (eV)
    energy_final : float
        Total energy of final structure (eV)
    n_images : int, default=10
        Number of interpolation points
    align : bool, default=True
        Whether to align structures before interpolation
    verbose : bool, default=False
        Print information about displacement and energy

    Returns
    -------
    Q_initial : NDArray[np.float64]
        Configuration coordinates for initial potential (amu^0.5·Å)
    E_initial : NDArray[np.float64]
        Energies for initial potential (eV), harmonic approximation around initial
    E_final : NDArray[np.float64]
        Energies for final potential (eV), harmonic approximation around final

    Notes
    -----
    This function provides a simple linear interpolation. For accurate results,
    use NEB calculations or single-point calculations along the path.

    Examples
    --------
    >>> Q_i, E_i, E_f = extract_cc_data_from_structures(
    ...     struct_0, struct_1, energy_0, energy_1, n_images=15
    ... )
    """
    _check_doped_available()

    if align:
        struct_final_aligned = orient_s2_like_s1(struct_initial, struct_final, verbose=verbose)
    else:
        struct_final_aligned = struct_final

    # Calculate total displacement
    dQ_total = _get_dQ_from_structures(struct_initial, struct_final_aligned)

    if verbose:
        print(f"Total displacement: ΔQ = {dQ_total:.2f} amu^0.5·Å")
        print(f"Energy difference: ΔE = {energy_final - energy_initial:.4f} eV")

    # Generate interpolated structures using doped's function
    path_structures = get_path_structures(struct_initial, struct_final_aligned, n_images=n_images)

    # Calculate Q values for each image
    Q_values = []
    for name in sorted(path_structures.keys()):
        struct = path_structures[name]
        dQ = _get_dQ_from_structures(struct_initial, struct)
        Q_values.append(dQ)

    Q_initial = np.array(Q_values)

    # Simple harmonic approximation for energies
    # E(Q) ≈ E0 + k/2 * (Q - Q0)^2
    # For initial state: centered at Q=0
    # For final state: centered at Q=dQ_total

    # Estimate force constant from energy difference and displacement
    # k ≈ 2 * ΔE / (ΔQ)^2  (crude approximation)
    k_approx = 2 * abs(energy_final - energy_initial) / (dQ_total ** 2) if dQ_total > 0 else 1.0

    E_initial = energy_initial + 0.5 * k_approx * (Q_initial - 0) ** 2
    E_final = energy_final + 0.5 * k_approx * (Q_initial - dQ_total) ** 2

    # Shift to put minimum at 0
    E_initial = E_initial - E_initial.min()
    E_final = E_final - E_final.min()

    if verbose:
        print(f"Generated {len(Q_initial)} interpolation points")
        print(f"Initial potential: E_min = {E_initial.min():.4f} eV at Q = {Q_initial[0]:.2f}")
        print(f"Final potential: E_min = {E_final.min():.4f} eV at Q = {Q_initial[-1]:.2f}")

    return Q_initial, E_initial, E_final


def create_potential_from_doped(
    defect_entry: Any,
    charge_state: int,
    Q_data: Optional[NDArray[np.float64]] = None,
    E_data: Optional[NDArray[np.float64]] = None,
    name: Optional[str] = None,
) -> "Potential":  # type: ignore  # Forward reference to avoid circular import
    """
    Create Potential object from doped DefectEntry and Q-E data.

    Parameters
    ----------
    defect_entry : DefectEntry
        Loaded defect entry from doped
    charge_state : int
        Charge state for this potential
    Q_data : NDArray[np.float64], optional
        Configuration coordinates (amu^0.5·Å).
        If None, must be provided later before fitting.
    E_data : NDArray[np.float64], optional
        Potential energies (eV).
        If None, must be provided later before fitting.
    name : str, optional
        Name for the potential. If None, generated from defect_entry name.

    Returns
    -------
    Potential
        Potential object with Q_data and E_data set

    Examples
    --------
    >>> pot_i = create_potential_from_doped(
    ...     defect, charge_state=0, Q_data=Q, E_data=E_initial
    ... )
    >>> pot_i.fit(fit_type="spline", order=4)
    >>> pot_i.solve(nev=60)
    """
    # Import here to avoid circular import
    from carriercapture.core.potential import Potential

    _check_doped_available()

    # Validate charge state
    validate_charge_states(defect_entry, charge_state, charge_state)

    # Generate name if not provided
    if name is None:
        defect_name = getattr(defect_entry, 'name', 'defect')
        name = f"{defect_name}_q{charge_state:+d}"

    # Create Potential object
    pot = Potential(Q_data=Q_data, E_data=E_data, name=name)

    # Store metadata
    pot.metadata = {
        'source': 'doped',
        'defect_name': getattr(defect_entry, 'name', None),
        'charge_state': charge_state,
    }

    return pot


__all__ = [
    "load_defect_entry",
    "get_available_charge_states",
    "validate_charge_states",
    "suggest_Q0",
    "load_path_calculations",
    "extract_cc_data_from_structures",
    "create_potential_from_doped",
]
