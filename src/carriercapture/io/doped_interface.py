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


# =============================================================================
# Workflow Helper Functions
# =============================================================================
# These functions provide a cleaner interface for the doped -> CarrierCapture
# workflow, using doped's CCD path generation functions properly.


def prepare_ccd_structures(
    defect_entry_initial: Any,
    defect_entry_final: Any,
    verbose: bool = False,
) -> Dict[str, Any]:
    """
    Prepare structures for CCD calculation from two DefectEntry objects.

    Uses doped's orient_s2_like_s1() to ensure structures are properly
    aligned for the shortest linear path between charge states, handling
    symmetry-equivalent configurations automatically.

    Parameters
    ----------
    defect_entry_initial : DefectEntry
        DefectEntry for the initial charge state (e.g., excited state)
    defect_entry_final : DefectEntry
        DefectEntry for the final charge state (e.g., ground state)
    verbose : bool, default=False
        Print alignment information and dQ values

    Returns
    -------
    dict
        Dictionary containing:
        - 'struct_initial': Structure - Initial state structure
        - 'struct_final': Structure - Final state structure (aligned)
        - 'struct_final_original': Structure - Original final structure (before alignment)
        - 'dQ': float - Mass-weighted displacement (amu^0.5*Angstrom)
        - 'energy_initial': float - Total energy of initial state (eV)
        - 'energy_final': float - Total energy of final state (eV)
        - 'dE': float - Energy difference E_final - E_initial (eV)
        - 'charge_initial': int - Initial charge state
        - 'charge_final': int - Final charge state
        - 'defect_name': str - Name of the defect

    Raises
    ------
    ValueError
        If defect entries are for different defects (by name)
    ImportError
        If doped package is not installed

    Examples
    --------
    >>> entry_0 = load_defect_entry("v_O_0.json.gz")
    >>> entry_1 = load_defect_entry("v_O_+1.json.gz")
    >>> ccd_data = prepare_ccd_structures(entry_0, entry_1, verbose=True)
    >>> print(f"dQ = {ccd_data['dQ']:.2f} amu^0.5*A")
    >>> print(f"dE = {ccd_data['dE']:.3f} eV")

    Notes
    -----
    The alignment process uses doped's orient_s2_like_s1() which:
    - Finds the symmetry-equivalent orientation of struct_final that
      minimizes the linear path distance to struct_initial
    - Ensures atomic indices match for proper displacement calculation
    - Returns the shortest-path configuration for NEB/CCD calculations
    """
    _check_doped_available()

    # Extract structures from DefectEntry objects
    # DefectEntry stores relaxed structure in sc_entry.structure
    if hasattr(defect_entry_initial, 'sc_entry') and hasattr(defect_entry_initial.sc_entry, 'structure'):
        struct_initial = defect_entry_initial.sc_entry.structure
    elif hasattr(defect_entry_initial, 'structure'):
        struct_initial = defect_entry_initial.structure
    else:
        raise ValueError(
            "Cannot extract structure from initial DefectEntry. "
            "Expected 'sc_entry.structure' or 'structure' attribute."
        )

    if hasattr(defect_entry_final, 'sc_entry') and hasattr(defect_entry_final.sc_entry, 'structure'):
        struct_final_original = defect_entry_final.sc_entry.structure
    elif hasattr(defect_entry_final, 'structure'):
        struct_final_original = defect_entry_final.structure
    else:
        raise ValueError(
            "Cannot extract structure from final DefectEntry. "
            "Expected 'sc_entry.structure' or 'structure' attribute."
        )

    # Extract energies
    if hasattr(defect_entry_initial, 'sc_entry') and hasattr(defect_entry_initial.sc_entry, 'energy'):
        energy_initial = defect_entry_initial.sc_entry.energy
    else:
        energy_initial = None
        if verbose:
            print("Warning: Could not extract energy from initial DefectEntry")

    if hasattr(defect_entry_final, 'sc_entry') and hasattr(defect_entry_final.sc_entry, 'energy'):
        energy_final = defect_entry_final.sc_entry.energy
    else:
        energy_final = None
        if verbose:
            print("Warning: Could not extract energy from final DefectEntry")

    # Extract charge states
    charge_initial = get_available_charge_states(defect_entry_initial)[0]
    charge_final = get_available_charge_states(defect_entry_final)[0]

    # Extract defect names and validate
    defect_name_initial = getattr(defect_entry_initial, 'name', None)
    defect_name_final = getattr(defect_entry_final, 'name', None)

    # Remove charge suffix for comparison (e.g., "v_O_0" -> "v_O")
    if defect_name_initial and defect_name_final:
        base_name_initial = defect_name_initial.rsplit('_', 1)[0] if '_' in defect_name_initial else defect_name_initial
        base_name_final = defect_name_final.rsplit('_', 1)[0] if '_' in defect_name_final else defect_name_final
        if base_name_initial != base_name_final:
            raise ValueError(
                f"DefectEntry objects appear to be for different defects: "
                f"'{defect_name_initial}' vs '{defect_name_final}'"
            )
        defect_name = base_name_initial
    else:
        defect_name = defect_name_initial or defect_name_final or "unknown"

    # Align final structure to initial using doped's orient_s2_like_s1
    if verbose:
        print(f"Aligning structures for {defect_name} ({charge_initial:+d} -> {charge_final:+d})...")

    struct_final = orient_s2_like_s1(struct_initial, struct_final_original, verbose=verbose)

    # Calculate mass-weighted displacement
    dQ = _get_dQ_from_structures(struct_initial, struct_final)
    dQ_original = _get_dQ_from_structures(struct_initial, struct_final_original)

    if verbose:
        print(f"dQ (aligned): {dQ:.4f} amu^0.5*A")
        if abs(dQ - dQ_original) > 0.01:
            print(f"dQ (original): {dQ_original:.4f} amu^0.5*A (alignment reduced dQ by {dQ_original - dQ:.4f})")
        if energy_initial is not None and energy_final is not None:
            print(f"dE: {energy_final - energy_initial:.4f} eV")

    result = {
        'struct_initial': struct_initial,
        'struct_final': struct_final,
        'struct_final_original': struct_final_original,
        'dQ': dQ,
        'energy_initial': energy_initial,
        'energy_final': energy_final,
        'dE': (energy_final - energy_initial) if (energy_initial is not None and energy_final is not None) else None,
        'charge_initial': charge_initial,
        'charge_final': charge_final,
        'defect_name': defect_name,
    }

    return result


def generate_ccd_path(
    ccd_data: Dict[str, Any],
    n_images: int = 11,
    displacements: Optional[NDArray[np.float64]] = None,
    output_dir: Optional[Union[str, Path]] = None,
    write_vasp: bool = False,
) -> Dict[str, Any]:
    """
    Generate interpolated structures along the configuration coordinate path.

    Uses doped's get_path_structures() to create linearly interpolated
    structures between initial and final states for single-point calculations.

    Parameters
    ----------
    ccd_data : dict
        Output from prepare_ccd_structures() containing aligned structures
    n_images : int, default=11
        Number of images along the path (including endpoints).
        Recommended: odd number so midpoint is included.
    displacements : NDArray[np.float64], optional
        Explicit fractional displacements (0.0 to 1.0) for path generation.
        If provided, overrides n_images. Example: np.linspace(0, 1, 11)
    output_dir : str or Path, optional
        Directory to write VASP input files. If None, only returns structures.
    write_vasp : bool, default=False
        Whether to write VASP POSCAR files to output_dir.
        Requires output_dir to be set.

    Returns
    -------
    dict
        Dictionary containing:
        - 'path_structures': Dict[float, Structure] - Interpolated structures
          keyed by fractional displacement (0.0 to 1.0)
        - 'Q_fractions': NDArray - Fractional positions along path (0 to 1)
        - 'Q_values': NDArray - Actual Q values (amu^0.5*A), 0 to dQ
        - 'n_images': int - Number of images generated
        - 'output_dir': Path or None - Where files were written (if any)

    Examples
    --------
    >>> path_data = generate_ccd_path(ccd_data, n_images=11)
    >>> print(f"Generated {path_data['n_images']} structures")
    >>> print(f"Q range: 0 to {path_data['Q_values'][-1]:.2f}")

    >>> # Write VASP inputs for single-point calculations
    >>> path_data = generate_ccd_path(
    ...     ccd_data, n_images=11,
    ...     output_dir="ccd_path/", write_vasp=True
    ... )

    Notes
    -----
    For accurate CCD calculations, single-point DFT calculations should be
    performed at each interpolated structure. The resulting energies can be
    loaded using load_path_calculations() and used to create Potential objects.
    """
    _check_doped_available()

    struct_initial = ccd_data['struct_initial']
    struct_final = ccd_data['struct_final']
    dQ = ccd_data['dQ']

    # Determine displacements
    if displacements is None:
        displacements = np.linspace(0, 1, n_images)
    else:
        displacements = np.asarray(displacements)
        n_images = len(displacements)

    # Generate path structures using doped
    # get_path_structures returns dict keyed by displacement fraction
    path_structures_raw = get_path_structures(
        struct_initial, struct_final,
        displacements=displacements
    )

    # Convert to more usable format
    path_structures = {}
    Q_fractions = []
    Q_values = []

    for frac in sorted(path_structures_raw.keys()):
        path_structures[frac] = path_structures_raw[frac]
        Q_fractions.append(frac)
        Q_values.append(frac * dQ)

    Q_fractions = np.array(Q_fractions)
    Q_values = np.array(Q_values)

    # Write VASP files if requested
    output_path = None
    if write_vasp and output_dir is not None:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        for i, frac in enumerate(sorted(path_structures.keys())):
            img_dir = output_path / f"image_{i:03d}"
            img_dir.mkdir(exist_ok=True)
            poscar_path = img_dir / "POSCAR"
            path_structures[frac].to(fmt="poscar", filename=str(poscar_path))

    result = {
        'path_structures': path_structures,
        'Q_fractions': Q_fractions,
        'Q_values': Q_values,
        'n_images': n_images,
        'output_dir': output_path,
    }

    return result


def estimate_phonon_frequency(
    Q_data: NDArray[np.float64],
    E_data: NDArray[np.float64],
    Q0: Optional[float] = None,
    method: str = "curvature",
) -> Dict[str, float]:
    """
    Estimate effective phonon frequency from Q-E data.

    The phonon frequency hw is related to the curvature of the potential
    energy surface near the minimum. This function estimates hw from DFT
    single-point calculations along the CCD path.

    Parameters
    ----------
    Q_data : NDArray[np.float64]
        Configuration coordinates (amu^0.5*Angstrom)
    E_data : NDArray[np.float64]
        Potential energies (eV)
    Q0 : float, optional
        Equilibrium position. If None, uses the position of minimum E_data.
    method : {"curvature", "harmonic_fit"}, default="curvature"
        Method for frequency estimation:
        - "curvature": Finite difference second derivative at minimum
        - "harmonic_fit": Least-squares fit of E = E0 + 0.5*k*(Q-Q0)^2

    Returns
    -------
    dict
        Dictionary containing:
        - 'hw': float - Phonon energy (eV)
        - 'hw_meV': float - Phonon energy (meV)
        - 'omega': float - Angular frequency (rad/fs)
        - 'frequency_THz': float - Frequency (THz)
        - 'Q0': float - Equilibrium position used
        - 'E0': float - Minimum energy
        - 'curvature': float - d²E/dQ² at minimum (eV/(amu*A²))
        - 'method': str - Method used for estimation

    Examples
    --------
    >>> Q, E = load_path_calculations("ccd_path/")
    >>> freq_data = estimate_phonon_frequency(Q, E)
    >>> print(f"Estimated phonon energy: {freq_data['hw_meV']:.1f} meV")

    Notes
    -----
    The relationship between curvature and phonon frequency is:
        hw = hbar * sqrt(k / m_eff)
    where k = d²E/dQ² and m_eff = 1 amu (since Q is mass-weighted).

    For accurate phonon frequencies, phonon calculations (DFPT/finite
    displacement) should be used. This estimate is useful for initial
    harmonic potential approximations.
    """
    from .._constants import HBAR

    Q_data = np.asarray(Q_data)
    E_data = np.asarray(E_data)

    # Find minimum
    min_idx = np.argmin(E_data)
    E0 = E_data[min_idx]

    if Q0 is None:
        Q0 = Q_data[min_idx]

    if method == "curvature":
        # Finite difference second derivative at minimum
        # Need at least 3 points near minimum
        if len(Q_data) < 3:
            raise ValueError("Need at least 3 data points for curvature estimation")

        # Sort by Q and interpolate
        sort_idx = np.argsort(Q_data)
        Q_sorted = Q_data[sort_idx]
        E_sorted = E_data[sort_idx]

        # Find closest points to Q0
        distances = np.abs(Q_sorted - Q0)
        closest_idx = np.argsort(distances)[:min(5, len(Q_sorted))]
        closest_idx = np.sort(closest_idx)

        if len(closest_idx) >= 3:
            # Use polynomial fit to estimate curvature
            Q_local = Q_sorted[closest_idx]
            E_local = E_sorted[closest_idx]
            coeffs = np.polyfit(Q_local - Q0, E_local, 2)
            curvature = 2 * coeffs[0]  # d²E/dQ² = 2*a for ax² + bx + c
        else:
            # Fall back to simple central difference
            dQ = Q_sorted[1] - Q_sorted[0] if len(Q_sorted) > 1 else 1.0
            if min_idx > 0 and min_idx < len(E_sorted) - 1:
                curvature = (E_sorted[min_idx + 1] - 2*E_sorted[min_idx] + E_sorted[min_idx - 1]) / (dQ**2)
            else:
                curvature = 1.0  # Default fallback

    elif method == "harmonic_fit":
        # Least-squares fit to harmonic potential
        # E = E0 + 0.5 * k * (Q - Q0)²
        # Linearize: E - E0 = 0.5 * k * (Q - Q0)²
        dQ_squared = (Q_data - Q0) ** 2
        dE = E_data - E0

        # Filter out points too far from minimum
        mask = dE < (E_data.max() - E0) * 0.5  # Use lower half of energy range
        if mask.sum() >= 2:
            dQ_sq_fit = dQ_squared[mask]
            dE_fit = dE[mask]
        else:
            dQ_sq_fit = dQ_squared
            dE_fit = dE

        # Linear fit: dE = 0.5 * k * dQ²
        # k = 2 * slope
        if len(dQ_sq_fit) > 1 and dQ_sq_fit.max() > 0:
            k_half = np.sum(dQ_sq_fit * dE_fit) / np.sum(dQ_sq_fit ** 2)
            curvature = 2 * k_half
        else:
            curvature = 1.0  # Default fallback

    else:
        raise ValueError(f"Unknown method: {method}. Must be 'curvature' or 'harmonic_fit'")

    # Ensure curvature is positive
    curvature = abs(curvature)

    # Convert curvature to phonon frequency
    # hw = hbar * sqrt(k / m_eff)
    # With Q in amu^0.5*A and E in eV:
    # k has units eV / (amu * A²)
    # m_eff = 1 amu (mass-weighted coordinate)
    #
    # hw (eV) = hbar (eV*s) * sqrt(k (eV/(amu*A²)) / m (amu))
    # Need to convert units properly

    # hbar in eV*fs, need angular frequency in rad/fs
    # k in eV/(amu*A²), convert A² to m² and amu to kg for SI
    # Actually simpler: use natural units

    # k in eV/(amu*A²) -> convert to eV/(eV/c² * A²) = c²/A²
    # omega² = k/m = k (eV/(amu*A²)) / (1 amu)
    # omega = sqrt(k) * (1/A) * sqrt(eV/amu)

    # Conversion: 1 amu = 931.5 MeV/c², 1 A = 1e-10 m
    # sqrt(eV/amu) = sqrt(1e-6 MeV / 931.5 MeV/c²) = sqrt(1.074e-9) c = 3.28e-5 c
    # In frequency: sqrt(eV/(amu*A²)) * A = sqrt(eV/amu) / A

    # Simpler approach using known conversion:
    # For harmonic oscillator: E_n = hw * (n + 0.5)
    # hw = hbar * omega = hbar * sqrt(k/m)
    #
    # With k in eV/A² and m in amu:
    # hw (eV) = 0.004136 * sqrt(k (eV/A²) / m (amu))
    # But our k is in eV/(amu*A²), so m_eff = 1:
    # hw (eV) = 0.004136 * sqrt(k (eV/(amu*A²)))

    conversion_factor = 0.004135665  # sqrt(hbar² / amu) in eV*A
    hw = conversion_factor * np.sqrt(curvature)

    # Angular frequency: omega = hw / hbar
    omega = hw / HBAR  # rad/s
    omega_rad_fs = omega * 1e-15  # rad/fs

    # Frequency in THz: f = omega / (2*pi) in Hz, then to THz
    frequency_Hz = omega / (2 * np.pi)
    frequency_THz = frequency_Hz * 1e-12

    result = {
        'hw': hw,
        'hw_meV': hw * 1000,
        'omega': omega_rad_fs,
        'frequency_THz': frequency_THz,
        'Q0': Q0,
        'E0': E0,
        'curvature': curvature,
        'method': method,
    }

    return result


def calculate_Q0_crossing(
    pot_initial: "Potential",  # type: ignore
    pot_final: "Potential",  # type: ignore
    method: str = "crossing",
    search_range: Optional[Tuple[float, float]] = None,
) -> Dict[str, float]:
    """
    Calculate optimal Q0 for overlap integral evaluation.

    Q0 determines where the (Q - Q0) operator is evaluated in the
    electron-phonon coupling matrix element. The optimal choice is
    typically near the crossing point of the two potential surfaces,
    where wavefunction overlap is maximized.

    Parameters
    ----------
    pot_initial : Potential
        Fitted initial state potential (must have fit_func or E array)
    pot_final : Potential
        Fitted final state potential (must have fit_func or E array)
    method : {"crossing", "midpoint", "minimum_barrier"}, default="crossing"
        Method for Q0 determination:
        - "crossing": Q where E_initial(Q) = E_final(Q)
        - "midpoint": Simple midpoint between Q0_initial and Q0_final
        - "minimum_barrier": Q that minimizes max(E_i, E_f) along path
    search_range : tuple[float, float], optional
        (Q_min, Q_max) range to search for crossing. If None, uses
        the overlap of the two potential Q grids.

    Returns
    -------
    dict
        Dictionary containing:
        - 'Q0': float - Recommended Q0 value (amu^0.5*A)
        - 'E_crossing': float - Energy at crossing point (eV), or None if no crossing
        - 'barrier_initial': float - Barrier from initial minimum to Q0 (eV)
        - 'barrier_final': float - Barrier from final minimum to Q0 (eV)
        - 'method': str - Method used
        - 'Q0_initial': float - Equilibrium Q of initial potential
        - 'Q0_final': float - Equilibrium Q of final potential

    Raises
    ------
    ValueError
        If no crossing point found in search range (for "crossing" method)
    ValueError
        If potentials don't have required data

    Examples
    --------
    >>> Q0_data = calculate_Q0_crossing(pot_i, pot_f, method="crossing")
    >>> print(f"Q0 = {Q0_data['Q0']:.2f} amu^0.5*A")
    >>> print(f"Barrier = {Q0_data['barrier_initial']:.3f} eV")

    Notes
    -----
    The crossing point method is physically motivated: the electron-phonon
    coupling is strongest where the electronic states are degenerate.
    However, for strongly asymmetric potentials, "minimum_barrier" may
    give better numerical convergence.
    """
    from carriercapture.core.potential import find_crossing

    # Get equilibrium positions
    Q0_initial = pot_initial.Q0
    Q0_final = pot_final.Q0
    E0_initial = pot_initial.E0
    E0_final = pot_final.E0

    # Determine search range
    if search_range is None:
        if pot_initial.Q is not None and pot_final.Q is not None:
            Q_min = max(pot_initial.Q.min(), pot_final.Q.min())
            Q_max = min(pot_initial.Q.max(), pot_final.Q.max())
        else:
            Q_min = min(Q0_initial, Q0_final) - abs(Q0_final - Q0_initial)
            Q_max = max(Q0_initial, Q0_final) + abs(Q0_final - Q0_initial)
        search_range = (Q_min, Q_max)

    E_crossing = None
    barrier_initial = None
    barrier_final = None

    if method == "crossing":
        try:
            Q0, E_crossing = find_crossing(pot_initial, pot_final)
            # Calculate barriers
            E_at_Q0_initial = pot_initial(Q0) if callable(pot_initial) else E_crossing
            E_at_Q0_final = pot_final(Q0) if callable(pot_final) else E_crossing
            barrier_initial = E_crossing - E0_initial
            barrier_final = E_crossing - E0_final
        except (ValueError, RuntimeError):
            # No crossing found, fall back to midpoint
            warnings.warn(
                "No crossing point found between potentials. "
                "Falling back to midpoint method."
            )
            Q0 = (Q0_initial + Q0_final) / 2
            method = "midpoint (fallback)"

    elif method == "midpoint":
        Q0 = (Q0_initial + Q0_final) / 2
        # Try to get energy at midpoint
        if callable(pot_initial) and callable(pot_final):
            E_initial_at_Q0 = pot_initial(Q0)
            E_final_at_Q0 = pot_final(Q0)
            barrier_initial = E_initial_at_Q0 - E0_initial
            barrier_final = E_final_at_Q0 - E0_final

    elif method == "minimum_barrier":
        # Find Q that minimizes the maximum barrier
        if pot_initial.Q is not None and pot_final.Q is not None:
            # Use common Q grid
            Q_common = np.linspace(search_range[0], search_range[1], 1000)

            # Evaluate potentials
            if callable(pot_initial) and callable(pot_final):
                E_i = np.array([pot_initial(q) for q in Q_common])
                E_f = np.array([pot_final(q) for q in Q_common])
            elif pot_initial.E is not None and pot_final.E is not None:
                E_i = np.interp(Q_common, pot_initial.Q, pot_initial.E)
                E_f = np.interp(Q_common, pot_final.Q, pot_final.E)
            else:
                raise ValueError("Potentials must have E array or be callable")

            # Maximum barrier at each Q
            max_barrier = np.maximum(E_i - E0_initial, E_f - E0_final)
            min_idx = np.argmin(max_barrier)
            Q0 = Q_common[min_idx]
            barrier_initial = E_i[min_idx] - E0_initial
            barrier_final = E_f[min_idx] - E0_final
        else:
            raise ValueError("Potentials must have Q grids for minimum_barrier method")

    else:
        raise ValueError(
            f"Unknown method: {method}. "
            f"Must be 'crossing', 'midpoint', or 'minimum_barrier'"
        )

    result = {
        'Q0': Q0,
        'E_crossing': E_crossing,
        'barrier_initial': barrier_initial,
        'barrier_final': barrier_final,
        'method': method,
        'Q0_initial': Q0_initial,
        'Q0_final': Q0_final,
    }

    return result


def create_ccd_from_defect_entries(
    defect_entry_initial: Union[Any, str, Path],
    defect_entry_final: Union[Any, str, Path],
    path_dir_initial: Optional[Union[str, Path]] = None,
    path_dir_final: Optional[Union[str, Path]] = None,
    fit_type: str = "spline",
    fit_kwargs: Optional[Dict[str, Any]] = None,
    nev_initial: int = 180,
    nev_final: int = 60,
    W: Optional[float] = None,
    degeneracy: int = 1,
    Q0_method: str = "auto",
    use_harmonic: bool = False,
    hw: Optional[float] = None,
    verbose: bool = False,
) -> Tuple[Any, Dict[str, Any]]:  # Returns (ConfigCoordinate, metadata)
    """
    Create ConfigCoordinate from two doped DefectEntry objects.

    This is the main convenience function for the doped -> CarrierCapture
    workflow. It handles structure alignment, Q-E data loading or generation,
    potential fitting, Schrodinger equation solving, and ConfigCoordinate
    creation in a single call.

    Parameters
    ----------
    defect_entry_initial : DefectEntry or str or Path
        DefectEntry for initial state, or path to JSON.GZ file
    defect_entry_final : DefectEntry or str or Path
        DefectEntry for final state, or path to JSON.GZ file
    path_dir_initial : str or Path, optional
        Directory with VASP single-point calculations for initial state.
        If None, uses harmonic approximation based on endpoint energies.
    path_dir_final : str or Path, optional
        Directory with VASP single-point calculations for final state.
        If None, uses harmonic approximation based on endpoint energies.
    fit_type : str, default="spline"
        Fitting method for potential: "spline", "harmonic", "morse", etc.
        Ignored if use_harmonic=True.
    fit_kwargs : dict, optional
        Additional arguments for potential fitting (order, smoothness, etc.)
    nev_initial : int, default=180
        Number of eigenvalues to compute for initial potential
    nev_final : int, default=60
        Number of eigenvalues to compute for final potential
    W : float, optional
        Electron-phonon coupling matrix element (eV).
        If None, must be set later before calculating capture coefficient.
    degeneracy : int, default=1
        Degeneracy factor for the capture process
    Q0_method : {"crossing", "midpoint", "auto"}, default="auto"
        Method for determining Q0. "auto" tries "crossing" first,
        falls back to "midpoint" if no crossing found.
    use_harmonic : bool, default=False
        Use simple harmonic potentials instead of fitting to path data.
        Useful for quick estimates or when path calculations unavailable.
    hw : float, optional
        Phonon energy for harmonic approximation (eV).
        If None and use_harmonic=True, estimates from structure displacement.
    verbose : bool, default=False
        Print progress information

    Returns
    -------
    cc : ConfigCoordinate
        Configured ConfigCoordinate with potentials solved and ready for
        calculate_overlap() and calculate_capture_coefficient()
    metadata : dict
        Dictionary with workflow details:
        - 'ccd_data': Output from prepare_ccd_structures()
        - 'Q0': float - Q0 value used/recommended
        - 'Q0_method': str - Method used for Q0
        - 'pot_initial': Potential - Reference to initial potential
        - 'pot_final': Potential - Reference to final potential
        - 'fit_type': str - Fitting type used
        - 'hw_estimated': float or None - Estimated phonon energy if available

    Examples
    --------
    >>> # Full workflow with path calculations
    >>> cc, meta = create_ccd_from_defect_entries(
    ...     "v_O_0.json.gz", "v_O_+1.json.gz",
    ...     path_dir_initial="ccd_path_0/",
    ...     path_dir_final="ccd_path_+1/",
    ...     fit_type="spline",
    ...     W=0.068,
    ...     verbose=True
    ... )
    >>> cc.calculate_overlap(Q0=meta['Q0'])
    >>> cc.calculate_capture_coefficient(volume=1e-21, temperature=temps)

    >>> # Quick harmonic estimate (no path calculations)
    >>> cc, meta = create_ccd_from_defect_entries(
    ...     entry_0, entry_1,
    ...     use_harmonic=True,
    ...     hw=0.008,  # 8 meV phonon
    ...     W=0.068,
    ... )

    Notes
    -----
    Workflow steps performed:
    1. Load DefectEntry objects (if paths provided)
    2. Align structures using orient_s2_like_s1()
    3. Load or generate Q-E data for both potentials
    4. Fit potentials using specified method
    5. Solve Schrodinger equation for both potentials
    6. Calculate optimal Q0
    7. Create ConfigCoordinate with all parameters

    For production calculations, path_dir_initial and path_dir_final should
    contain single-point DFT calculations. The harmonic approximation is
    useful for quick screening but may not capture anharmonic effects.
    """
    from carriercapture.core.potential import Potential
    from carriercapture.core.config_coord import ConfigCoordinate

    _check_doped_available()

    if fit_kwargs is None:
        fit_kwargs = {}

    # Step 1: Load DefectEntry objects if paths provided
    if isinstance(defect_entry_initial, (str, Path)):
        if verbose:
            print(f"Loading initial DefectEntry from {defect_entry_initial}")
        defect_entry_initial = load_defect_entry(defect_entry_initial)

    if isinstance(defect_entry_final, (str, Path)):
        if verbose:
            print(f"Loading final DefectEntry from {defect_entry_final}")
        defect_entry_final = load_defect_entry(defect_entry_final)

    # Step 2: Prepare and align structures
    if verbose:
        print("\nPreparing CCD structures...")
    ccd_data = prepare_ccd_structures(
        defect_entry_initial, defect_entry_final, verbose=verbose
    )

    dQ = ccd_data['dQ']
    defect_name = ccd_data['defect_name']
    charge_initial = ccd_data['charge_initial']
    charge_final = ccd_data['charge_final']

    # Step 3: Load or generate Q-E data
    hw_estimated = None

    if use_harmonic:
        # Use harmonic approximation
        if verbose:
            print("\nUsing harmonic approximation for potentials...")

        if hw is None:
            # Estimate hw from displacement (crude approximation)
            # Typical phonon energies: 5-50 meV
            hw = 0.010  # Default 10 meV
            if verbose:
                print(f"Using default phonon energy: {hw*1000:.1f} meV")
        hw_estimated = hw

        # Create harmonic potentials
        Q_range = (-dQ * 0.5, dQ * 1.5)

        # Initial state centered at Q=0
        E0_initial = ccd_data['dE'] if ccd_data['dE'] is not None else 0.5
        pot_initial = Potential.from_harmonic(
            hw=hw, Q0=0.0, E0=E0_initial,
            Q_range=Q_range, npoints=3000
        )
        pot_initial.name = f"{defect_name}_q{charge_initial:+d}"

        # Final state centered at Q=dQ
        pot_final = Potential.from_harmonic(
            hw=hw, Q0=dQ, E0=0.0,
            Q_range=Q_range, npoints=3000
        )
        pot_final.name = f"{defect_name}_q{charge_final:+d}"

        actual_fit_type = "harmonic"

    else:
        # Load from path calculations or use interpolation
        if path_dir_initial is not None and path_dir_final is not None:
            if verbose:
                print(f"\nLoading path calculations from {path_dir_initial} and {path_dir_final}...")

            Q_i, E_i = load_path_calculations(path_dir_initial, verbose=verbose)
            Q_f, E_f = load_path_calculations(path_dir_final, verbose=verbose)

            # Estimate phonon frequency from data
            try:
                freq_data = estimate_phonon_frequency(Q_i, E_i)
                hw_estimated = freq_data['hw']
                if verbose:
                    print(f"Estimated phonon energy: {hw_estimated*1000:.1f} meV")
            except Exception as e:
                warnings.warn(f"Could not estimate phonon frequency: {e}")

        else:
            # Use structure-based interpolation
            if verbose:
                print("\nGenerating Q-E data from structure interpolation...")

            if ccd_data['energy_initial'] is None or ccd_data['energy_final'] is None:
                raise ValueError(
                    "Cannot generate Q-E data without energies. "
                    "Either provide path_dir_initial/path_dir_final or use use_harmonic=True"
                )

            Q_i, E_i, E_f = extract_cc_data_from_structures(
                ccd_data['struct_initial'],
                ccd_data['struct_final'],
                ccd_data['energy_initial'],
                ccd_data['energy_final'],
                n_images=15,
                align=False,  # Already aligned
                verbose=verbose
            )
            Q_f = Q_i  # Same Q grid

        # Create Potential objects
        pot_initial = Potential(
            Q_data=Q_i, E_data=E_i,
            name=f"{defect_name}_q{charge_initial:+d}"
        )
        pot_final = Potential(
            Q_data=Q_f, E_data=E_f,
            name=f"{defect_name}_q{charge_final:+d}"
        )

        # Step 4: Fit potentials
        if verbose:
            print(f"\nFitting potentials with {fit_type}...")

        pot_initial.fit(fit_type=fit_type, **fit_kwargs)
        pot_final.fit(fit_type=fit_type, **fit_kwargs)

        actual_fit_type = fit_type

    # Step 5: Solve Schrodinger equation
    if verbose:
        print(f"\nSolving Schrodinger equation (nev_i={nev_initial}, nev_f={nev_final})...")

    pot_initial.solve(nev=nev_initial)
    pot_final.solve(nev=nev_final)

    if verbose:
        print(f"Initial potential: {len(pot_initial.eigenvalues)} eigenvalues")
        print(f"Final potential: {len(pot_final.eigenvalues)} eigenvalues")

    # Step 6: Calculate Q0
    if verbose:
        print(f"\nCalculating Q0 using method='{Q0_method}'...")

    if Q0_method == "auto":
        Q0_method = "crossing"

    Q0_data = calculate_Q0_crossing(pot_initial, pot_final, method=Q0_method)
    Q0 = Q0_data['Q0']

    if verbose:
        print(f"Q0 = {Q0:.4f} amu^0.5*A")
        if Q0_data['barrier_initial'] is not None:
            print(f"Barrier (initial): {Q0_data['barrier_initial']:.4f} eV")

    # Step 7: Create ConfigCoordinate
    if verbose:
        print("\nCreating ConfigCoordinate...")

    cc = ConfigCoordinate(
        pot_i=pot_initial,
        pot_f=pot_final,
        W=W,
        degeneracy=degeneracy,
    )

    # Prepare metadata
    metadata = {
        'ccd_data': ccd_data,
        'Q0': Q0,
        'Q0_method': Q0_data['method'],
        'Q0_data': Q0_data,
        'pot_initial': pot_initial,
        'pot_final': pot_final,
        'fit_type': actual_fit_type,
        'hw_estimated': hw_estimated,
    }

    if verbose:
        print("\nConfigCoordinate created successfully!")
        print(f"Next steps:")
        print(f"  cc.calculate_overlap(Q0={Q0:.2f}, cutoff=0.25, sigma=0.01)")
        print(f"  cc.calculate_capture_coefficient(volume=..., temperature=...)")

    return cc, metadata


__all__ = [
    "load_defect_entry",
    "get_available_charge_states",
    "validate_charge_states",
    "suggest_Q0",
    "load_path_calculations",
    "extract_cc_data_from_structures",
    "create_potential_from_doped",
    # Workflow helper functions
    "prepare_ccd_structures",
    "generate_ccd_path",
    "estimate_phonon_frequency",
    "calculate_Q0_crossing",
    "create_ccd_from_defect_entries",
]
