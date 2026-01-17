"""
File writers for saving results.

Supports exporting to:
- JSON (human-readable, full precision)
- YAML (configuration files)
- CSV (tabular data)
- NumPy NPZ (compressed arrays)
- Plain text (simple Q-E data)
"""

from pathlib import Path
from typing import Union, Dict, Any
import numpy as np
from numpy.typing import NDArray


def write_json(
    data: Dict[str, Any],
    filepath: Union[str, Path],
    indent: int = 2,
) -> None:
    """
    Write data to JSON file.

    Parameters
    ----------
    data : dict
        Data to write
    filepath : str or Path
        Output file path
    indent : int, default=2
        Indentation level for pretty printing

    Examples
    --------
    >>> pot_data = pot.to_dict()
    >>> write_json(pot_data, "potential.json")
    """
    import json

    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)

    with open(filepath, 'w') as f:
        json.dump(data, f, indent=indent)


def write_yaml(
    data: Dict[str, Any],
    filepath: Union[str, Path],
) -> None:
    """
    Write data to YAML file.

    Parameters
    ----------
    data : dict
        Data to write
    filepath : str or Path
        Output file path

    Examples
    --------
    >>> config = {'potential': {...}, 'capture': {...}}
    >>> write_yaml(config, "config.yaml")
    """
    import yaml

    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)

    with open(filepath, 'w') as f:
        yaml.dump(data, f, default_flow_style=False, sort_keys=False)


def write_csv(
    Q_data: NDArray[np.float64],
    E_data: NDArray[np.float64],
    filepath: Union[str, Path],
    header: bool = True,
) -> None:
    """
    Write Q-E data to CSV file.

    Parameters
    ----------
    Q_data : NDArray[np.float64]
        Configuration coordinates
    E_data : NDArray[np.float64]
        Potential energies
    filepath : str or Path
        Output file path
    header : bool, default=True
        Include header row

    Examples
    --------
    >>> write_csv(Q, E, "potential.csv")
    """
    import pandas as pd

    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)

    df = pd.DataFrame({'Q': Q_data, 'E': E_data})
    df.to_csv(filepath, index=False, header=header)


def write_npz(
    data: Dict[str, NDArray],
    filepath: Union[str, Path],
    compressed: bool = True,
) -> None:
    """
    Write arrays to NumPy NPZ file.

    Parameters
    ----------
    data : dict
        Dictionary of arrays to save
    filepath : str or Path
        Output file path
    compressed : bool, default=True
        Use compression

    Examples
    --------
    >>> data = {'Q': Q_data, 'E': E_data, 'eigenvalues': eigs}
    >>> write_npz(data, "results.npz")
    """
    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)

    if compressed:
        np.savez_compressed(filepath, **data)
    else:
        np.savez(filepath, **data)


def write_potential_data(
    Q_data: NDArray[np.float64],
    E_data: NDArray[np.float64],
    filepath: Union[str, Path],
    header: str = None,
    fmt: str = "%.10e",
) -> None:
    """
    Write Q-E data to plain text file.

    Parameters
    ----------
    Q_data : NDArray[np.float64]
        Configuration coordinates
    E_data : NDArray[np.float64]
        Potential energies
    filepath : str or Path
        Output file path
    header : str, optional
        Header comment line
    fmt : str, default="%.10e"
        Format string for numbers

    Examples
    --------
    >>> write_potential_data(Q, E, "potential.dat", header="Q(amu^0.5·Å)  E(eV)")
    """
    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)

    data = np.column_stack([Q_data, E_data])

    if header is not None:
        header = f"# {header}"

    np.savetxt(filepath, data, fmt=fmt, header=header or "")


def write_capture_results(
    config_coord,
    filepath: Union[str, Path],
    file_format: str = "json",
    include_partial: bool = False,
) -> None:
    """
    Write capture coefficient results to file.

    Parameters
    ----------
    config_coord : ConfigCoordinate
        Configuration coordinate with computed results
    filepath : str or Path
        Output file path
    file_format : str, default="json"
        Output format: json, yaml, csv, npz
    include_partial : bool, default=False
        Include partial capture coefficients (large arrays)

    Examples
    --------
    >>> write_capture_results(cc, "results.json")
    >>> write_capture_results(cc, "results.csv", file_format="csv")
    """
    filepath = Path(filepath)

    if file_format == "json":
        data = config_coord.to_dict()
        if not include_partial and "partial_capture_coefficient" in data:
            del data["partial_capture_coefficient"]
        write_json(data, filepath)

    elif file_format == "yaml":
        data = config_coord.to_dict()
        if not include_partial and "partial_capture_coefficient" in data:
            del data["partial_capture_coefficient"]
        write_yaml(data, filepath)

    elif file_format == "csv":
        # Write temperature and capture coefficient to CSV
        if config_coord.temperature is None or config_coord.capture_coefficient is None:
            raise ValueError("No capture coefficient data to write")

        write_csv(
            config_coord.temperature,
            config_coord.capture_coefficient,
            filepath,
            header=True
        )

    elif file_format == "npz":
        data = {}
        if config_coord.temperature is not None:
            data['temperature'] = config_coord.temperature
        if config_coord.capture_coefficient is not None:
            data['capture_coefficient'] = config_coord.capture_coefficient
        if include_partial and config_coord.partial_capture_coefficient is not None:
            data['partial_capture_coefficient'] = config_coord.partial_capture_coefficient

        write_npz(data, filepath)

    else:
        raise ValueError(f"Unsupported format: {file_format}")


def save_potential(
    potential,
    filepath: Union[str, Path],
    file_format: str = None,
) -> None:
    """
    Save Potential object to file (auto-detect format).

    Parameters
    ----------
    potential : Potential
        Potential object to save
    filepath : str or Path
        Output file path
    file_format : str, optional
        Force specific format (json, yaml, npz, dat, csv)
        If None, auto-detect from extension

    Examples
    --------
    >>> save_potential(pot, "potential.json")
    >>> save_potential(pot, "potential.npz")
    """
    filepath = Path(filepath)

    # Auto-detect format from extension
    if file_format is None:
        ext = filepath.suffix.lower()
        format_map = {
            '.json': 'json',
            '.yaml': 'yaml',
            '.yml': 'yaml',
            '.npz': 'npz',
            '.dat': 'dat',
            '.txt': 'dat',
            '.csv': 'csv',
        }
        file_format = format_map.get(ext, 'json')

    # Save based on format
    if file_format == 'json':
        write_json(potential.to_dict(), filepath)
    elif file_format == 'yaml':
        write_yaml(potential.to_dict(), filepath)
    elif file_format == 'npz':
        data = {}
        if potential.Q is not None:
            data['Q'] = potential.Q
        if potential.E is not None:
            data['E'] = potential.E
        if potential.eigenvalues is not None:
            data['eigenvalues'] = potential.eigenvalues
        if potential.eigenvectors is not None:
            data['eigenvectors'] = potential.eigenvectors
        write_npz(data, filepath)
    elif file_format == 'dat':
        if potential.Q_data is None or potential.E_data is None:
            raise ValueError("No Q_data/E_data to write. Use json/npz for full potential.")
        write_potential_data(
            potential.Q_data,
            potential.E_data,
            filepath,
            header="Q(amu^0.5·Å)  E(eV)"
        )
    elif file_format == 'csv':
        if potential.Q_data is None or potential.E_data is None:
            raise ValueError("No Q_data/E_data to write. Use json/npz for full potential.")
        write_csv(potential.Q_data, potential.E_data, filepath)
    else:
        raise ValueError(f"Unsupported format: {file_format}")
