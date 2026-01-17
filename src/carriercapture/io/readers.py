"""
File readers for various data formats.

Supports loading potential data, configuration files, and results from:
- Plain text (DAT, TXT)
- CSV
- JSON
- YAML
- NumPy NPZ
"""

from pathlib import Path
from typing import Union, Tuple, Dict, Any
import numpy as np
from numpy.typing import NDArray


def read_potential_data(
    filepath: Union[str, Path],
    delimiter: str = None,
    skip_header: int = 0,
) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    """
    Read potential energy surface data from text file.

    Expects two-column format: Q (amu^0.5·Å), E (eV)

    Parameters
    ----------
    filepath : str or Path
        Path to data file
    delimiter : str, optional
        Column delimiter (auto-detected if None)
    skip_header : int, default=0
        Number of header lines to skip

    Returns
    -------
    Q_data : NDArray[np.float64]
        Configuration coordinates (amu^0.5·Å)
    E_data : NDArray[np.float64]
        Potential energies (eV)

    Raises
    ------
    ValueError
        If file doesn't have exactly 2 columns

    Examples
    --------
    >>> Q, E = read_potential_data("excited.dat")
    >>> Q.shape, E.shape
    ((100,), (100,))
    """
    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(f"File not found: {filepath}")

    # Read data
    data = np.loadtxt(filepath, delimiter=delimiter, skiprows=skip_header)

    if data.ndim != 2 or data.shape[1] != 2:
        raise ValueError(
            f"Expected 2 columns (Q, E), got shape {data.shape}. "
            "File should have format: Q(amu^0.5·Å)  E(eV)"
        )

    Q_data = data[:, 0]
    E_data = data[:, 1]

    return Q_data, E_data


def read_json(filepath: Union[str, Path]) -> Dict[str, Any]:
    """
    Read JSON file.

    Parameters
    ----------
    filepath : str or Path
        Path to JSON file

    Returns
    -------
    data : dict
        Loaded JSON data

    Examples
    --------
    >>> data = read_json("potential.json")
    >>> pot = Potential.from_dict(data)
    """
    import json

    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(f"File not found: {filepath}")

    with open(filepath, 'r') as f:
        data = json.load(f)

    return data


def read_yaml(filepath: Union[str, Path]) -> Dict[str, Any]:
    """
    Read YAML configuration file.

    Parameters
    ----------
    filepath : str or Path
        Path to YAML file

    Returns
    -------
    config : dict
        Loaded YAML configuration

    Examples
    --------
    >>> config = read_yaml("config.yaml")
    >>> pot_config = config['potential_initial']
    """
    import yaml

    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(f"File not found: {filepath}")

    with open(filepath, 'r') as f:
        config = yaml.safe_load(f)

    return config


def read_csv(
    filepath: Union[str, Path],
    has_header: bool = True,
) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    """
    Read CSV file with potential data.

    Parameters
    ----------
    filepath : str or Path
        Path to CSV file
    has_header : bool, default=True
        Whether file has header row

    Returns
    -------
    Q_data : NDArray[np.float64]
        Configuration coordinates
    E_data : NDArray[np.float64]
        Potential energies

    Examples
    --------
    >>> Q, E = read_csv("data.csv")
    """
    import pandas as pd

    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(f"File not found: {filepath}")

    if has_header:
        df = pd.read_csv(filepath)
    else:
        df = pd.read_csv(filepath, header=None, names=["Q", "E"])

    if df.shape[1] < 2:
        raise ValueError(f"Expected at least 2 columns, got {df.shape[1]}")

    Q_data = df.iloc[:, 0].values
    E_data = df.iloc[:, 1].values

    return Q_data, E_data


def read_npz(filepath: Union[str, Path]) -> Dict[str, Any]:
    """
    Read NumPy NPZ file.

    Parameters
    ----------
    filepath : str or Path
        Path to NPZ file

    Returns
    -------
    data : dict
        Dictionary of arrays from NPZ file

    Examples
    --------
    >>> data = read_npz("results.npz")
    >>> Q = data['Q_data']
    >>> E = data['E_data']
    """
    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(f"File not found: {filepath}")

    data = np.load(filepath)
    return dict(data)


def load_potential_from_file(
    filepath: Union[str, Path],
    file_format: str = None,
) -> Dict[str, Any]:
    """
    Load potential from file (auto-detect format).

    Supports: .json, .yaml, .yml, .npz, .dat, .txt, .csv

    Parameters
    ----------
    filepath : str or Path
        Path to file
    file_format : str, optional
        Force specific format (json, yaml, npz, dat, csv)
        If None, auto-detect from extension

    Returns
    -------
    data : dict
        Potential data dictionary

    Examples
    --------
    >>> data = load_potential_from_file("potential.json")
    >>> pot = Potential.from_dict(data)
    """
    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(f"File not found: {filepath}")

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
        file_format = format_map.get(ext, 'dat')

    # Load based on format
    if file_format == 'json':
        return read_json(filepath)
    elif file_format == 'yaml':
        return read_yaml(filepath)
    elif file_format == 'npz':
        return read_npz(filepath)
    elif file_format == 'dat':
        Q_data, E_data = read_potential_data(filepath)
        return {'Q_data': Q_data, 'E_data': E_data}
    elif file_format == 'csv':
        Q_data, E_data = read_csv(filepath)
        return {'Q_data': Q_data, 'E_data': E_data}
    else:
        raise ValueError(f"Unsupported format: {file_format}")
