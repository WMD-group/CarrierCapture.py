# API Reference

Complete API documentation for CarrierCapture.py.

## Overview

CarrierCapture is organized into several modules:

- **[Core](core.md)** - Fundamental classes for potentials, capture calculations, and Schrödinger solver
- **[Analysis](analysis.md)** - Parameter scanning for high-throughput screening
- **[I/O](io.md)** - File I/O, readers, writers, and doped integration
- **[Visualization](visualization.md)** - Plotting functions and interactive dashboards
- **[CLI](cli.md)** - Command-line interface

## Quick Links

### Most Used Classes

| Class | Module | Description |
|-------|--------|-------------|
| [`Potential`](core.md#potential) | `core.potential` | Potential energy surface |
| [`ConfigCoordinate`](core.md#configcoordinate) | `core.config_coord` | Capture coefficient calculation |
| [`ParameterScanner`](analysis.md#parameterscanner) | `analysis.parameter_scan` | Parameter sweeping |

### Quick Import Examples

```python
# Core functionality
from carriercapture.core import Potential, ConfigCoordinate

# Parameter scanning
from carriercapture.analysis import ParameterScanner, ScanParameters

# Visualization
from carriercapture.visualization import plot_potential, plot_capture_coefficient

# I/O
from carriercapture.io import load_potential, save_results

# doped integration
from carriercapture.io.doped_interface import load_defect_entry, create_potential_from_doped
```

## Package Structure

```
carriercapture/
├── core/              # Core calculation engine
│   ├── potential      # Potential class
│   ├── config_coord   # ConfigCoordinate class
│   ├── transfer_coord # TransferCoordinate class (experimental)
│   └── schrodinger    # Schrödinger solver
├── analysis/          # High-throughput tools
│   └── parameter_scan # ParameterScanner
├── io/                # Input/output
│   ├── readers        # File readers
│   ├── writers        # File writers
│   └── doped_interface # doped integration
├── visualization/     # Plotting
│   ├── static         # Static Plotly plots
│   ├── interactive    # Dash dashboard
│   └── themes         # Styling
└── cli/               # Command-line interface
    └── commands/      # CLI commands
```

## Documentation Style

All API documentation follows the [NumPy docstring convention](https://numpydoc.readthedocs.io/en/latest/format.html):

- **Parameters**: All function/method parameters with types
- **Returns**: Return values with types
- **Examples**: Usage examples
- **Notes**: Implementation details
- **References**: Literature citations

## Type Hints

CarrierCapture uses type hints throughout:

```python
from typing import Optional, Tuple
import numpy as np
from numpy.typing import NDArray

def solve(self, nev: int, maxiter: Optional[int] = None) -> Tuple[NDArray, NDArray]:
    """Solve Schrödinger equation."""
    ...
```

## Navigation

Browse the API documentation by module:

- **[Core API →](core.md)**
- **[Analysis API →](analysis.md)**
- **[I/O API →](io.md)**
- **[Visualization API →](visualization.md)**
- **[CLI Reference →](cli.md)**
