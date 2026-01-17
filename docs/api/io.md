# I/O API

File input/output and integration with external tools.

## Readers

Functions for loading potential data from various file formats.

::: carriercapture.io.readers
    options:
      show_root_heading: true
      show_source: false
      heading_level: 3

## Writers

Functions for exporting results to various formats.

::: carriercapture.io.writers
    options:
      show_root_heading: true
      show_source: false
      heading_level: 3

## doped Interface

Integration with the [doped](https://github.com/SMTG-Bham/doped) package for defect calculations.

!!! info "Optional Dependency"
    The doped interface requires `pip install carriercapture[doped]`

::: carriercapture.io.doped_interface
    options:
      show_root_heading: true
      show_source: false
      heading_level: 3

## Usage Examples

### Loading Data

```python
from carriercapture.io import load_potential, read_csv_data
from carriercapture.core import Potential

# Load from CSV file
Q_data, E_data = read_csv_data('potential_data.csv')

# Create potential from data
pot = Potential(Q_data=Q_data, E_data=E_data)
pot.fit(fit_type='spline', order=4, smoothness=0.001)
```

### Saving Results

```python
from carriercapture.io import save_results
from carriercapture.core import ConfigCoordinate

# After calculating capture coefficient
cc = ConfigCoordinate(...)
cc.calculate_capture_coefficient(...)

# Save to JSON
save_results(cc, 'results.json', format='json')

# Save to HDF5
save_results(cc, 'results.h5', format='hdf5')
```

### doped Integration

```python
from carriercapture.io.doped_interface import (
    load_defect_entry,
    create_potential_from_doped
)

# Load defect data from doped
defect = load_defect_entry('defect.json.gz')

# Create potentials for charge state transition
pot_initial = create_potential_from_doped(defect, charge_state=0)
pot_final = create_potential_from_doped(defect, charge_state=+1)

# Continue with standard workflow
pot_initial.solve(nev=180)
pot_final.solve(nev=60)
# ...
```

## File Formats

### Supported Input Formats

| Format | Extension | Description |
|--------|-----------|-------------|
| CSV | `.csv`, `.dat` | Comma or space-separated values |
| JSON | `.json` | Structured JSON with metadata |
| NPZ | `.npz` | NumPy compressed arrays |
| HDF5 | `.h5`, `.hdf5` | Hierarchical data format |
| doped | `.json.gz` | doped DefectEntry files |

### Supported Output Formats

| Format | Extension | Use Case |
|--------|-----------|----------|
| JSON | `.json` | Human-readable, portable |
| NPZ | `.npz` | Fast, compact, NumPy-native |
| HDF5 | `.h5` | Large datasets, hierarchical |
| CSV | `.csv` | Spreadsheet software |

## See Also

- **[User Guide: doped Integration](../user-guide/doped-integration.md)** - Complete doped workflow
- **[Getting Started: First Calculation](../getting-started/first-calculation.md)** - File I/O examples
