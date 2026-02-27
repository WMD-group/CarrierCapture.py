# Analysis API

High-throughput parameter scanning for materials screening.

## ParameterScanner

Main class for running parameter sweeps over ΔQ, ΔE, and phonon frequencies.

::: carriercapture.analysis.parameter_scan.ParameterScanner
    options:
      show_root_heading: true
      show_source: false
      members_order: source
      heading_level: 3

## ScanParameters

Configuration dataclass for defining parameter scan ranges.

::: carriercapture.analysis.parameter_scan.ScanParameters
    options:
      show_root_heading: true
      show_source: false
      members_order: source
      heading_level: 3

## ScanResult

Container for scan results with save/load functionality.

::: carriercapture.analysis.parameter_scan.ScanResult
    options:
      show_root_heading: true
      show_source: false
      members_order: source
      heading_level: 3

## Usage Example

```python
from carriercapture.analysis import ParameterScanner, ScanParameters
import numpy as np

# Define scan parameters
params = ScanParameters(
    dQ_range=(0, 25, 25),    # ΔQ: 0-25 amu^0.5·Å, 25 points
    dE_range=(0, 2.5, 10),   # ΔE: 0-2.5 eV, 10 points
    hbar_omega_i=0.008,      # Initial phonon energy
    hbar_omega_f=0.008,      # Final phonon energy
    temperature=300.0,       # Temperature (K)
    volume=1e-21,            # Supercell volume (cm³)
    degeneracy=1,
    nev_initial=180,
    nev_final=60,
)

# Create scanner
scanner = ParameterScanner(params, verbose=True)

# Run scan (use all CPU cores)
results = scanner.run_harmonic_scan(n_jobs=-1, show_progress=True)

# Save results
results.save('scan_results.npz', format='npz')

# Analyze results
print(f"Scanned {results.capture_coefficients.size} parameter combinations")
print(f"Max capture: {np.nanmax(results.capture_coefficients):.3e} cm³/s")
print(f"Min capture: {np.nanmin(results.capture_coefficients):.3e} cm³/s")
```

## See Also

- **[User Guide: Parameter Scanning](../user-guide/parameter-scanning.md)** - Detailed usage guide
- **[Example Notebook: Parameter Scan](https://github.com/WMD-group/CarrierCapture.py/blob/main/examples/notebooks/03_parameter_scan.ipynb)** - Step-by-step example
