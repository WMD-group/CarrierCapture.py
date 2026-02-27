# Quick Start

Get started with CarrierCapture in 5 minutes!

## Installation

If you haven't installed CarrierCapture yet:

```bash
pip install carriercapture  # Coming soon
# or from source
git clone https://github.com/WMD-group/CarrierCapture.py.git
cd CarrierCapture.py
pip install -e .
```

## Your First Calculation

### Python API

Calculate the capture coefficient for a simple harmonic oscillator model:

```python
import numpy as np
from carriercapture.core.potential import Potential
from carriercapture.core.config_coord import ConfigCoordinate

# Create two displaced harmonic potentials
pot_initial = Potential.from_harmonic(
    hw=0.008,  # Phonon energy (eV)
    Q0=0.0,    # Equilibrium position (amu^0.5·Å)
    E0=0.5     # Energy offset (eV)
)

pot_final = Potential.from_harmonic(
    hw=0.008,
    Q0=10.5,   # Displaced by 10.5 amu^0.5·Å
    E0=0.0
)

# Solve Schrödinger equation for both potentials
pot_initial.solve(nev=180)  # 180 eigenvalues
pot_final.solve(nev=60)     # 60 eigenvalues

# Calculate capture coefficient
cc = ConfigCoordinate(
    pot_i=pot_initial,
    pot_f=pot_final,
    W=0.068  # Electron-phonon coupling (eV)
)

# Calculate overlap matrix
cc.calculate_overlap(Q0=5.0, cutoff=0.25, sigma=0.01)

# Calculate capture coefficient vs temperature
temperatures = np.linspace(100, 500, 50)  # 100-500 K
cc.calculate_capture_coefficient(
    volume=1e-21,      # Supercell volume (cm³)
    temperature=temperatures
)

# Results
print(f"Capture coefficient at 300K: {cc.capture_coefficient[20]:.3e} cm³/s")
print(f"Capture coefficient at 100K: {cc.capture_coefficient[0]:.3e} cm³/s")
print(f"Capture coefficient at 500K: {cc.capture_coefficient[-1]:.3e} cm³/s")
```

**Expected output:**

```
Capture coefficient at 300K: 6.824e-14 cm³/s
Capture coefficient at 100K: 3.142e-20 cm³/s
Capture coefficient at 500K: 2.135e-11 cm³/s
```

### Command-Line Interface

The same calculation using the CLI with a configuration file:

**config.yaml:**

```yaml
potential_initial:
  type: harmonic
  hw: 0.008
  Q0: 0.0
  E0: 0.5
  nev: 180

potential_final:
  type: harmonic
  hw: 0.008
  Q0: 10.5
  E0: 0.0
  nev: 60

capture:
  W: 0.068
  Q0: 5.0
  cutoff: 0.25
  sigma: 0.01
  volume: 1.0e-21
  temperature:
    min: 100
    max: 500
    points: 50
```

**Run the calculation:**

```bash
carriercapture capture config.yaml -O results.json
```

**Visualize results:**

```bash
carriercapture plot results.json --show
```

## Understanding the Results

### What We Calculated

- **Initial state**: Excited electronic state with higher energy
- **Final state**: Ground electronic state
- **Capture**: Non-radiative transition from initial → final via phonon emission

### Physical Interpretation

The capture coefficient tells us:

- **At low T (100K)**: Very slow capture (10⁻²⁰ cm³/s) - needs thermal activation
- **At room T (300K)**: Moderate capture (10⁻¹⁴ cm³/s) - typical for defects
- **At high T (500K)**: Fast capture (10⁻¹¹ cm³/s) - thermally enhanced

### Temperature Dependence

The exponential increase with temperature is characteristic of:

1. **Thermal activation** over barriers
2. **Phonon emission** processes
3. **Huang-Rhys** multiphonon theory

## What's Next?

### Learn the Concepts

- **[Basic Concepts](basic-concepts.md)** - Understand potentials and CC diagrams
- **[First Calculation](first-calculation.md)** - Detailed step-by-step tutorial

### Explore Features

- **[Parameter Scanning](../user-guide/parameter-scanning.md)** - High-throughput screening
- **[Visualization](../user-guide/visualization.md)** - Interactive dashboards
- **[doped Integration](../user-guide/doped-integration.md)** - Real defect calculations

### Try the Example Notebooks

- **[Notebook 1: Sn in ZnO](https://github.com/WMD-group/CarrierCapture.py/blob/main/examples/notebooks/01_harmonic_sn_zn.ipynb)** - Detailed harmonic example
- **[Notebook 3: Parameter Scan](https://github.com/WMD-group/CarrierCapture.py/blob/main/examples/notebooks/03_parameter_scan.ipynb)** - Materials screening

## Common Next Steps

### Visualize Your Results

```python
from carriercapture.visualization import plot_capture_coefficient, plot_potential

# Plot capture coefficient vs temperature
fig = plot_capture_coefficient(cc, title="Harmonic Oscillator")
fig.show()

# Plot potential energy surfaces
fig = plot_potential(pot_initial, show_wavefunctions=True)
fig.show()
```

### Run a Parameter Scan

```python
from carriercapture.analysis import ParameterScanner, ScanParameters

# Define scan ranges
params = ScanParameters(
    dQ_range=(0, 25, 25),   # ΔQ: 0-25, 25 points
    dE_range=(0, 2.5, 10),  # ΔE: 0-2.5 eV, 10 points
    hbar_omega_i=0.008,
    hbar_omega_f=0.008,
    temperature=300.0,
    volume=1e-21
)

# Run scan
scanner = ParameterScanner(params)
results = scanner.run_harmonic_scan(n_jobs=-1)  # Use all CPU cores

print(f"Scanned {results.capture_coefficients.size} parameter combinations")
```

### Launch Interactive Dashboard

```bash
carriercapture viz --port 8050
```

Then open http://localhost:8050 in your browser for real-time exploration!

## Troubleshooting

### Convergence Issues

If you see warnings about partition function convergence:

```python
# Increase number of eigenvalues
pot_initial.solve(nev=200)  # Was 180
pot_final.solve(nev=80)     # Was 60
```

### Slow Calculations

For faster calculations:

```python
# Use coarser Q grid
pot = Potential.from_harmonic(hw=0.008, Q0=0.0, E0=0.5, npoints=2000)  # Was 5000

# Reduce eigenvalues if convergence allows
pot.solve(nev=50)  # Was 180
```

### Memory Issues

For large parameter scans:

```python
# Use fewer grid points
params = ScanParameters(
    dQ_range=(0, 25, 15),   # Reduced from 25
    dE_range=(0, 2.5, 8),   # Reduced from 10
    ...
)

# Or run serially
results = scanner.run_harmonic_scan(n_jobs=1)
```

## Getting Help

- **Documentation**: You're reading it!
- **Examples**: Check [examples/](https://github.com/WMD-group/CarrierCapture.py/tree/main/examples)
- **Issues**: [GitHub Issues](https://github.com/WMD-group/CarrierCapture.py/issues)
- **Discussions**: [GitHub Discussions](https://github.com/WMD-group/CarrierCapture.py/discussions)
