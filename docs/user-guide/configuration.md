# Configuration Files

YAML configuration for reproducible carrier capture calculations.

## Overview

Configuration files (YAML format) provide a clean way to specify all parameters for capture coefficient calculations. They are:

- **Reproducible** - Document exact parameters used
- **Shareable** - Easy to share with collaborators
- **Version-controlled** - Track parameter changes with git
- **Reusable** - Run same calculation with different data

---

## Basic Structure

```yaml
# config.yaml - Basic capture calculation

potential_initial:
  file: excited_solved.json

potential_final:
  file: ground_solved.json

capture:
  W: 0.205              # Electron-phonon coupling (eV)
  degeneracy: 1         # Degeneracy factor
  volume: 1.0e-21       # Supercell volume (cm³)
  Q0: 10.0              # Coordinate shift (amu^0.5·Å)
  cutoff: 0.25          # Energy cutoff (eV)
  sigma: 0.025          # Gaussian width (eV)
  temperature:
    min: 100            # Minimum temperature (K)
    max: 500            # Maximum temperature (K)
    n_points: 50        # Number of points
```

**Usage:**
```bash
carriercapture capture config.yaml -o results.json -v
```

---

## Complete Example

```yaml
# Complete configuration with all options

# Metadata (optional, for documentation)
metadata:
  name: "Sn_Zn Carrier Capture"
  material: "ZnO"
  defect: "Sn on Zn site"
  transition: "q=0 → q=+1"
  date: "2024-01-15"
  author: "J. Doe"
  notes: "DFT with PBE+U, U_Zn=10eV, 3×3×2 supercell"

# Initial state potential
potential_initial:
  file: Sn_Zn_excited_solved.json
  # Alternative: specify data directly
  # Q_data: [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20]
  # E_data: [0.0, 0.05, 0.18, 0.35, 0.48, 0.55, 0.58, 0.59, 0.60, 0.61, 0.62]
  # fit_type: spline
  # fit_params:
  #   order: 4
  #   smoothness: 0.001
  # nev: 180  # Solve if not already solved

# Final state potential
potential_final:
  file: Sn_Zn_ground_solved.json

# Capture calculation parameters
capture:
  # Coupling parameters
  W: 0.205              # Electron-phonon coupling (eV)
  degeneracy: 1         # Degeneracy factor (g)

  # System parameters
  volume: 1.08e-21      # Supercell volume (cm³)
  Q0: 10.0              # Configuration coordinate shift (amu^0.5·Å)

  # Numerical parameters
  cutoff: 0.25          # Energy cutoff for overlaps (eV)
  sigma: 0.025          # Gaussian delta width (eV)

  # Temperature range
  temperature:
    min: 100            # Minimum temperature (K)
    max: 500            # Maximum temperature (K)
    n_points: 50        # Number of temperature points

  # Optional: specific temperatures
  # temperature: [100, 200, 300, 400, 500]  # List of temperatures

# Output options (optional)
output:
  format: json          # Output format: json, yaml, csv, npz
  save_overlaps: true   # Save overlap matrix
  save_potentials: true # Include potentials in output
  precision: 6          # Number of significant figures
```

---

## Section Reference

### metadata (Optional)

Document calculation details:

```yaml
metadata:
  name: "String - calculation name"
  material: "String - host material"
  defect: "String - defect description"
  transition: "String - charge state transition"
  date: "String - YYYY-MM-DD"
  author: "String - your name"
  notes: "String - any additional notes"
  tags: ["tag1", "tag2"]  # List of tags
```

**Not used by calculation** - purely for documentation.

### potential_initial

Specify initial state (before carrier capture):

```yaml
potential_initial:
  # Option 1: Load from file (recommended)
  file: excited_solved.json

  # Option 2: Inline data
  Q_data: [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20]
  E_data: [0.0, 0.05, 0.18, 0.35, 0.48, 0.55, 0.58, 0.59, 0.60, 0.61, 0.62]

  # If providing inline data, also specify fitting
  fit_type: spline  # spline, harmonic, morse, polynomial
  fit_params:
    order: 4          # Spline order
    smoothness: 0.001 # Smoothing parameter

  # Solving parameters (if not already solved)
  nev: 180            # Number of eigenvalues
  Q_range: [-20, 20]  # Q grid range (optional)
  npoints: 3001       # Grid points (optional)
```

### potential_final

Specify final state (after carrier capture):

```yaml
potential_final:
  file: ground_solved.json

  # Or inline data (same options as potential_initial)
  # Q_data: [...]
  # E_data: [...]
  # fit_type: spline
  # nev: 60
```

### capture (Required)

Core capture calculation parameters:

```yaml
capture:
  # ---- Required parameters ----
  W: 0.205              # Electron-phonon coupling (eV)
                        # Typical range: 0.1 - 0.5 eV

  volume: 1.0e-21       # Supercell volume (cm³)
                        # From DFT: V_ang3 * 1e-24

  Q0: 10.0              # Configuration coordinate shift (amu^0.5·Å)
                        # Estimate: ΔQ/2 or scan to optimize

  # ---- Optional parameters ----
  degeneracy: 1         # Degeneracy factor g (default: 1)
                        # Examples: 1 (non-degenerate)
                        #           2 (doublet)
                        #           3 (triplet)

  cutoff: 0.25          # Energy cutoff (eV, default: 0.25)
                        # Only calculate overlaps for |E_i - E_f| < cutoff

  sigma: 0.025          # Gaussian delta width (eV, default: 0.025)
                        # Broadening of energy-conserving delta

  # ---- Temperature options ----

  # Option 1: Range (most common)
  temperature:
    min: 100            # K
    max: 500            # K
    n_points: 50        # Number of points

  # Option 2: Explicit list
  # temperature: [100, 200, 300, 400, 500]

  # Option 3: Single temperature
  # temperature: 300
```

### output (Optional)

Control output format:

```yaml
output:
  format: json          # json, yaml, csv, npz (default: json)
  save_overlaps: true   # Include overlap matrix (default: false)
  save_potentials: true # Include full potentials (default: false)
  precision: 6          # Significant figures (default: 6)
  compress: true        # Compress output (npz only, default: true)
```

---

## Parameter Guidelines

### Choosing W (Coupling Strength)

```yaml
capture:
  W: 0.205  # eV

  # Guidelines:
  # - Typical range: 0.1 - 0.5 eV
  # - From theory: ⟨ψ_e|∂V/∂Q|ψ_h⟩
  # - If unknown, use 0.2 eV as default
  # - Can fit to experimental data
  # - Scales capture coefficient quadratically: C ∝ W²
```

### Choosing volume (Supercell Volume)

```yaml
capture:
  volume: 1.08e-21  # cm³

  # From DFT supercell:
  # 1. Get volume in Ų: V_ang3 = structure.volume
  # 2. Convert: V_cm3 = V_ang3 * 1e-24
  #
  # Examples:
  # - 2×2×2 cubic (a=5.2Å): V ~ 1.1e-21 cm³
  # - 3×3×3 cubic: V ~ 3.8e-21 cm³
  # - 4×4×4 cubic: V ~ 8.9e-21 cm³
```

### Choosing Q0 (Coordinate Shift)

```yaml
capture:
  Q0: 10.0  # amu^0.5·Å

  # Guidelines:
  # 1. From structures:
  #    Q0 ≈ sqrt(M) * Δr
  #    where M is reduced mass, Δr is displacement
  #
  # 2. From potential minima:
  #    Q0 ≈ |Q_min_initial - Q_min_final|
  #
  # 3. Rule of thumb:
  #    Q0 ≈ (ΔQ_total) / 2
  #
  # 4. Optimize:
  #    Scan Q0 values: [8, 10, 12, 14]
  #    Choose Q0 that maximizes capture
```

### Choosing cutoff and sigma

```yaml
capture:
  cutoff: 0.25   # eV
  sigma: 0.025   # eV

  # cutoff: Energy window for overlap calculation
  # - Standard: 0.25 eV (~10 kT at 300K)
  # - Strict: 0.1 eV
  # - Broad: 0.5 eV
  # - Affects: Computational cost, physically reasonable

  # sigma: Gaussian broadening of delta function
  # - Standard: 0.025 eV (~1 kT at 300K)
  # - Strict: 0.01 eV
  # - Broad: 0.05 eV
  # - Rule of thumb: sigma ≈ (phonon spacing) / 2
```

### Choosing Temperature Range

```yaml
capture:
  temperature:
    min: 100
    max: 500
    n_points: 50

  # Guidelines:
  # - min: 100K (below this, quantum effects important)
  # - max: 500-1000K (depending on material)
  # - n_points: 50 for smooth Arrhenius plot
  #
  # For specific comparison with experiment:
  # temperature: [150, 200, 250, 300, 350, 400]
```

---

## Multiple Configurations

### Example 1: Parameter Scan Config

```yaml
# scan_config.yaml - For parameter sensitivity analysis

capture_base:
  W: 0.205
  volume: 1.0e-21
  cutoff: 0.25
  sigma: 0.025
  temperature:
    min: 100
    max: 500
    n_points: 50

# Scan over Q0 values
Q0_values: [8, 10, 12, 14, 16]

# Scan over W values
W_values: [0.15, 0.20, 0.25, 0.30]
```

**Usage script:**
```python
import yaml
from carriercapture.core import ConfigCoordinate, Potential
from carriercapture.io import load_potential_from_file
import numpy as np

# Load config
with open('scan_config.yaml') as f:
    config = yaml.safe_load(f)

# Load potentials
pot_i = Potential.from_dict(load_potential_from_file('excited.json'))
pot_f = Potential.from_dict(load_potential_from_file('ground.json'))

# Scan Q0
results = {}
for Q0 in config['Q0_values']:
    cc = ConfigCoordinate(pot_i, pot_f, W=config['capture_base']['W'], degeneracy=1)
    cc.calculate_overlap(Q0=Q0, cutoff=0.25, sigma=0.025)
    cc.calculate_capture_coefficient(volume=1e-21, temperature=np.array([300]))
    results[Q0] = cc.capture_coefficient[0]

print("Q0 scan results:")
for Q0, C in results.items():
    print(f"  Q0={Q0:4.1f}: C(300K)={C:.3e} cm³/s")
```

### Example 2: Multiple Defects

```yaml
# defects_config.yaml - Configuration for multiple defects

defaults:
  capture:
    W: 0.205
    degeneracy: 1
    volume: 1.08e-21
    cutoff: 0.25
    sigma: 0.025
    temperature:
      min: 100
      max: 500
      n_points: 50

defects:
  Sn_Zn:
    potential_initial: data/Sn_Zn_excited.json
    potential_final: data/Sn_Zn_ground.json
    Q0: 10.0
    notes: "Sn substituting Zn in ZnO"

  V_Zn:
    potential_initial: data/V_Zn_excited.json
    potential_final: data/V_Zn_ground.json
    Q0: 8.5
    W: 0.18  # Override default
    notes: "Zn vacancy"

  Sn_O:
    potential_initial: data/Sn_O_excited.json
    potential_final: data/Sn_O_ground.json
    Q0: 12.0
    notes: "Sn substituting O in ZnO"
```

---

## Templates

### Template 1: Standard Calculation

```yaml
# template_standard.yaml

metadata:
  name: "DEFECT_NAME"
  material: "MATERIAL"
  date: "YYYY-MM-DD"

potential_initial:
  file: excited_solved.json

potential_final:
  file: ground_solved.json

capture:
  W: 0.205
  volume: 1.0e-21
  Q0: 10.0
  degeneracy: 1
  cutoff: 0.25
  sigma: 0.025
  temperature:
    min: 100
    max: 500
    n_points: 50
```

### Template 2: With Inline Data

```yaml
# template_inline.yaml

metadata:
  name: "Test Calculation"
  notes: "Inline Q-E data for testing"

potential_initial:
  Q_data: [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20]
  E_data: [0.0, 0.05, 0.18, 0.35, 0.48, 0.55, 0.58, 0.59, 0.60, 0.61, 0.62]
  fit_type: spline
  fit_params:
    order: 4
    smoothness: 0.001
  nev: 180

potential_final:
  Q_data: [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20]
  E_data: [1.5, 1.48, 1.42, 1.32, 1.18, 1.02, 0.88, 0.78, 0.72, 0.70, 0.70]
  fit_type: spline
  fit_params:
    order: 4
    smoothness: 0.001
  nev: 60

capture:
  W: 0.205
  volume: 1.0e-21
  Q0: 10.0
  temperature:
    min: 100
    max: 500
    n_points: 50
```

### Template 3: High-Precision

```yaml
# template_high_precision.yaml

metadata:
  name: "High Precision Calculation"
  notes: "Dense grid, tight convergence"

potential_initial:
  file: excited.json
  nev: 200
  npoints: 5001
  Q_range: [-25, 25]

potential_final:
  file: ground.json
  nev: 100
  npoints: 5001
  Q_range: [-25, 25]

capture:
  W: 0.205
  volume: 1.0e-21
  Q0: 10.0
  cutoff: 0.15      # Stricter
  sigma: 0.01       # Narrower
  temperature:
    min: 50
    max: 600
    n_points: 100   # More points

output:
  format: npz
  save_overlaps: true
  save_potentials: true
  precision: 8
```

---

## Validation

### Validating Config Files

```python
from carriercapture.io import validate_config

# Check if config is valid
try:
    validate_config('config.yaml')
    print("✓ Config is valid")
except ValueError as e:
    print(f"✗ Config error: {e}")
```

### Schema

Required fields:
- `potential_initial` (dict)
  - `file` (str) OR (`Q_data` + `E_data` + `fit_type`)
- `potential_final` (dict)
  - `file` (str) OR (`Q_data` + `E_data` + `fit_type`)
- `capture` (dict)
  - `W` (float)
  - `volume` (float)
  - `Q0` (float)

Optional fields:
- `metadata` (dict)
- `capture.degeneracy` (int, default: 1)
- `capture.cutoff` (float, default: 0.25)
- `capture.sigma` (float, default: 0.025)
- `capture.temperature` (dict or list or float, default: {min: 100, max: 500, n_points: 50})
- `output` (dict)

---

## Best Practices

### 1. Version Control

```bash
# Initialize git repo
git init
git add config.yaml data/
git commit -m "Initial configuration"

# Track changes
git add config.yaml
git commit -m "Update Q0 from 10 to 12 based on structure analysis"
```

### 2. Documentation

```yaml
metadata:
  name: "Sn_Zn Capture Calculation v3"
  notes: |
    Version 3 updates:
    - Increased nev_initial from 150 to 180 for better convergence
    - Adjusted Q0 from 10.0 to 10.5 based on structure displacement
    - Changed W from 0.20 to 0.205 to match Ref. [1]

    References:
    [1] Smith et al., Phys. Rev. B 95, 144102 (2023)

  changelog:
    - date: "2024-01-10"
      changes: "Initial version"
    - date: "2024-01-12"
      changes: "Updated Q0 and W parameters"
```

### 3. Reproducibility

```yaml
# Include all details needed to reproduce
metadata:
  dft_code: "VASP 6.3.0"
  dft_functional: "PBE+U"
  dft_parameters:
    U_Zn: 10.0
    ENCUT: 500
    KPOINTS: "4×4×3 Γ-centered"
  supercell: "3×3×2 of ZnO wurtzite"
  relaxation_threshold: 0.01  # eV/Å

potential_initial:
  file: excited_solved.json
  source: "VASP path calculation, 11 images"
  fit_quality:
    RMSE: 0.003  # eV
    method: "spline, order=4, smoothness=0.001"
```

---

## See Also

- **[CLI Reference: capture](../api/cli.md#capture---calculate-capture-coefficient)** - Command-line usage
- **[User Guide: CLI Usage](cli-usage.md)** - Practical CLI examples
- **[Examples](../examples/notebooks.md)** - Example configurations
