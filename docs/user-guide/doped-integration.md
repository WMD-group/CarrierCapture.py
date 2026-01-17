# doped Integration

Seamless integration with the [doped](https://github.com/SMTG-Bham/doped) package for automated defect calculations.

## Overview

The [doped](https://github.com/SMTG-Bham/doped) package automates DFT defect calculations and provides tools for:
- Generating defect structures
- Running VASP calculations
- Parsing and analyzing defect energetics
- Configuration coordinate diagram setup

CarrierCapture integrates directly with doped to:
1. **Load defect data** from doped DefectEntry objects
2. **Extract Q-E data** from VASP path calculations
3. **Create potentials** for different charge states
4. **Calculate capture rates** using multiphonon theory

This provides an end-to-end workflow from DFT to carrier capture rates.

---

## Installation

Install CarrierCapture with doped support:

```bash
pip install carriercapture[doped]
```

This installs:
- `doped` - Defect calculation automation
- `pymatgen` - Materials analysis
- `monty` - Serialization utilities

Verify installation:

```python
from carriercapture.io.doped_interface import load_defect_entry
print("✓ doped integration available")
```

---

## Quick Start

```python
from carriercapture.io.doped_interface import (
    load_defect_entry,
    load_path_calculations,
    create_potential_from_doped
)
from carriercapture.core import ConfigCoordinate
import numpy as np

# 1. Load defect entry
defect = load_defect_entry('Sn_Zn_defect.json.gz')

# 2. Load VASP path calculations for two charge states
Q_i, E_i = load_path_calculations('path_q0/')  # Neutral
Q_f, E_f = load_path_calculations('path_q1/')  # +1 charged

# 3. Create potentials
pot_i = create_potential_from_doped(defect, charge_state=0, Q_data=Q_i, E_data=E_i)
pot_f = create_potential_from_doped(defect, charge_state=+1, Q_data=Q_f, E_data=E_f)

# 4. Fit and solve
pot_i.fit(fit_type='spline', order=4, smoothness=0.001)
pot_f.fit(fit_type='spline', order=4, smoothness=0.001)

pot_i.solve(nev=180)
pot_f.solve(nev=60)

# 5. Calculate capture coefficient
cc = ConfigCoordinate(pot_i, pot_f, W=0.205, degeneracy=1)
cc.calculate_overlap(Q0=10.0, cutoff=0.25, sigma=0.025)

temperature = np.linspace(100, 500, 50)
cc.calculate_capture_coefficient(volume=1e-21, temperature=temperature)

print(f"C(300K) = {cc.capture_coefficient[20]:.3e} cm³/s")
```

---

## doped Package Overview

### What is doped?

[doped](https://github.com/SMTG-Bham/doped) automates defect calculations:

```python
from doped import DefectsGenerator

# Generate all defects for a structure
defect_gen = DefectsGenerator.from_structures(
    bulk_structure,
    charge_states=[-2, -1, 0, +1, +2]
)

# Write VASP input files
defect_gen.write_files()

# After VASP runs, parse results
from doped import DefectThermodynamics

thermo = DefectThermodynamics.from_directories()
thermo.get_equilibrium_concentrations(temperature=300, fermi_level=0.5)
```

### Configuration Coordinate Diagrams in doped

doped can generate CC diagram input:

```python
from doped.utils.configurations import get_cc_structures

# Generate interpolated structures between two charge states
structures = get_cc_structures(
    defect_entry,
    charge_state_initial=0,
    charge_state_final=+1,
    n_images=11  # Number of interpolation points
)

# Write VASP input for each structure
# Run VASP calculations
# Parse energies → Q-E data for CarrierCapture
```

**CarrierCapture picks up here** - loading the Q-E data and calculating capture rates.

---

## Complete Workflow

### Step 1: doped Defect Setup

```python
# In your doped workflow script

from doped import DefectsGenerator
from pymatgen.core import Structure

# Load bulk structure
bulk = Structure.from_file('POSCAR_bulk')

# Generate defects
defect_gen = DefectsGenerator.from_structures(
    bulk,
    extrinsic_elements=['Sn'],  # Dopant
    charge_states=[-1, 0, +1, +2],
    oxidation_states={'Zn': +2, 'O': -2}
)

# Write input files
defect_gen.write_files(output_path='defects')

# ... Run VASP calculations for defect energies ...
# ... Analyze with DefectThermodynamics ...
```

### Step 2: Generate CC Diagram Structures

```python
from doped.utils.configurations import get_cc_structures

# After identifying interesting defect (e.g., Sn_Zn)
structures_q0 = get_cc_structures(
    defect_entry,
    charge_state_initial=0,
    charge_state_final=0,  # Ground state path
    n_images=11
)

structures_q1 = get_cc_structures(
    defect_entry,
    charge_state_initial=+1,
    charge_state_final=+1,  # Excited state path
    n_images=11
)

# Write VASP input files for path calculations
# Directory structure:
# path_q0/
#   ├── image_00/POSCAR
#   ├── image_01/POSCAR
#   ├── ...
#   └── image_10/POSCAR
# path_q1/
#   ├── image_00/POSCAR
#   ├── ...
#   └── image_10/POSCAR
```

### Step 3: Run VASP Path Calculations

```bash
# For each image directory, run VASP
cd path_q0/image_00
vasp_std  # or your VASP command

cd path_q0/image_01
vasp_std

# ... repeat for all images in both paths ...
```

**VASP Settings:**
- Use same settings as defect calculations
- Single-point energy calculation (NSW=0)
- High accuracy (PREC=Accurate)
- Dense k-point mesh

### Step 4: CarrierCapture Analysis

```python
from carriercapture.io.doped_interface import (
    load_defect_entry,
    load_path_calculations,
    create_potential_from_doped
)
from carriercapture.core import ConfigCoordinate
from carriercapture.visualization import plot_configuration_coordinate
import numpy as np

# Load defect entry
defect = load_defect_entry('defects/Sn_Zn_0.json.gz')

# Load Q-E data from VASP calculations
print("Loading VASP path calculations...")
Q_i, E_i = load_path_calculations('path_q0/', verbose=True)
Q_f, E_f = load_path_calculations('path_q1/', verbose=True)

print(f"✓ Initial state: {len(Q_i)} images")
print(f"✓ Final state: {len(Q_f)} images")

# Create potentials
pot_i = create_potential_from_doped(
    defect,
    charge_state=0,
    Q_data=Q_i,
    E_data=E_i,
    name="Sn_Zn q=0"
)

pot_f = create_potential_from_doped(
    defect,
    charge_state=+1,
    Q_data=Q_f,
    E_data=E_f,
    name="Sn_Zn q=+1"
)

# Fit potentials
print("\nFitting potentials...")
pot_i.fit(fit_type='spline', order=4, smoothness=0.001)
pot_f.fit(fit_type='spline', order=4, smoothness=0.001)
print("✓ Fitted")

# Solve Schrödinger equation
print("\nSolving Schrödinger equation...")
pot_i.solve(nev=180)
pot_f.solve(nev=60)
print(f"✓ Initial: {len(pot_i.eigenvalues)} eigenvalues")
print(f"✓ Final: {len(pot_f.eigenvalues)} eigenvalues")

# Visualize CC diagram
fig = plot_configuration_coordinate(
    pot_i,
    pot_f,
    Q0=10.0,
    show_crossing=True,
    show_wavefunctions=True,
    title="Sn_Zn Configuration Coordinate Diagram"
)
fig.write_html('cc_diagram.html')
print("✓ CC diagram saved to cc_diagram.html")

# Calculate capture coefficient
print("\nCalculating capture coefficient...")
cc = ConfigCoordinate(pot_i, pot_f, W=0.205, degeneracy=1)
cc.calculate_overlap(Q0=10.0, cutoff=0.25, sigma=0.025)

temperature = np.linspace(100, 500, 50)
cc.calculate_capture_coefficient(volume=1e-21, temperature=temperature)

print(f"\nResults:")
print(f"  C(100K) = {cc.capture_coefficient[0]:.3e} cm³/s")
print(f"  C(300K) = {cc.capture_coefficient[24]:.3e} cm³/s")
print(f"  C(500K) = {cc.capture_coefficient[-1]:.3e} cm³/s")

# Save results
from carriercapture.io import save_potential, write_capture_results

save_potential(pot_i, 'Sn_Zn_q0_solved.json')
save_potential(pot_f, 'Sn_Zn_q1_solved.json')
write_capture_results(cc, 'Sn_Zn_capture_results.json')
print("\n✓ Results saved")
```

---

## doped Interface Functions

### Loading Defect Data

```python
from carriercapture.io.doped_interface import load_defect_entry

# Load from doped JSON.GZ file
defect = load_defect_entry('defect_entry.json.gz')

# Access defect properties
print(f"Defect: {defect.name}")
print(f"Charge: {defect.charge_state:+d}")
print(f"Structure: {defect.structure}")

# Check available charge states
from carriercapture.io.doped_interface import get_available_charge_states
charges = get_available_charge_states(defect)
print(f"Available charges: {charges}")
```

### Loading VASP Path Calculations

```python
from carriercapture.io.doped_interface import load_path_calculations

# Load from directory of VASP calculations
Q_data, E_data = load_path_calculations(
    'path_calculations/',
    image_pattern='image_*',  # Glob pattern for subdirectories
    verbose=True              # Print progress
)

print(f"Loaded {len(Q_data)} images")
print(f"Q range: {Q_data.min():.2f} to {Q_data.max():.2f} amu^0.5·Å")
print(f"E range: {E_data.min():.3f} to {E_data.max():.3f} eV")

# Q_data: Configuration coordinates (amu^0.5·Å)
# E_data: Energies relative to first image (eV)
```

**Directory structure expected:**
```
path_calculations/
├── image_00/
│   ├── vasprun.xml
│   └── OUTCAR
├── image_01/
│   ├── vasprun.xml
│   └── OUTCAR
├── ...
└── image_10/
    ├── vasprun.xml
    └── OUTCAR
```

### Creating Potentials from doped

```python
from carriercapture.io.doped_interface import create_potential_from_doped

# Create potential with Q-E data
pot = create_potential_from_doped(
    defect_entry,
    charge_state=0,
    Q_data=Q_array,    # From load_path_calculations()
    E_data=E_array,    # From load_path_calculations()
    name="Custom Name"  # Optional
)

# Potential is ready for fitting and solving
pot.fit(fit_type='spline')
pot.solve(nev=180)
```

### Suggesting Q₀ from Structures

```python
from carriercapture.io.doped_interface import suggest_Q0
from pymatgen.core import Structure

# Load relaxed structures for two charge states
struct_i = Structure.from_file('path_q0/image_00/CONTCAR')
struct_f = Structure.from_file('path_q1/image_00/CONTCAR')

# Get suggested Q0 (typically ~midpoint between structures)
Q0_suggested = suggest_Q0(struct_i, struct_f, align=True)
print(f"Suggested Q0: {Q0_suggested:.2f} amu^0.5·Å")

# Use in overlap calculation
cc.calculate_overlap(Q0=Q0_suggested, cutoff=0.25, sigma=0.025)
```

---

## Command-Line Interface

Run from command line:

```bash
# Calculate capture from doped data
carriercapture capture --doped Sn_Zn_0.json.gz \
  --charge-i 0 --charge-f +1 \
  --doped-path-i path_q0/ \
  --doped-path-f path_q1/ \
  -W 0.205 -V 1e-21 \
  --temp-range 100 500 50 \
  --auto-Q0 \
  -o capture_results.json \
  --plot --plot-output arrhenius.png \
  -vv
```

**Options:**
- `--doped` - Path to DefectEntry JSON.GZ
- `--charge-i` - Initial charge state
- `--charge-f` - Final charge state
- `--doped-path-i` - VASP path directory (initial)
- `--doped-path-f` - VASP path directory (final)
- `--auto-Q0` - Automatically suggest Q0 from structures
- `-W` - Electron-phonon coupling
- `-V` - Supercell volume
- `--plot` - Generate Arrhenius plot

---

## Advanced Usage

### Multiple Charge State Transitions

```python
# Calculate capture for multiple transitions
# e.g., 0 → +1 and +1 → +2

transitions = [
    {'charge_i': 0, 'charge_f': +1, 'path_i': 'path_q0/', 'path_f': 'path_q1/'},
    {'charge_i': +1, 'charge_f': +2, 'path_i': 'path_q1/', 'path_f': 'path_q2/'},
]

results = {}

for trans in transitions:
    print(f"\nTransition: q={trans['charge_i']:+d} → q={trans['charge_f']:+d}")

    # Load data
    Q_i, E_i = load_path_calculations(trans['path_i'])
    Q_f, E_f = load_path_calculations(trans['path_f'])

    # Create potentials
    pot_i = create_potential_from_doped(defect, trans['charge_i'], Q_i, E_i)
    pot_f = create_potential_from_doped(defect, trans['charge_f'], Q_f, E_f)

    # Fit, solve, calculate
    pot_i.fit('spline', order=4, smoothness=0.001)
    pot_f.fit('spline', order=4, smoothness=0.001)
    pot_i.solve(nev=180)
    pot_f.solve(nev=60)

    cc = ConfigCoordinate(pot_i, pot_f, W=0.205, degeneracy=1)
    cc.calculate_overlap(Q0=10.0)
    cc.calculate_capture_coefficient(volume=1e-21, temperature=np.array([300]))

    label = f"q{trans['charge_i']:+d}→q{trans['charge_f']:+d}"
    results[label] = cc.capture_coefficient[0]
    print(f"  C(300K) = {cc.capture_coefficient[0]:.3e} cm³/s")

print("\nSummary:")
for label, C in results.items():
    print(f"  {label}: {C:.3e} cm³/s")
```

### Systematic Defect Screening

```python
# Screen all defects in a material
from pathlib import Path

defect_files = list(Path('defects').glob('*_0.json.gz'))

screening_results = []

for defect_file in defect_files:
    defect_name = defect_file.stem
    print(f"\nProcessing {defect_name}...")

    try:
        # Load defect
        defect = load_defect_entry(defect_file)

        # Check if path calculations exist
        path_i = Path(f'paths/{defect_name}_q0')
        path_f = Path(f'paths/{defect_name}_q1')

        if not (path_i.exists() and path_f.exists()):
            print(f"  Skipping: path calculations not found")
            continue

        # Load and calculate
        Q_i, E_i = load_path_calculations(path_i)
        Q_f, E_f = load_path_calculations(path_f)

        pot_i = create_potential_from_doped(defect, 0, Q_i, E_i)
        pot_f = create_potential_from_doped(defect, +1, Q_f, E_f)

        pot_i.fit('spline', order=4, smoothness=0.001)
        pot_f.fit('spline', order=4, smoothness=0.001)
        pot_i.solve(nev=180)
        pot_f.solve(nev=60)

        cc = ConfigCoordinate(pot_i, pot_f, W=0.205, degeneracy=1)
        cc.calculate_overlap(Q0=10.0)
        cc.calculate_capture_coefficient(volume=1e-21, temperature=np.array([300]))

        C_300 = cc.capture_coefficient[0]
        screening_results.append({
            'defect': defect_name,
            'C_300K': C_300
        })

        print(f"  ✓ C(300K) = {C_300:.3e} cm³/s")

    except Exception as e:
        print(f"  ✗ Error: {e}")
        continue

# Sort by capture rate
screening_results.sort(key=lambda x: x['C_300K'], reverse=True)

print("\n" + "="*50)
print("Screening Results (sorted by C at 300K):")
print("="*50)
for i, result in enumerate(screening_results, 1):
    print(f"{i:2d}. {result['defect']:20s}: {result['C_300K']:.3e} cm³/s")
```

---

## Best Practices

### 1. VASP Convergence

```python
# Check that VASP calculations converged
from carriercapture.io.doped_interface import check_vasp_convergence

convergence = check_vasp_convergence('path_q0/', check_forces=True)

if not convergence['all_converged']:
    print("⚠️  Warning: Not all VASP calculations converged")
    for img, status in convergence['images'].items():
        if not status['converged']:
            print(f"  {img}: {status['reason']}")
else:
    print("✓ All VASP calculations converged")
```

### 2. Energy Alignment

```python
# Ensure potentials are aligned relative to same reference
# doped automatically handles this for DefectEntry objects

# If manually aligning:
E_i_aligned = E_i - E_i[0]  # Relative to initial image
E_f_aligned = E_f - E_f[0]  # Relative to initial image

# Then offset by charge state energy difference
dE = defect_thermodynamics.get_formation_energy(charge=+1, fermi_level=0)
E_f_aligned += dE
```

### 3. Q₀ Optimization

```python
# Scan over Q0 to find optimal value
Q0_values = np.linspace(8, 14, 7)
C_values = []

for Q0 in Q0_values:
    cc_temp = ConfigCoordinate(pot_i, pot_f, W=0.205, degeneracy=1)
    cc_temp.calculate_overlap(Q0=Q0)
    cc_temp.calculate_capture_coefficient(volume=1e-21, temperature=np.array([300]))
    C_values.append(cc_temp.capture_coefficient[0])

# Find Q0 that maximizes capture
i_max = np.argmax(C_values)
Q0_optimal = Q0_values[i_max]
print(f"Optimal Q0: {Q0_optimal:.2f} amu^0.5·Å")
print(f"C(300K) = {C_values[i_max]:.3e} cm³/s")
```

---

## Troubleshooting

### Import errors

```bash
# If you get ImportError for doped or pymatgen:
pip install --upgrade carriercapture[doped]

# Verify installation:
python -c "from carriercapture.io.doped_interface import load_defect_entry; print('OK')"
```

### VASP file not found

```python
# Check directory structure
from pathlib import Path

path_dir = Path('path_q0')
image_dirs = sorted(path_dir.glob('image_*'))

print(f"Found {len(image_dirs)} image directories:")
for img_dir in image_dirs:
    has_vasprun = (img_dir / 'vasprun.xml').exists()
    has_outcar = (img_dir / 'OUTCAR').exists()
    print(f"  {img_dir.name}: vasprun={has_vasprun}, OUTCAR={has_outcar}")
```

### Q-E data issues

```python
# Check Q-E data quality
import matplotlib.pyplot as plt

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Plot Q vs image index
ax1.plot(range(len(Q_i)), Q_i, 'o-', label='Initial')
ax1.plot(range(len(Q_f)), Q_f, 's-', label='Final')
ax1.set_xlabel('Image index')
ax1.set_ylabel('Q (amu$^{0.5}$·Å)')
ax1.legend()
ax1.set_title('Configuration Coordinate')

# Plot E vs Q
ax2.plot(Q_i, E_i, 'o-', label='Initial')
ax2.plot(Q_f, E_f, 's-', label='Final')
ax2.set_xlabel('Q (amu$^{0.5}$·Å)')
ax2.set_ylabel('E (eV)')
ax2.legend()
ax2.set_title('Potential Energy Surface')

plt.tight_layout()
plt.show()

# Check for issues:
# - Q should be monotonic
# - E should be smooth (no spikes)
# - E should have clear minimum
```

---

## Complete Example

See **[Tutorial 2: DX Center from doped](../tutorials/02-dx-center.md)** for a complete worked example.

---

## See Also

- **[doped Documentation](https://doped.readthedocs.io/)** - Full doped package documentation
- **[API Reference: doped Interface](../api/io.md#doped-interface)** - Complete API documentation
- **[CLI Reference: capture --doped](../api/cli.md#capture---calculate-capture-coefficient)** - Command-line usage
- **[Tutorial 2: DX Center](../tutorials/02-dx-center.md)** - Complete example with doped
- **[Getting Started: Quick Start](../getting-started/quick-start.md)** - doped integration example
