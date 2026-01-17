# Working with Potentials

Complete guide to creating, fitting, and solving potential energy surfaces in CarrierCapture.

## Overview

The `Potential` class is the foundation of all carrier capture calculations. It represents 1D potential energy surfaces along the configuration coordinate Q, with capabilities for:

- **Data input** from various sources (files, arrays, theoretical models)
- **Fitting** to smooth functional forms
- **Solving** the 1D Schrödinger equation for phonon states
- **Serialization** for saving and loading

---

## Creating Potentials

### From Data Arrays

The most common way to create a potential is from Q-E data:

```python
import numpy as np
from carriercapture.core import Potential

# Your Q-E data from DFT calculations
Q_data = np.array([0, 2, 4, 6, 8, 10])  # amu^0.5·Å
E_data = np.array([0.0, 0.1, 0.3, 0.5, 0.6, 0.65])  # eV

# Create potential
pot = Potential(Q_data=Q_data, E_data=E_data, name="My Potential")
```

**Parameters:**
- `Q_data` - Configuration coordinates (amu^0.5·Å)
- `E_data` - Energies (eV)
- `name` (optional) - Descriptive name
- `Q0` (optional) - Coordinate offset
- `E0` (optional) - Energy offset

### From File

Load potential data from files:

```python
from carriercapture.io import load_potential, read_potential_data

# Option 1: Read Q-E data from CSV/DAT
Q, E = read_potential_data('potential.csv')
pot = Potential(Q_data=Q, E_data=E)

# Option 2: Load fitted potential from JSON/YAML/NPZ
pot_dict = load_potential('fitted_potential.json')
pot = Potential.from_dict(pot_dict)
```

**Supported input formats:**
- CSV/DAT files (two columns: Q, E)
- JSON/YAML with metadata
- NPZ with NumPy arrays

### From Harmonic Model

Create harmonic oscillator potential analytically:

```python
from carriercapture.core import Potential

# Create harmonic potential: E(Q) = 0.5 * k * (Q - Q0)^2 + E0
# where k = mass * omega^2, and hw = ℏω
pot = Potential.from_harmonic(
    hw=0.008,      # Phonon energy (eV) - 8 meV
    Q0=0.0,        # Equilibrium position (amu^0.5·Å)
    E0=0.0,        # Energy offset (eV)
    Q_range=(-20, 20),  # Range for evaluation
    npoints=3001   # Grid points
)

# Already fitted and ready to solve
print(f"Fit type: {pot.fit_type}")  # 'harmonic'
```

**Use cases:**
- Quick testing and benchmarking
- Initial guess for anharmonic potentials
- Simple defects with harmonic character

### From doped Package

Extract potentials from [doped](https://github.com/SMTG-Bham/doped) defect calculations:

```python
from carriercapture.io.doped_interface import (
    load_defect_entry,
    load_path_calculations,
    create_potential_from_doped
)

# Load DefectEntry
defect = load_defect_entry('Sn_Zn.json.gz')

# Load VASP path calculations for two charge states
Q_i, E_i = load_path_calculations('path_q0/')
Q_f, E_f = load_path_calculations('path_q1/')

# Create potentials
pot_initial = create_potential_from_doped(
    defect, charge_state=0, Q_data=Q_i, E_data=E_i
)
pot_final = create_potential_from_doped(
    defect, charge_state=+1, Q_data=Q_f, E_data=E_f
)
```

See [doped Integration](doped-integration.md) for complete workflow.

---

## Fitting Potentials

Raw Q-E data must be fitted to a smooth function before solving the Schrödinger equation.

### Spline Fitting (Recommended)

Cubic spline interpolation with smoothing - best for anharmonic potentials:

```python
pot = Potential(Q_data=Q_data, E_data=E_data)

# Fit with spline
pot.fit(
    fit_type='spline',
    order=4,           # Spline order (3=cubic, 4=quartic, 5=quintic)
    smoothness=0.001   # Smoothing parameter (0=exact interpolation)
)

# Check fit quality
Q_fit = np.linspace(Q_data.min(), Q_data.max(), 500)
E_fit = pot(Q_fit)  # Evaluate fitted potential

import matplotlib.pyplot as plt
plt.plot(Q_data, E_data, 'o', label='Data')
plt.plot(Q_fit, E_fit, '-', label='Fit')
plt.legend()
plt.show()
```

**Parameters:**
- `order` - Spline order (typically 3-5)
- `smoothness` - Controls smoothing vs. exact fit
  - `0` = Exact interpolation through all points
  - `0.001` = Light smoothing (recommended for DFT data)
  - `0.01` = Heavy smoothing (for noisy data)

**When to use:**
- DFT-calculated potential energy surfaces
- Anharmonic potentials
- Data with moderate noise

### Harmonic Fitting

Fit to harmonic oscillator: $E(Q) = \frac{1}{2}k(Q-Q_0)^2 + E_0$

```python
pot = Potential(Q_data=Q_data, E_data=E_data)

# Fit to harmonic form
pot.fit(fit_type='harmonic')

# Access fit parameters
print(f"Frequency: {pot.fit_params['hw']:.4f} eV")
print(f"Q0: {pot.fit_params['Q0']:.3f} amu^0.5·Å")
print(f"E0: {pot.fit_params['E0']:.3f} eV")
```

**When to use:**
- Near-harmonic potentials
- Quick fits for simple systems
- Initial approximation

### Morse Potential Fitting

Fit to Morse potential: $E(Q) = D_e[1 - e^{-a(Q-Q_e)}]^2$

```python
pot.fit(
    fit_type='morse',
    initial_guess={'De': 1.0, 'a': 0.5, 'Qe': 5.0}  # optional
)

# Access fit parameters
print(f"Dissociation energy: {pot.fit_params['De']:.3f} eV")
print(f"Width parameter: {pot.fit_params['a']:.3f}")
```

**When to use:**
- Bond-breaking processes
- Asymmetric potentials
- Anharmonic wells with well-defined minimum

### Polynomial Fitting

Fit to polynomial: $E(Q) = \sum_{i=0}^{n} c_i Q^i$

```python
pot.fit(
    fit_type='polynomial',
    degree=6  # Polynomial degree
)
```

**When to use:**
- Symmetric potentials near minimum
- Quick analytical form
- Not recommended for large Q ranges (oscillations)

### Comparing Fit Quality

```python
# Try multiple fitting methods
methods = ['spline', 'harmonic', 'polynomial']
fits = {}

for method in methods:
    pot_copy = Potential(Q_data=Q_data, E_data=E_data)
    if method == 'polynomial':
        pot_copy.fit(fit_type=method, degree=6)
    elif method == 'spline':
        pot_copy.fit(fit_type=method, order=4, smoothness=0.001)
    else:
        pot_copy.fit(fit_type=method)
    fits[method] = pot_copy

# Evaluate RMSE
Q_test = np.linspace(Q_data.min(), Q_data.max(), 100)
for method, pot_fit in fits.items():
    E_pred = pot_fit(Q_test)
    # Interpolate data for comparison
    E_true = np.interp(Q_test, Q_data, E_data)
    rmse = np.sqrt(np.mean((E_pred - E_true)**2))
    print(f"{method:12s}: RMSE = {rmse:.6f} eV")
```

---

## Solving the Schrödinger Equation

After fitting, solve the 1D Schrödinger equation to find phonon eigenvalues and eigenvectors.

### Basic Solving

```python
# Solve for 60 eigenvalues
pot.solve(nev=60)

# Access results
eigenvalues = pot.eigenvalues    # Energy levels (eV)
eigenvectors = pot.eigenvectors  # Wavefunctions (normalized)

print(f"Ground state energy: {eigenvalues[0]:.6f} eV")
print(f"First excited state: {eigenvalues[1]:.6f} eV")
print(f"Level spacing: {eigenvalues[1] - eigenvalues[0]:.6f} eV")
```

**Parameters:**
- `nev` - Number of eigenvalues to compute
- `maxiter` (optional) - Maximum iterations for eigensolver

**Typical values:**
- Initial (neutral) state: `nev=180` (requires more states)
- Final (charged) state: `nev=60` (fewer states needed)

### Grid Considerations

The Q grid affects solution accuracy:

```python
# Option 1: Let the potential determine the grid
pot.solve(nev=60)  # Uses existing Q grid

# Option 2: Set custom grid before solving
pot.Q = np.linspace(-15, 15, 3001)  # 3001 points from -15 to 15
pot.E = pot.fit_func(pot.Q)         # Evaluate potential on new grid
pot.solve(nev=60)

# Option 3: Use denser grid for higher accuracy
pot.Q = np.linspace(-20, 20, 5001)  # 5001 points
pot.E = pot.fit_func(pot.Q)
pot.solve(nev=180)
```

**Grid best practices:**
- **Range**: Cover full potential well + classically forbidden regions
- **Density**: At least 2000-3000 points for accurate results
- **Symmetry**: Center grid around potential minimum for symmetric wells

### Accessing Wavefunctions

```python
pot.solve(nev=60)

# Ground state wavefunction
psi_0 = pot.eigenvectors[0, :]  # Shape: (npoints,)
Q = pot.Q

# Plot ground and first excited states
import matplotlib.pyplot as plt

fig, ax = plt.subplots()
ax.plot(Q, pot.E, 'k-', label='Potential', linewidth=2)

for i in range(3):
    psi = pot.eigenvectors[i, :]
    # Scale and shift for visualization
    psi_scaled = 0.2 * psi + pot.eigenvalues[i]
    ax.plot(Q, psi_scaled, label=f'ψ_{i}')
    ax.axhline(pot.eigenvalues[i], linestyle='--', alpha=0.3)

ax.set_xlabel('Q (amu$^{0.5}$·Å)')
ax.set_ylabel('E (eV)')
ax.legend()
plt.show()
```

### Checking Solver Convergence

```python
# Solve and check spacing
pot.solve(nev=60)

spacings = np.diff(pot.eigenvalues)
print(f"Average spacing: {spacings.mean():.6f} eV")
print(f"Min spacing: {spacings.min():.6f} eV")
print(f"Max spacing: {spacings.max():.6f} eV")

# For harmonic oscillator, spacing should be constant ≈ ℏω
# Increasing spacing → anharmonicity
```

---

## Saving and Loading

### Save Potential

```python
from carriercapture.io import save_potential

# Save to JSON (human-readable)
save_potential(pot, 'potential.json')

# Save to NPZ (compact, fast)
save_potential(pot, 'potential.npz')

# Save to YAML (human-readable, includes metadata)
save_potential(pot, 'potential.yaml')
```

### Load Potential

```python
from carriercapture.io import load_potential_from_file

# Load from any format
data = load_potential_from_file('potential.json')
pot = Potential.from_dict(data)

# Check what's included
print(f"Has fit: {pot.fit_func is not None}")
print(f"Has eigenvalues: {pot.eigenvalues is not None}")
```

### Serialization Workflow

```python
# Complete workflow with saving
from carriercapture.core import Potential
from carriercapture.io import save_potential, load_potential_from_file

# 1. Create and fit
pot = Potential(Q_data=Q_data, E_data=E_data, name="Excited State")
pot.fit(fit_type='spline', order=4, smoothness=0.001)
save_potential(pot, 'excited_fitted.json')

# 2. Load and solve
pot = Potential.from_dict(load_potential_from_file('excited_fitted.json'))
pot.solve(nev=180)
save_potential(pot, 'excited_solved.json')

# 3. Use in capture calculation
pot = Potential.from_dict(load_potential_from_file('excited_solved.json'))
# Ready for ConfigCoordinate calculation
```

---

## Best Practices

### Choosing Fit Type

| Situation | Recommended Fit | Reason |
|-----------|----------------|--------|
| DFT data, anharmonic | `spline` | Flexible, accurate |
| DFT data, near-harmonic | `spline` or `harmonic` | Compare both |
| Noisy data | `spline` with smoothness~0.01 | Reduces noise |
| Symmetric well | `polynomial` or `spline` | Both work well |
| Bond dissociation | `morse` | Physical form |
| Quick testing | `harmonic` | Analytical |

### Choosing Number of Eigenvalues

```python
# Rule of thumb based on energy range
E_range = E_data.max() - E_data.min()  # Full potential depth

# Initial state (typically needs more)
nev_initial = int(E_range / 0.01) + 50  # ~50-200 states

# Final state (typically needs fewer)
nev_final = int(E_range / 0.02) + 30   # ~30-80 states

# Minimum values
nev_initial = max(nev_initial, 60)
nev_final = max(nev_final, 30)

print(f"Recommended nev_i: {nev_initial}")
print(f"Recommended nev_f: {nev_final}")
```

### Checking Grid Adequacy

```python
# Check that highest eigenvalue is well below potential maximum
pot.solve(nev=60)

E_max = pot.eigenvalues[-1]
V_max = pot.E.max()

if E_max > 0.8 * V_max:
    print(f"⚠️  Warning: Highest eigenvalue ({E_max:.3f} eV) is too close to")
    print(f"   potential maximum ({V_max:.3f} eV).")
    print(f"   → Increase Q range or reduce nev")
else:
    print(f"✓ Grid adequate: E_max = {E_max:.3f} eV, V_max = {V_max:.3f} eV")
```

### Troubleshooting

**Problem: `ArpackNoConvergence` error**
```python
# Solution 1: Reduce nev
pot.solve(nev=40)  # Instead of 60

# Solution 2: Increase maxiter
pot.solve(nev=60, maxiter=10000)

# Solution 3: Improve grid
pot.Q = np.linspace(-20, 20, 5001)  # Denser grid
pot.E = pot.fit_func(pot.Q)
pot.solve(nev=60)
```

**Problem: Unphysical eigenvalues**
```python
# Check for negative spacing
spacings = np.diff(pot.eigenvalues)
if np.any(spacings < 0):
    print("⚠️  Eigenvalues not monotonic - check fit quality")

# Visualize potential and eigenvalues
from carriercapture.visualization import plot_potential
fig = plot_potential(pot, show_wavefunctions=True)
fig.show()
```

**Problem: Fit oscillates**
```python
# Use more smoothing for spline
pot.fit(fit_type='spline', order=4, smoothness=0.01)

# Or reduce polynomial degree
pot.fit(fit_type='polynomial', degree=4)  # Instead of 6
```

---

## Advanced Usage

### Custom Potential Function

```python
# Define custom potential
def double_well(Q, a=1.0, b=0.5):
    """Quartic double-well potential."""
    return a * Q**4 - b * Q**2

# Create potential
Q = np.linspace(-3, 3, 3001)
E = double_well(Q, a=0.1, b=1.0)

pot = Potential(Q_data=Q, E_data=E, name="Double Well")
pot.fit(fit_type='spline')
pot.solve(nev=100)
```

### Potential Arithmetic

```python
# Shift potential energy
pot_shifted = pot.copy()
pot_shifted.E = pot_shifted.E + 0.5  # Add 0.5 eV

# Shift coordinate
pot_shifted = pot.copy()
pot_shifted.Q = pot_shifted.Q - 5.0  # Shift by 5 amu^0.5·Å
```

### Extrapolation Warning

```python
# Potentials should not be evaluated outside fitted range
Q_fit_range = (pot.Q.min(), pot.Q.max())

# Don't do this:
Q_extrap = np.linspace(-30, 30, 1000)  # Outside fit range!
E_extrap = pot(Q_extrap)  # May give unphysical results

# Instead, refit with wider range:
pot.Q = Q_extrap
pot.E = <your_data_or_model>
pot.fit(fit_type='spline')
```

---

## Complete Example

```python
import numpy as np
from carriercapture.core import Potential
from carriercapture.io import save_potential
from carriercapture.visualization import plot_potential

# 1. Load data (from DFT or model)
Q_data = np.linspace(0, 20, 15)
E_data = 0.5 * 0.008 * (Q_data - 10)**2  # Harmonic example

# 2. Create potential
pot = Potential(
    Q_data=Q_data,
    E_data=E_data,
    name="Example: Sn_Zn Ground State"
)

# 3. Fit potential
pot.fit(fit_type='spline', order=4, smoothness=0.001)
print(f"✓ Fitted with {pot.fit_type} method")

# 4. Solve Schrödinger equation
pot.solve(nev=60)
print(f"✓ Solved for {len(pot.eigenvalues)} states")
print(f"  Ground state: {pot.eigenvalues[0]:.6f} eV")
print(f"  Level spacing: {pot.eigenvalues[1] - pot.eigenvalues[0]:.6f} eV")

# 5. Visualize
fig = plot_potential(
    pot,
    show_wavefunctions=True,
    n_states=10
)
fig.show()

# 6. Save for later use
save_potential(pot, 'ground_state_solved.json')
print(f"✓ Saved to ground_state_solved.json")
```

---

## See Also

- **[API Reference: Potential](../api/core.md#potential)** - Complete API documentation
- **[CLI Reference: fit](../api/cli.md#fit---fit-potential-energy-surface)** - Fitting from command line
- **[CLI Reference: solve](../api/cli.md#solve---solve-schr%C3%B6dinger-equation)** - Solving from command line
- **[Theory: Configuration Coordinates](../theory/configuration-coordinates.md)** - Theoretical background
- **[Tutorial 1: Harmonic Oscillator](../tutorials/01-harmonic-oscillator.md)** - Step-by-step tutorial
