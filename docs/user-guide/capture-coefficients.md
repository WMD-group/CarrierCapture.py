# Calculating Capture Coefficients

Complete guide to computing carrier capture rates using configuration coordinate theory.

## Overview

The `ConfigCoordinate` class calculates non-radiative carrier capture coefficients using multiphonon emission theory. Given two potential energy surfaces (initial and final states), it:

1. **Computes wavefunction overlaps** (Franck-Condon factors)
2. **Applies energy-conserving delta function** (with Gaussian broadening)
3. **Sums over phonon states** (Boltzmann-weighted)
4. **Calculates temperature-dependent capture rate** C(T)

---

## Quick Start

```python
from carriercapture.core import ConfigCoordinate
import numpy as np

# Assume pot_i and pot_f are already created, fitted, and solved
# (see Potentials guide)

# Create configuration coordinate
cc = ConfigCoordinate(
    pot_i=pot_initial,  # Initial state (e.g., neutral defect)
    pot_f=pot_final,    # Final state (e.g., charged defect)
    W=0.205,            # Electron-phonon coupling (eV)
    degeneracy=1        # Degeneracy factor
)

# Calculate overlaps
cc.calculate_overlap(
    Q0=10.0,       # Configuration coordinate shift (amu^0.5·Å)
    cutoff=0.25,   # Energy cutoff (eV)
    sigma=0.025    # Gaussian delta width (eV)
)

# Calculate capture coefficient
temperature = np.linspace(100, 500, 50)  # K
cc.calculate_capture_coefficient(
    volume=1e-21,          # Supercell volume (cm³)
    temperature=temperature
)

# Results
print(f"C(300K) = {cc.capture_coefficient[20]:.3e} cm³/s")
```

---

## Theory Background

### Capture Coefficient Formula

The capture rate follows from Fermi's Golden Rule:

$$
C(T) = \frac{2\pi g}{\hbar V} \sum_{i,f} P_i |\langle \psi_f | \hat{Q} | \psi_i \rangle|^2 W^2 \delta(E_i - E_f)
$$

Where:
- $g$ = degeneracy factor
- $V$ = supercell volume
- $P_i$ = Boltzmann occupation of initial state $i$
- $\langle \psi_f | \hat{Q} | \psi_i \rangle$ = overlap integral (Franck-Condon factor)
- $W$ = electron-phonon coupling matrix element
- $\delta(E_i - E_f)$ = energy-conserving delta function

### Key Concepts

**1. Franck-Condon Principle**
Electronic transition is fast compared to nuclear motion → vertical transition in configuration coordinate diagram

**2. Overlap Integrals**
Measure wavefunction overlap at shifted coordinate Q₀:
$$
\langle \psi_f | \hat{Q} | \psi_i \rangle = \int \psi_f(Q) \cdot Q \cdot \psi_i(Q - Q_0) \, dQ
$$

**3. Energy Conservation**
Only transitions conserving energy contribute (broadened by $\sigma$):
$$
\delta_\sigma(E_i - E_f) = \frac{1}{\sqrt{2\pi}\sigma} \exp\left(-\frac{(E_i - E_f)^2}{2\sigma^2}\right)
$$

**4. Boltzmann Weighting**
Temperature dependence from initial state population:
$$
P_i(T) = \frac{e^{-E_i/k_B T}}{Z(T)}, \quad Z(T) = \sum_i e^{-E_i/k_B T}
$$

---

## Step-by-Step Workflow

### Step 1: Prepare Potentials

Both potentials must be fitted and solved:

```python
from carriercapture.core import Potential
from carriercapture.io import load_potential_from_file

# Load solved potentials
pot_i = Potential.from_dict(load_potential_from_file('excited_solved.json'))
pot_f = Potential.from_dict(load_potential_from_file('ground_solved.json'))

# Verify they're ready
assert pot_i.eigenvalues is not None, "Initial potential not solved!"
assert pot_f.eigenvalues is not None, "Final potential not solved!"

print(f"Initial state: {len(pot_i.eigenvalues)} states")
print(f"Final state: {len(pot_f.eigenvalues)} states")
```

**Requirements:**
- Both potentials must have `.eigenvalues` and `.eigenvectors`
- Potentials should be on compatible Q grids (same range and spacing)

### Step 2: Create ConfigCoordinate

```python
from carriercapture.core import ConfigCoordinate

cc = ConfigCoordinate(
    pot_i=pot_i,
    pot_f=pot_f,
    W=0.205,        # Coupling strength (eV)
    degeneracy=1    # g = 1 for non-degenerate states
)
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `pot_i` | Potential | Initial state (before capture) |
| `pot_f` | Potential | Final state (after capture) |
| `W` | float | Electron-phonon coupling (eV) |
| `degeneracy` | int | Degeneracy factor $g$ |

**Determining W:**
- From theory: $W = \langle \Psi_e | \partial V/\partial Q | \Psi_h \rangle$
- Typical values: 0.1 - 0.5 eV
- Can be fitted to experimental data
- Default if unknown: 0.2 eV

### Step 3: Calculate Overlaps

```python
cc.calculate_overlap(
    Q0=10.0,       # Configuration coordinate shift
    cutoff=0.25,   # Energy cutoff (eV)
    sigma=0.025    # Gaussian broadening (eV)
)

# Check results
print(f"Overlap matrix shape: {cc.overlap_matrix.shape}")
print(f"Non-zero elements: {np.count_nonzero(cc.overlap_matrix)}")
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `Q0` | float | required | Coordinate shift (amu^0.5·Å) - **most critical parameter** |
| `cutoff` | float | 0.25 | Energy cutoff for overlaps (eV) |
| `sigma` | float | 0.025 | Gaussian delta width (eV) |

**Choosing Q₀:**

Q₀ represents the shift in configuration coordinate between initial and final states.

```python
# Option 1: From DFT relaxed structures
# Q0 ≈ sqrt(M) * Δr, where M is reduced mass, Δr is displacement

# Option 2: From potential minima
Q_min_i = pot_i.Q[np.argmin(pot_i.E)]
Q_min_f = pot_f.Q[np.argmin(pot_f.E)]
Q0_estimate = abs(Q_min_f - Q_min_i)
print(f"Estimated Q0: {Q0_estimate:.2f} amu^0.5·Å")

# Option 3: Scan over Q0 values
Q0_values = np.linspace(5, 15, 11)
for Q0 in Q0_values:
    cc_temp = ConfigCoordinate(pot_i, pot_f, W=0.205, degeneracy=1)
    cc_temp.calculate_overlap(Q0=Q0, cutoff=0.25, sigma=0.025)
    cc_temp.calculate_capture_coefficient(volume=1e-21, temperature=np.array([300]))
    print(f"Q0={Q0:5.1f}: C(300K)={cc_temp.capture_coefficient[0]:.3e} cm³/s")
```

**Energy Cutoff:**

Limits which transitions are included based on energy conservation:

```python
# Only calculate overlaps for |E_i - E_f| < cutoff
# Reduces computational cost, focuses on relevant transitions

# Typical values:
# - cutoff = 0.25 eV (standard, ~10 kT at 300K)
# - cutoff = 0.5 eV (include more transitions)
# - cutoff = 0.1 eV (strict energy conservation)

# Check energy differences
E_diffs = pot_i.eigenvalues[:, None] - pot_f.eigenvalues[None, :]
print(f"Energy differences range: {E_diffs.min():.3f} to {E_diffs.max():.3f} eV")
print(f"Within cutoff: {np.sum(np.abs(E_diffs) < 0.25)} out of {E_diffs.size}")
```

**Gaussian Width (σ):**

Controls broadening of energy-conserving delta function:

```python
# Smaller σ → sharper energy conservation
# Larger σ → more transitions included

# Typical values:
# - sigma = 0.025 eV (standard, ~1 kT at 300K)
# - sigma = 0.01 eV (strict)
# - sigma = 0.05 eV (broader)

# Rule of thumb: σ ≈ average phonon spacing
spacing = np.diff(pot_i.eigenvalues).mean()
sigma_suggested = spacing / 2
print(f"Suggested σ: {sigma_suggested:.4f} eV")
```

### Step 4: Calculate Capture Coefficient

```python
import numpy as np

# Define temperature range
temperature = np.linspace(100, 500, 50)  # 100-500 K in 50 points

# Calculate C(T)
cc.calculate_capture_coefficient(
    volume=1e-21,          # Supercell volume (cm³)
    temperature=temperature
)

# Access results
C_T = cc.capture_coefficient  # cm³/s
T = cc.temperature            # K

# Print key values
print(f"C(100K) = {C_T[0]:.3e} cm³/s")
print(f"C(300K) = {C_T[len(T)//2]:.3e} cm³/s")
print(f"C(500K) = {C_T[-1]:.3e} cm³/s")
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `volume` | float | Supercell volume (cm³) |
| `temperature` | array | Temperature array (K) |

**Determining Volume:**

The supercell volume V converts the rate from quantum states to concentration units:

```python
# Option 1: From DFT supercell
from pymatgen.core import Structure
structure = Structure.from_file('POSCAR')
volume_ang3 = structure.volume  # Ų
volume_cm3 = volume_ang3 * 1e-24  # cm³
print(f"Supercell volume: {volume_cm3:.3e} cm³")

# Option 2: Typical values
# 2×2×2 supercell of ZnO: ~1e-21 cm³
# 3×3×3 supercell: ~3e-21 cm³
# 4×4×4 supercell: ~8e-21 cm³

# Option 3: From lattice parameters
a, b, c = 5.2, 5.2, 5.2  # Å (cubic example)
volume_cm3 = (a * b * c) * 1e-24
```

---

## Analyzing Results

### Arrhenius Plot

```python
import matplotlib.pyplot as plt

# Arrhenius plot: log(C) vs 1000/T
fig, ax = plt.subplots(figsize=(8, 6))

x = 1000.0 / temperature
y = np.log10(cc.capture_coefficient)

ax.plot(x, y, 'o-', linewidth=2, markersize=4)
ax.set_xlabel('1000/T (K$^{-1}$)', fontsize=12)
ax.set_ylabel('log$_{10}$(C) [cm$^3$/s]', fontsize=12)
ax.set_title('Capture Coefficient (Arrhenius Plot)', fontsize=14)
ax.grid(True, alpha=0.3)

# Add temperature labels on top
ax2 = ax.twiny()
temps = [100, 200, 300, 400, 500]
ax2.set_xlim(ax.get_xlim())
ax2.set_xticks([1000/t for t in temps])
ax2.set_xticklabels([f'{t}K' for t in temps])

plt.tight_layout()
plt.show()
```

Or use built-in visualization:

```python
from carriercapture.visualization import plot_capture_coefficient

fig = plot_capture_coefficient(
    cc,
    arrhenius=True,
    show_activation_energy=True
)
fig.show()
```

### Temperature Dependence

```python
# Identify temperature dependence
T_low = temperature[temperature < 200]
T_high = temperature[temperature > 350]

C_low = cc.capture_coefficient[temperature < 200]
C_high = cc.capture_coefficient[temperature > 350]

# Low-T: typically C ∝ exp(-E_a/kT) if activated
# High-T: may plateau or decrease if barrier-less

# Extract apparent activation energy (if activated)
if len(T_low) > 5:
    x = 1.0 / T_low
    y = np.log(C_low)
    # Linear fit: ln(C) = ln(C0) - Ea/kT
    p = np.polyfit(x, y, 1)
    E_a = -p[0] * 8.617e-5  # eV
    print(f"Apparent activation energy: {E_a:.3f} eV")
```

### Sensitivity Analysis

Check sensitivity to key parameters:

```python
# Vary Q0
Q0_values = [8, 10, 12]
results_Q0 = {}

for Q0 in Q0_values:
    cc_temp = ConfigCoordinate(pot_i, pot_f, W=0.205, degeneracy=1)
    cc_temp.calculate_overlap(Q0=Q0, cutoff=0.25, sigma=0.025)
    cc_temp.calculate_capture_coefficient(volume=1e-21, temperature=np.array([300]))
    results_Q0[Q0] = cc_temp.capture_coefficient[0]

print("\nQ0 sensitivity:")
for Q0, C in results_Q0.items():
    print(f"  Q0 = {Q0:4.1f} amu^0.5·Å: C(300K) = {C:.3e} cm³/s")

# Vary W
W_values = [0.15, 0.205, 0.25]
results_W = {}

for W in W_values:
    cc_temp = ConfigCoordinate(pot_i, pot_f, W=W, degeneracy=1)
    cc_temp.calculate_overlap(Q0=10.0, cutoff=0.25, sigma=0.025)
    cc_temp.calculate_capture_coefficient(volume=1e-21, temperature=np.array([300]))
    results_W[W] = cc_temp.capture_coefficient[0]

print("\nW sensitivity:")
for W, C in results_W.items():
    print(f"  W = {W:.3f} eV: C(300K) = {C:.3e} cm³/s")
```

---

## Advanced Usage

### Inspecting Overlap Matrix

```python
# Overlap matrix has shape (n_initial, n_final)
overlap = cc.overlap_matrix
print(f"Shape: {overlap.shape}")

# Find dominant transitions
max_overlap = np.max(np.abs(overlap))
i_max, f_max = np.unravel_index(np.argmax(np.abs(overlap)), overlap.shape)

print(f"Maximum overlap: {max_overlap:.6f}")
print(f"  Between initial state {i_max} (E={pot_i.eigenvalues[i_max]:.4f} eV)")
print(f"  and final state {f_max} (E={pot_f.eigenvalues[f_max]:.4f} eV)")

# Visualize overlap matrix
from carriercapture.visualization import plot_overlap_matrix

fig = plot_overlap_matrix(cc, log_scale=True)
fig.show()
```

### Configuration Coordinate Diagram

```python
from carriercapture.visualization import plot_configuration_coordinate

fig = plot_configuration_coordinate(
    pot_i,
    pot_f,
    Q0=10.0,
    show_crossing=True,
    show_wavefunctions=True,
    n_states=10
)
fig.show()
```

### Custom Temperature Ranges

```python
# Low temperature only
T_low = np.linspace(10, 100, 30)
cc.calculate_capture_coefficient(volume=1e-21, temperature=T_low)

# High temperature only
T_high = np.linspace(300, 1000, 50)
cc.calculate_capture_coefficient(volume=1e-21, temperature=T_high)

# Non-uniform spacing (focus on interesting region)
T_custom = np.concatenate([
    np.linspace(100, 250, 10),   # Sparse at low T
    np.linspace(250, 350, 30),   # Dense near 300K
    np.linspace(350, 500, 10)    # Sparse at high T
])
cc.calculate_capture_coefficient(volume=1e-21, temperature=T_custom)
```

### Multiple Degeneracy Factors

```python
# Example: spin degeneracy for different transitions
# g=1: singlet-singlet
# g=2: doublet-doublet
# g=3: triplet-triplet

degeneracies = [1, 2, 3]
results = {}

for g in degeneracies:
    cc_temp = ConfigCoordinate(pot_i, pot_f, W=0.205, degeneracy=g)
    cc_temp.calculate_overlap(Q0=10.0, cutoff=0.25, sigma=0.025)
    cc_temp.calculate_capture_coefficient(volume=1e-21, temperature=np.array([300]))
    results[g] = cc_temp.capture_coefficient[0]

print("\nDegeneracy effect:")
for g, C in results.items():
    print(f"  g = {g}: C(300K) = {C:.3e} cm³/s")

# Note: C scales linearly with g
```

---

## Best Practices

### 1. Convergence Checks

```python
# Check nev convergence
nev_values = [40, 60, 80, 100]
C_values = []

for nev in nev_values:
    # Re-solve with different nev
    pot_f_temp = Potential.from_dict(load_potential_from_file('ground_fitted.json'))
    pot_f_temp.solve(nev=nev)

    cc_temp = ConfigCoordinate(pot_i, pot_f_temp, W=0.205, degeneracy=1)
    cc_temp.calculate_overlap(Q0=10.0)
    cc_temp.calculate_capture_coefficient(volume=1e-21, temperature=np.array([300]))
    C_values.append(cc_temp.capture_coefficient[0])

print("\nnev convergence:")
for nev, C in zip(nev_values, C_values):
    print(f"  nev = {nev:3d}: C(300K) = {C:.3e} cm³/s")

# Check relative change
for i in range(1, len(C_values)):
    rel_change = abs(C_values[i] - C_values[i-1]) / C_values[i-1]
    print(f"  {nev_values[i-1]}→{nev_values[i]}: {rel_change*100:.2f}% change")
```

### 2. Physical Reasonableness

```python
# Typical capture coefficient ranges
# Fast: 1e-6 to 1e-8 cm³/s
# Moderate: 1e-10 to 1e-12 cm³/s
# Slow: 1e-14 to 1e-16 cm³/s

C_300 = cc.capture_coefficient[temperature == 300][0]

if C_300 > 1e-6:
    print(f"⚠️  C(300K) = {C_300:.3e} cm³/s is very fast - check parameters")
elif C_300 < 1e-20:
    print(f"⚠️  C(300K) = {C_300:.3e} cm³/s is very slow - check parameters")
else:
    print(f"✓ C(300K) = {C_300:.3e} cm³/s is reasonable")
```

### 3. Comparing with Experiments

```python
# If experimental data available
T_exp = np.array([200, 250, 300, 350, 400])  # K
C_exp = np.array([1e-12, 5e-12, 2e-11, 7e-11, 2e-10])  # cm³/s

# Calculate at same temperatures
cc.calculate_capture_coefficient(volume=1e-21, temperature=T_exp)
C_calc = cc.capture_coefficient

# Compare
fig, ax = plt.subplots()
ax.plot(1000/T_exp, np.log10(C_exp), 'o', label='Experiment', markersize=8)
ax.plot(1000/T_exp, np.log10(C_calc), 's', label='Calculation', markersize=8)
ax.set_xlabel('1000/T (K$^{-1}$)')
ax.set_ylabel('log$_{10}$(C) [cm$^3$/s]')
ax.legend()
ax.grid(True, alpha=0.3)
plt.show()

# Quantify agreement
rmse = np.sqrt(np.mean((np.log10(C_calc) - np.log10(C_exp))**2))
print(f"RMSE (log scale): {rmse:.3f}")
```

---

## Troubleshooting

### Problem: C(T) is too large/small

**Possible causes:**
1. **Wrong volume** → Check supercell volume units (should be cm³, typically ~1e-21)
2. **Wrong W** → Adjust electron-phonon coupling
3. **Wrong Q₀** → Scan over Q₀ values
4. **Insufficient nev** → Increase number of eigenvalues

**Debug:**
```python
print(f"Volume: {cc.volume:.3e} cm³ (typical: 1e-21)")
print(f"W: {cc.W:.3f} eV (typical: 0.1-0.5)")
print(f"Q0: {cc.Q0:.2f} amu^0.5·Å")
print(f"States: {len(pot_i.eigenvalues)} initial, {len(pot_f.eigenvalues)} final")
```

### Problem: Strange temperature dependence

**Possible causes:**
1. **Too few phonon states** → Increase nev
2. **Wrong sigma** → Try sigma = average level spacing / 2
3. **Potentials on incompatible grids** → Check Q grids match

**Debug:**
```python
# Check if Q grids are compatible
print(f"Initial Q: {pot_i.Q.min():.2f} to {pot_i.Q.max():.2f} ({len(pot_i.Q)} points)")
print(f"Final Q: {pot_f.Q.min():.2f} to {pot_f.Q.max():.2f} ({len(pot_f.Q)} points)")

if not np.allclose(pot_i.Q, pot_f.Q):
    print("⚠️  Warning: Q grids do not match!")
```

### Problem: Very few non-zero overlaps

**Possible causes:**
1. **Cutoff too small** → Increase cutoff
2. **Poor energy alignment** → Check potential energy offsets
3. **Q₀ far from optimal** → Scan Q₀

**Debug:**
```python
n_nonzero = np.count_nonzero(cc.overlap_matrix)
n_total = cc.overlap_matrix.size
print(f"Non-zero overlaps: {n_nonzero} / {n_total} ({100*n_nonzero/n_total:.1f}%)")

if n_nonzero < n_total * 0.01:
    print("⚠️  Very sparse overlap matrix - consider:")
    print("   1. Increase cutoff")
    print("   2. Adjust Q0")
    print("   3. Check energy alignment")
```

---

## Complete Example

```python
import numpy as np
from carriercapture.core import Potential, ConfigCoordinate
from carriercapture.io import load_potential_from_file
from carriercapture.visualization import (
    plot_capture_coefficient,
    plot_configuration_coordinate
)

# 1. Load solved potentials
pot_i = Potential.from_dict(load_potential_from_file('excited_solved.json'))
pot_f = Potential.from_dict(load_potential_from_file('ground_solved.json'))
print(f"✓ Loaded potentials")

# 2. Create configuration coordinate
cc = ConfigCoordinate(
    pot_i=pot_i,
    pot_f=pot_f,
    W=0.205,        # eV
    degeneracy=1
)
print(f"✓ Created ConfigCoordinate")

# 3. Calculate overlaps
cc.calculate_overlap(
    Q0=10.0,        # amu^0.5·Å
    cutoff=0.25,    # eV
    sigma=0.025     # eV
)
n_nonzero = np.count_nonzero(cc.overlap_matrix)
print(f"✓ Calculated overlaps: {n_nonzero} non-zero elements")

# 4. Calculate capture coefficient
temperature = np.linspace(100, 500, 50)
cc.calculate_capture_coefficient(
    volume=1e-21,   # cm³
    temperature=temperature
)
print(f"✓ Calculated C(T)")

# 5. Print results
print(f"\nResults:")
print(f"  C(100K) = {cc.capture_coefficient[0]:.3e} cm³/s")
print(f"  C(300K) = {cc.capture_coefficient[24]:.3e} cm³/s")
print(f"  C(500K) = {cc.capture_coefficient[-1]:.3e} cm³/s")

# 6. Visualize
fig1 = plot_capture_coefficient(cc, arrhenius=True)
fig1.show()

fig2 = plot_configuration_coordinate(pot_i, pot_f, Q0=10.0, show_crossing=True)
fig2.show()

print(f"\n✓ Complete!")
```

---

## See Also

- **[API Reference: ConfigCoordinate](../api/core.md#configcoordinate)** - Complete API documentation
- **[Theory: Multiphonon Theory](../theory/multiphonon-theory.md)** - Theoretical background
- **[Theory: Configuration Coordinates](../theory/configuration-coordinates.md)** - CC diagram explanation
- **[Tutorial 1: Sn in ZnO](../tutorials/01-harmonic-oscillator.md)** - Complete example
- **[CLI Reference: capture](../api/cli.md#capture---calculate-capture-coefficient)** - Command-line usage
