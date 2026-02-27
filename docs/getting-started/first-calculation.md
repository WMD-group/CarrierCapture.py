# First Calculation: Step-by-Step

This tutorial walks you through a complete carrier capture calculation from start to finish.

## Overview

We'll calculate the capture coefficient for the **Sn_Zn** defect in ZnO - a well-studied system from the literature.

**System**: Sn substituting Zn in ZnO
**Transition**: Neutral (0) → Positive (+1) charge state
**Method**: Harmonic approximation

## Prerequisites

Ensure CarrierCapture is installed:

```bash
python -c "import carriercapture; print('✓ CarrierCapture installed')"
```

## Step 1: Import Libraries

```python
import numpy as np
import matplotlib.pyplot as plt
from carriercapture.core import Potential, ConfigCoordinate
from carriercapture.visualization import plot_potential, plot_capture_coefficient
```

## Step 2: Create Potential Energy Surfaces

We'll create two displaced harmonic potentials representing the neutral and charged states.

### Initial State (Neutral, q=0)

```python
# Excited state (higher energy)
pot_initial = Potential.from_harmonic(
    hw=0.008,    # 8 meV phonon (typical for oxides)
    Q0=0.0,      # Centered at origin
    E0=0.5,      # 0.5 eV above ground state
    Q_range=(-30, 30),  # Wide range for wavefunctions
    npoints=5000         # Fine grid for accuracy
)

print("Initial potential created:")
print(f"  Phonon energy: {0.008*1000:.1f} meV")
print(f"  Equilibrium: Q₀ = {0.0} amu^0.5·Å")
print(f"  Energy offset: {0.5} eV")
```

### Final State (Charged, q=+1)

```python
# Ground state (lower energy, displaced)
pot_final = Potential.from_harmonic(
    hw=0.008,    # Same phonon frequency
    Q0=10.5,     # Displaced by 10.5 amu^0.5·Å
    E0=0.0,      # Ground state reference
    Q_range=(-30, 30),
    npoints=5000
)

print("Final potential created:")
print(f"  Displacement: ΔQ = {10.5} amu^0.5·Å")
print(f"  Energy difference: ΔE = {0.5} eV")
```

**Physical interpretation:**

- **ΔQ = 10.5**: Significant lattice relaxation upon charge state change
- **ΔE = 0.5 eV**: Transition level between charge states

## Step 3: Solve Schrödinger Equation

Find the vibrational states (eigenvalues and wavefunctions) for each potential.

```python
# Solve for initial state
print("\nSolving initial state...")
pot_initial.solve(nev=180, maxiter=90000, tol=0)

print(f"  Found {len(pot_initial.eigenvalues)} eigenvalues")
print(f"  Ground state: ε₀ = {pot_initial.eigenvalues[0]:.4f} eV")
print(f"  First excited: ε₁ = {pot_initial.eigenvalues[1]:.4f} eV")
print(f"  Spacing: Δε = {pot_initial.eigenvalues[1] - pot_initial.eigenvalues[0]:.4f} eV")
print(f"  Expected (ℏω): {0.008:.4f} eV ✓")

# Solve for final state
print("\nSolving final state...")
pot_final.solve(nev=60, maxiter=30000, tol=0)

print(f"  Found {len(pot_final.eigenvalues)} eigenvalues")
print(f"  Ground state: ε₀ = {pot_final.eigenvalues[0]:.4f} eV")
```

**What's happening:**

- ARPACK eigenvalue solver finds lowest 180 (initial) and 60 (final) states
- For harmonic oscillator: $\varepsilon_n = E_0 + \hbar\omega(n + \frac{1}{2})$
- Wavefunctions are Hermite polynomials

**Why different nev?**

- Initial state is higher in energy → more thermally accessible states
- Need enough states for partition function convergence

## Step 4: Visualize Potentials

Let's see what we created:

```python
# Plot both potentials with wavefunctions
fig = plot_potential(
    pot_initial,
    title="Initial State (Neutral)",
    show_wavefunctions=True,
    n_states=10  # Show first 10 states
)
fig.show()

fig = plot_potential(
    pot_final,
    title="Final State (Charged)",
    show_wavefunctions=True,
    n_states=10
)
fig.show()
```

**What to observe:**

- Ground state wavefunction is Gaussian-like (n=0)
- Higher states have more nodes (n=1 has 1 node, n=2 has 2, ...)
- Wavefunctions extend further for higher energies
- Displacement between potentials is visible

## Step 5: Calculate Overlap Matrix

Calculate overlaps between initial and final state wavefunctions.

```python
# Create ConfigCoordinate object
cc = ConfigCoordinate(
    pot_i=pot_initial,
    pot_f=pot_final,
    W=0.068,  # Electron-phonon coupling (eV)
    name="Sn_Zn"
)

# Calculate overlap matrix
print("\nCalculating overlaps...")
cc.calculate_overlap(
    Q0=5.0,      # Coupling coordinate (midpoint)
    cutoff=0.25,  # Energy cutoff for delta function (eV)
    sigma=0.01    # Gaussian broadening (eV)
)

print(f"  Overlap matrix shape: {cc.overlap_matrix.shape}")
print(f"  Non-zero overlaps: {np.count_nonzero(cc.overlap_matrix)}")
print(f"  Max overlap: {np.max(np.abs(cc.overlap_matrix)):.4f}")
```

**Parameters explained:**

- **Q0**: Point where coupling is evaluated (often midpoint between minima)
- **cutoff**: Only include transitions within this energy range
- **sigma**: Broadening for energy-conserving delta function

**Physical meaning:**

- Large overlaps → easy transitions
- Energy cutoff enforces approximate energy conservation
- Gaussian broadening accounts for finite phonon lifetime

## Step 6: Calculate Capture Coefficient

Finally, calculate how capture coefficient varies with temperature.

```python
# Temperature range: 100 to 500 K
temperatures = np.linspace(100, 500, 50)

print("\nCalculating capture coefficient...")
cc.calculate_capture_coefficient(
    volume=1e-21,      # Supercell volume (cm³)
    temperature=temperatures,
    degeneracy=1       # No spin degeneracy
)

print(f"  Calculated for {len(temperatures)} temperatures")
print(f"  Range: {temperatures[0]:.0f} - {temperatures[-1]:.0f} K")
```

**What's being calculated:**

$$C(T) = \frac{V \cdot 2\pi}{\hbar} \cdot g \cdot W^2 \cdot \sum_{i,j} p_i |S_{ij}|^2 \delta(\varepsilon_i - \varepsilon_j)$$

- Sum over all possible $i \to j$ transitions
- Weight by thermal occupation $p_i$
- Include overlap factors $S_{ij}$
- Enforce energy conservation $\delta(\varepsilon_i - \varepsilon_j)$

## Step 7: Analyze Results

### Print Key Values

```python
# Find indices for specific temperatures
idx_100K = 0
idx_300K = np.argmin(np.abs(temperatures - 300))
idx_500K = -1

print("\n" + "="*50)
print("RESULTS: Capture Coefficient")
print("="*50)
print(f"  T = 100 K: C = {cc.capture_coefficient[idx_100K]:.3e} cm³/s")
print(f"  T = 300 K: C = {cc.capture_coefficient[idx_300K]:.3e} cm³/s")
print(f"  T = 500 K: C = {cc.capture_coefficient[idx_500K]:.3e} cm³/s")
print("="*50)

# Calculate Arrhenius slope (if linear region exists)
# C ~ exp(-Ea/kT) → ln(C) = const - Ea/kT
log_C = np.log(cc.capture_coefficient)
inv_T = 1 / temperatures

# Fit linear region (e.g., 200-400 K)
mask = (temperatures >= 200) & (temperatures <= 400)
slope, intercept = np.polyfit(inv_T[mask], log_C[mask], 1)
Ea = -slope * 8.617e-5  # Convert to eV (kB in eV/K)

print(f"\nArrhenius analysis (200-400 K):")
print(f"  Activation energy: Ea ≈ {Ea:.3f} eV")
```

### Plot Temperature Dependence

```python
# Arrhenius plot (log C vs 1000/T)
fig = plot_capture_coefficient(
    cc,
    title="Sn_Zn Capture Coefficient",
    arrhenius=True  # Log scale
)
fig.show()
```

**What to look for:**

1. **Low T region**: Exponential increase (thermally activated)
2. **High T region**: May plateau or decrease (multiphonon cascade)
3. **Arrhenius slope**: Related to activation barrier

## Step 8: Physical Interpretation

### Temperature Dependence

```python
# Ratio: how much faster at 500K vs 100K?
ratio = cc.capture_coefficient[-1] / cc.capture_coefficient[0]
print(f"\nCapture is {ratio:.1e}× faster at 500K than 100K")
```

### Typical values for Sn_Zn

Expected results (literature comparison):

- **100 K**: ~10⁻²⁰ cm³/s (very slow, frozen)
- **300 K**: ~10⁻¹⁴ cm³/s (moderate, typical for defects)
- **500 K**: ~10⁻¹¹ cm³/s (fast, thermally enhanced)

### Defect Lifetime

If you know the carrier density:

```python
n_carrier = 1e16  # cm⁻³ (typical doping)
C_300K = cc.capture_coefficient[idx_300K]

capture_rate = C_300K * n_carrier  # s⁻¹
lifetime = 1 / capture_rate  # s

print(f"\nAt 300K with n = {n_carrier:.1e} cm⁻³:")
print(f"  Capture rate: {capture_rate:.3e} s⁻¹")
print(f"  Defect lifetime: {lifetime:.3e} s")
```

## Step 9: Save Results

### Save to File

```python
# Save calculation for later use
cc.to_dict()  # Convert to dictionary

# Or export specific data
np.savetxt('capture_vs_temperature.txt',
           np.column_stack([temperatures, cc.capture_coefficient]),
           header='Temperature(K)  CaptureCoefficient(cm3/s)')
```

### Export Figure

```python
fig = plot_capture_coefficient(cc)
fig.write_html('capture_coefficient.html')
fig.write_image('capture_coefficient.png', width=800, height=600)
```

## Complete Code

Here's the full calculation in one script:

```python
import numpy as np
from carriercapture.core import Potential, ConfigCoordinate
from carriercapture.visualization import plot_capture_coefficient

# 1. Create potentials
pot_initial = Potential.from_harmonic(hw=0.008, Q0=0.0, E0=0.5)
pot_final = Potential.from_harmonic(hw=0.008, Q0=10.5, E0=0.0)

# 2. Solve Schrödinger equation
pot_initial.solve(nev=180)
pot_final.solve(nev=60)

# 3. Calculate capture coefficient
cc = ConfigCoordinate(pot_i=pot_initial, pot_f=pot_final, W=0.068)
cc.calculate_overlap(Q0=5.0, cutoff=0.25, sigma=0.01)

temperatures = np.linspace(100, 500, 50)
cc.calculate_capture_coefficient(volume=1e-21, temperature=temperatures)

# 4. Results
idx_300K = np.argmin(np.abs(temperatures - 300))
print(f"Capture coefficient at 300K: {cc.capture_coefficient[idx_300K]:.3e} cm³/s")

# 5. Visualize
fig = plot_capture_coefficient(cc)
fig.show()
```

## Troubleshooting

### "Partition function not converged"

Increase `nev` for the initial state:

```python
pot_initial.solve(nev=200)  # Was 180
```

### Very slow calculation

Reduce grid resolution:

```python
pot_initial = Potential.from_harmonic(
    hw=0.008, Q0=0.0, E0=0.5,
    npoints=2000  # Was 5000
)
```

### ARPACK convergence failure

Increase maxiter or adjust tolerance:

```python
pot_initial.solve(nev=180, maxiter=100000)
```

## What's Next?

### Try Different Parameters

Explore how results change:

```python
# Vary displacement
pot_final = Potential.from_harmonic(hw=0.008, Q0=15.0, E0=0.0)  # Was 10.5

# Vary coupling
cc = ConfigCoordinate(pot_i=pot_initial, pot_f=pot_final, W=0.1)  # Was 0.068

# Vary phonon energy
pot = Potential.from_harmonic(hw=0.016, Q0=0.0, E0=0.5)  # Was 0.008
```

### Real Systems

Move beyond harmonic:

- **[User Guide: Potentials](../user-guide/potentials.md)** - Spline and Morse fitting
- **[doped Integration](../user-guide/doped-integration.md)** - Load from defect calculations

### Parameter Scanning

Screen many materials:

- **[Example Notebook: Parameter Scan](https://github.com/WMD-group/CarrierCapture.py/blob/main/examples/notebooks/03_parameter_scan.ipynb)** - High-throughput workflow
- **[User Guide: Parameter Scanning](../user-guide/parameter-scanning.md)** - Detailed guide

### Interactive Exploration

- **[User Guide: Visualization](../user-guide/visualization.md)** - Dashboard and plotting

## Summary

You've learned:

✓ Create harmonic potentials
✓ Solve Schrödinger equation
✓ Calculate overlap integrals
✓ Compute capture coefficients
✓ Interpret temperature dependence
✓ Visualize and export results

You're now ready to apply CarrierCapture to real systems!
