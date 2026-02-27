# Parameter Scanning

High-throughput screening of carrier capture rates across parameter space.

## Overview

The `ParameterScanner` class enables systematic exploration of capture coefficients as a function of key materials parameters:

- **ΔQ** - Configuration coordinate offset (structural displacement)
- **ΔE** - Energy offset between charge states
- **ℏω** - Phonon frequencies

This is useful for:
- **Materials screening** - Identify promising defect-material combinations
- **Sensitivity analysis** - Understand parameter dependencies
- **Design principles** - Find trends for defect engineering
- **Benchmark comparisons** - Test against databases

---

## Quick Start

```python
from carriercapture.analysis import ParameterScanner, ScanParameters
import numpy as np

# Define scan parameters
params = ScanParameters(
    dQ_range=(0, 25, 25),      # ΔQ: 0-25 amu^0.5·Å, 25 points
    dE_range=(0, 2.5, 10),     # ΔE: 0-2.5 eV, 10 points
    hbar_omega_i=0.008,        # Initial state phonon (eV)
    hbar_omega_f=0.008,        # Final state phonon (eV)
    temperature=300.0,         # Temperature (K)
    volume=1e-21,              # Supercell volume (cm³)
    degeneracy=1,
    nev_initial=180,
    nev_final=60,
)

# Create scanner
scanner = ParameterScanner(params, verbose=True)

# Run scan (parallel execution)
results = scanner.run_harmonic_scan(n_jobs=-1, show_progress=True)

# Analyze results
print(f"Scanned {results.capture_coefficients.size} parameter combinations")
print(f"Max C: {np.nanmax(results.capture_coefficients):.3e} cm³/s")
print(f"Min C: {np.nanmin(results.capture_coefficients):.3e} cm³/s")

# Save results
results.save('scan_results.npz')
```

---

## Parameter Space

### ScanParameters Dataclass

Configure all scan parameters:

```python
from carriercapture.analysis import ScanParameters

params = ScanParameters(
    dQ_range=(Q_min, Q_max, n_points),  # ΔQ scan range
    dE_range=(E_min, E_max, n_points),  # ΔE scan range
    hbar_omega_i=0.008,                 # Initial phonon (eV)
    hbar_omega_f=0.008,                 # Final phonon (eV)
    temperature=300.0,                  # Temperature (K)
    volume=1e-21,                       # Volume (cm³)
    degeneracy=1,                       # Degeneracy
    sigma=0.01,                         # Gaussian width (eV)
    cutoff=0.25,                        # Energy cutoff (eV)
    nev_initial=180,                    # Initial eigenvalues
    nev_final=60,                       # Final eigenvalues
)
```

**Key Parameters:**

| Parameter | Type | Description | Typical Range |
|-----------|------|-------------|---------------|
| `dQ_range` | tuple | (min, max, points) for ΔQ | (0, 30, 25-50) |
| `dE_range` | tuple | (min, max, points) for ΔE | (0, 3, 10-20) |
| `hbar_omega_i` | float | Initial phonon energy (eV) | 0.005-0.015 |
| `hbar_omega_f` | float | Final phonon energy (eV) | 0.005-0.015 |
| `temperature` | float | Temperature (K) | 100-500 |
| `volume` | float | Supercell volume (cm³) | 1e-21 - 1e-20 |
| `nev_initial` | int | Initial state eigenvalues | 100-200 |
| `nev_final` | int | Final state eigenvalues | 40-80 |

### Choosing Scan Ranges

**ΔQ Range:**

Configuration coordinate offset representing structural distortion.

```python
# Typical values for semiconductors:
# Small displacement: ΔQ ~ 5-10 amu^0.5·Å
# Medium displacement: ΔQ ~ 10-20 amu^0.5·Å
# Large displacement: ΔQ ~ 20-30 amu^0.5·Å

# Example: ZnO defects
dQ_range = (0, 25, 25)  # 0-25 amu^0.5·Å, 25 points

# Example: High-resolution near interesting region
dQ_range = (8, 15, 50)  # Dense sampling 8-15 amu^0.5·Å

# Example: Wide exploratory scan
dQ_range = (0, 40, 30)  # Broader range
```

**ΔE Range:**

Energy difference between charge states.

```python
# Typical values:
# Shallow defects: ΔE ~ 0.1-0.5 eV
# Mid-gap defects: ΔE ~ 0.5-1.5 eV
# Deep defects: ΔE ~ 1.5-3.0 eV

# Example: Screen full band gap
dE_range = (0, 2.5, 10)  # 0-2.5 eV, 10 points

# Example: Focus on mid-gap region
dE_range = (0.5, 2.0, 20)  # 0.5-2.0 eV, 20 points

# Example: Fine scan near threshold
dE_range = (0, 0.5, 25)  # 0-0.5 eV, high resolution
```

**Phonon Energies:**

Characteristic phonon frequencies for each charge state.

```python
# Typical values for semiconductors:
# Acoustic phonons: 0.002-0.005 eV (2-5 meV)
# Optical phonons: 0.005-0.015 eV (5-15 meV)
# High-frequency modes: 0.015-0.030 eV (15-30 meV)

# Example: ZnO (optical phonons ~8 meV)
hbar_omega_i = 0.008  # eV
hbar_omega_f = 0.008  # eV

# Example: Different phonons for each state
hbar_omega_i = 0.010  # Initial state
hbar_omega_f = 0.008  # Final state

# From DFT phonon calculations:
from carriercapture.analysis import estimate_phonon_energy
# hbar_omega = estimate_phonon_energy(phonopy_yaml_file)
```

### Grid Resolution

Trade-off between accuracy and computational cost:

```python
# Coarse scan (fast, ~250 points)
params_coarse = ScanParameters(
    dQ_range=(0, 25, 10),   # 10 ΔQ points
    dE_range=(0, 2.5, 5),   # 5 ΔE points → 50 total
    # ... other params
)

# Medium scan (standard, ~500 points)
params_medium = ScanParameters(
    dQ_range=(0, 25, 25),   # 25 ΔQ points
    dE_range=(0, 2.5, 10),  # 10 ΔE points → 250 total
)

# Fine scan (high-res, ~2000 points)
params_fine = ScanParameters(
    dQ_range=(0, 25, 50),   # 50 ΔQ points
    dE_range=(0, 2.5, 20),  # 20 ΔE points → 1000 total
)

# Estimate computation time
n_points = 25 * 10  # 250
time_per_point = 0.5  # seconds (typical)
total_time = n_points * time_per_point / 60  # minutes
print(f"Estimated time: {total_time:.1f} minutes (serial)")
print(f"With 8 cores: {total_time/8:.1f} minutes")
```

---

## Running Scans

### Basic Harmonic Scan

Use harmonic potentials for fast screening:

```python
from carriercapture.analysis import ParameterScanner, ScanParameters

# Configure scan
params = ScanParameters(
    dQ_range=(0, 25, 25),
    dE_range=(0, 2.5, 10),
    hbar_omega_i=0.008,
    hbar_omega_f=0.008,
    temperature=300.0,
    volume=1e-21,
)

# Create scanner
scanner = ParameterScanner(params, verbose=True)

# Run scan (serial)
results = scanner.run_harmonic_scan(n_jobs=1, show_progress=True)

print(f"✓ Scan complete: {results.capture_coefficients.shape}")
```

**Why harmonic?**
- **Fast**: Analytical potential, no fitting needed
- **Sufficient**: Captures main trends for screening
- **Scalable**: Enables large parameter sweeps

### Parallel Execution

Leverage multiple CPU cores:

```python
# Use all available cores
results = scanner.run_harmonic_scan(n_jobs=-1, show_progress=True)

# Use specific number of cores
results = scanner.run_harmonic_scan(n_jobs=4, show_progress=True)

# Serial (for debugging)
results = scanner.run_harmonic_scan(n_jobs=1, show_progress=True)
```

**Speedup:**
```python
import multiprocessing
n_cores = multiprocessing.cpu_count()
print(f"Available cores: {n_cores}")

# Typical speedup: 0.7-0.9 × n_cores (due to overhead)
# 8 cores → ~6× speedup
# 16 cores → ~12× speedup
```

### Progress Monitoring

```python
# With progress bar (default)
results = scanner.run_harmonic_scan(show_progress=True)

# Without progress bar (for scripts/logs)
results = scanner.run_harmonic_scan(show_progress=False)

# Verbose output
scanner = ParameterScanner(params, verbose=True)
results = scanner.run_harmonic_scan(n_jobs=-1, show_progress=True)
# Prints:
# [ParameterScanner] Scanning 25×10 parameter grid...
# [ParameterScanner] Using 8 parallel workers
# 100%|██████████| 250/250 [02:15<00:00, 1.85it/s]
```

---

## Analyzing Results

### ScanResult Object

Results are stored in a `ScanResult` object:

```python
# Access scan data
dQ_grid = results.dQ_grid              # 1D array: ΔQ values
dE_grid = results.dE_grid              # 1D array: ΔE values
C_matrix = results.capture_coefficients # 2D array: C(ΔQ, ΔE)

print(f"ΔQ grid: {dQ_grid.shape}")      # (25,)
print(f"ΔE grid: {dE_grid.shape}")      # (10,)
print(f"C matrix: {C_matrix.shape}")    # (25, 10)

# Access parameters
print(f"Temperature: {results.temperature} K")
print(f"Phonon ℏω_i: {results.hbar_omega_i} eV")
print(f"Volume: {results.volume:.2e} cm³")
```

### Statistical Analysis

```python
import numpy as np

# Basic statistics
C = results.capture_coefficients
valid_mask = ~np.isnan(C)

print(f"Valid calculations: {np.sum(valid_mask)} / {C.size}")
print(f"Max C: {np.nanmax(C):.3e} cm³/s")
print(f"Min C: {np.nanmin(C):.3e} cm³/s")
print(f"Mean C: {np.nanmean(C):.3e} cm³/s")
print(f"Median C: {np.nanmedian(C):.3e} cm³/s")

# Find maximum
i_max, j_max = np.unravel_index(np.nanargmax(C), C.shape)
dQ_max = results.dQ_grid[i_max]
dE_max = results.dE_grid[j_max]
C_max = C[i_max, j_max]

print(f"\nMaximum capture:")
print(f"  ΔQ = {dQ_max:.2f} amu^0.5·Å")
print(f"  ΔE = {dE_max:.2f} eV")
print(f"  C = {C_max:.3e} cm³/s")

# Find minimum (non-zero)
C_nonzero = C[valid_mask & (C > 0)]
C_min_nonzero = np.min(C_nonzero)
print(f"Minimum (non-zero) C: {C_min_nonzero:.3e} cm³/s")
```

### Slicing Results

```python
# Fix ΔE, vary ΔQ
dE_fixed = 1.0  # eV
j = np.argmin(np.abs(results.dE_grid - dE_fixed))
C_vs_dQ = results.capture_coefficients[:, j]

import matplotlib.pyplot as plt
plt.figure()
plt.semilogy(results.dQ_grid, C_vs_dQ, 'o-')
plt.xlabel('ΔQ (amu$^{0.5}$·Å)')
plt.ylabel('C (cm$^3$/s)')
plt.title(f'Capture vs ΔQ at ΔE = {results.dE_grid[j]:.2f} eV')
plt.grid(True, alpha=0.3)
plt.show()

# Fix ΔQ, vary ΔE
dQ_fixed = 15.0  # amu^0.5·Å
i = np.argmin(np.abs(results.dQ_grid - dQ_fixed))
C_vs_dE = results.capture_coefficients[i, :]

plt.figure()
plt.semilogy(results.dE_grid, C_vs_dE, 's-')
plt.xlabel('ΔE (eV)')
plt.ylabel('C (cm$^3$/s)')
plt.title(f'Capture vs ΔE at ΔQ = {results.dQ_grid[i]:.2f} amu$^{{0.5}}$·Å')
plt.grid(True, alpha=0.3)
plt.show()
```

### Identifying Trends

```python
# Where is capture fastest?
C = results.capture_coefficients
threshold = 1e-10  # cm³/s

fast_capture_mask = C > threshold
n_fast = np.sum(fast_capture_mask)
print(f"Fast capture (C > {threshold:.0e}): {n_fast} / {C.size} points")

# Get (ΔQ, ΔE) coordinates for fast capture
dQ_fast = results.dQ_grid[np.any(fast_capture_mask, axis=1)]
dE_fast = results.dE_grid[np.any(fast_capture_mask, axis=0)]
print(f"ΔQ range for fast capture: {dQ_fast.min():.1f} - {dQ_fast.max():.1f} amu^0.5·Å")
print(f"ΔE range for fast capture: {dE_fast.min():.2f} - {dE_fast.max():.2f} eV")

# Identify "sweet spots"
i_indices, j_indices = np.where(fast_capture_mask)
for i, j in zip(i_indices[:5], j_indices[:5]):  # Print first 5
    dQ = results.dQ_grid[i]
    dE = results.dE_grid[j]
    C_val = C[i, j]
    print(f"  ΔQ={dQ:5.1f}, ΔE={dE:4.2f} → C={C_val:.2e} cm³/s")
```

---

## Visualization

### 2D Heatmap

```python
from carriercapture.visualization import plot_scan_heatmap

# Basic heatmap
fig = plot_scan_heatmap(results, log_scale=True)
fig.show()

# Customized
fig = plot_scan_heatmap(
    results,
    title=f"Capture Coefficient at {results.temperature}K",
    log_scale=True,
    colorscale='Viridis',
    width=900,
    height=700
)
fig.write_html('scan_heatmap.html')
```

### Custom Plotting with Plotly

```python
import plotly.graph_objects as go

# Prepare data
dQ_2d, dE_2d = np.meshgrid(results.dQ_grid, results.dE_grid, indexing='ij')
C_log = np.log10(results.capture_coefficients + 1e-30)

# Create heatmap
fig = go.Figure(data=go.Heatmap(
    x=results.dE_grid,
    y=results.dQ_grid,
    z=C_log,
    colorscale='Viridis',
    colorbar=dict(title='log₁₀(C) [cm³/s]')
))

fig.update_layout(
    title='Parameter Scan: Capture Coefficient',
    xaxis_title='ΔE (eV)',
    yaxis_title='ΔQ (amu<sup>0.5</sup>·Å)',
    width=800,
    height=700
)

fig.show()
```

### Contour Plot

```python
import plotly.graph_objects as go

# Contour levels
levels = [-16, -14, -12, -10, -8, -6]  # log₁₀(C)

fig = go.Figure(data=go.Contour(
    x=results.dE_grid,
    y=results.dQ_grid,
    z=C_log,
    colorscale='Viridis',
    contours=dict(
        start=levels[0],
        end=levels[-1],
        size=2,
        showlabels=True,
    ),
    colorbar=dict(title='log₁₀(C) [cm³/s]')
))

fig.update_layout(
    title='Capture Coefficient Contours',
    xaxis_title='ΔE (eV)',
    yaxis_title='ΔQ (amu<sup>0.5</sup>·Å)',
)

fig.show()
```

---

## Saving and Loading

### Save Results

```python
# Save to NPZ (NumPy compressed, fast)
results.save('scan_results.npz', format='npz')

# Save to HDF5 (for very large datasets)
results.save('scan_results.h5', format='hdf5')

# File sizes (typical 25×10 grid):
# NPZ: ~50 KB
# HDF5: ~30 KB
```

### Load Results

```python
from carriercapture.analysis import ScanResult

# Load from NPZ
results = ScanResult.load('scan_results.npz', format='npz')

# Load from HDF5
results = ScanResult.load('scan_results.h5', format='hdf5')

# Access data immediately
print(f"Loaded scan: {results.capture_coefficients.shape}")
print(f"Temperature: {results.temperature} K")
```

### Combining Multiple Scans

```python
# Scan at different temperatures
temperatures = [200, 300, 400, 500]
all_results = {}

for T in temperatures:
    params = ScanParameters(
        dQ_range=(0, 25, 25),
        dE_range=(0, 2.5, 10),
        temperature=T,
        # ... other params
    )
    scanner = ParameterScanner(params, verbose=True)
    results = scanner.run_harmonic_scan(n_jobs=-1)
    all_results[T] = results
    results.save(f'scan_T{T}K.npz')

# Compare captures at different T
import matplotlib.pyplot as plt
fig, ax = plt.subplots()

for T, res in all_results.items():
    # Take diagonal slice
    n = min(len(res.dQ_grid), len(res.dE_grid))
    diagonal_C = np.diag(res.capture_coefficients[:n, :n])
    ax.semilogy(range(n), diagonal_C, 'o-', label=f'{T}K')

ax.set_xlabel('Grid index')
ax.set_ylabel('C (cm$^3$/s)')
ax.legend()
ax.grid(True, alpha=0.3)
plt.show()
```

---

## Advanced Usage

### Multi-Temperature Scan

Scan over (ΔQ, ΔE, T) space:

```python
temperatures = np.linspace(100, 500, 5)
results_multi_T = {}

for T in temperatures:
    params = ScanParameters(
        dQ_range=(0, 25, 25),
        dE_range=(0, 2.5, 10),
        temperature=T,
        # ... other params
    )
    scanner = ParameterScanner(params, verbose=False)
    results = scanner.run_harmonic_scan(n_jobs=-1, show_progress=False)
    results_multi_T[T] = results
    print(f"✓ Completed T = {T:.0f}K")

# Create 3D visualization
# (ΔQ, ΔE, T) → C(ΔQ, ΔE, T)
```

### Adaptive Refinement

Refine grid near interesting regions:

```python
# 1. Coarse scan
params_coarse = ScanParameters(
    dQ_range=(0, 30, 10),  # Coarse: 10 points
    dE_range=(0, 3, 6),    # Coarse: 6 points
    # ... other params
)
scanner_coarse = ParameterScanner(params_coarse)
results_coarse = scanner_coarse.run_harmonic_scan(n_jobs=-1)

# 2. Find maximum
C = results_coarse.capture_coefficients
i_max, j_max = np.unravel_index(np.nanargmax(C), C.shape)
dQ_max = results_coarse.dQ_grid[i_max]
dE_max = results_coarse.dE_grid[j_max]

print(f"Coarse maximum at ΔQ={dQ_max:.1f}, ΔE={dE_max:.2f}")

# 3. Fine scan near maximum
dQ_range_fine = (max(0, dQ_max-5), dQ_max+5, 25)
dE_range_fine = (max(0, dE_max-0.5), dE_max+0.5, 15)

params_fine = ScanParameters(
    dQ_range=dQ_range_fine,
    dE_range=dE_range_fine,
    # ... other params
)
scanner_fine = ParameterScanner(params_fine)
results_fine = scanner_fine.run_harmonic_scan(n_jobs=-1)

print(f"Fine scan complete: {results_fine.capture_coefficients.shape}")
```

### Custom Post-Processing

```python
# Calculate additional quantities
C = results.capture_coefficients

# Effective Huang-Rhys parameter
S_eff = (results.dQ_grid[:, None]**2) / (2 * results.hbar_omega_i)
# Shape: (n_dQ, 1) → broadcasts to (n_dQ, n_dE)

# Activation energy (approximate, for large ΔE)
E_a_approx = (results.dE_grid[None, :] - results.dQ_grid[:, None]**2 /
              (8 * results.hbar_omega_i)) / 4
E_a_approx[E_a_approx < 0] = 0

# Plot activation energy map
import plotly.graph_objects as go
fig = go.Figure(data=go.Heatmap(
    x=results.dE_grid,
    y=results.dQ_grid,
    z=E_a_approx,
    colorscale='RdYlBu_r',
    colorbar=dict(title='E_a (eV)')
))
fig.update_layout(
    title='Approximate Activation Energy',
    xaxis_title='ΔE (eV)',
    yaxis_title='ΔQ (amu<sup>0.5</sup>·Å)',
)
fig.show()
```

---

## Command-Line Interface

Run scans from the command line:

```bash
# Basic scan
carriercapture scan \
  --dQ-min 0 --dQ-max 25 --dQ-points 25 \
  --dE-min 0 --dE-max 2.5 --dE-points 10 \
  -o scan_results.npz

# Parallel with 8 cores
carriercapture scan \
  --dQ-min 0 --dQ-max 25 --dQ-points 25 \
  --dE-min 0 --dE-max 2.5 --dE-points 10 \
  -j 8 -o scan_results.npz -v

# Custom parameters
carriercapture scan \
  --dQ-min 0 --dQ-max 25 --dQ-points 25 \
  --dE-min 0 --dE-max 2.5 --dE-points 10 \
  --hbar-omega-i 0.010 --hbar-omega-f 0.008 \
  -T 500 -V 2e-21 -g 2 \
  -j -1 -o scan_high_T.npz

# Visualize results
carriercapture scan-plot scan_results.npz --log-scale --show

# Save plot
carriercapture scan-plot scan_results.npz --log-scale -o heatmap.html
```

---

## Best Practices

### 1. Start with Coarse Scan

```python
# Phase 1: Coarse exploratory scan (fast)
params_coarse = ScanParameters(
    dQ_range=(0, 30, 10),
    dE_range=(0, 3, 6),
    # ... other params
)
results_coarse = ParameterScanner(params_coarse).run_harmonic_scan(n_jobs=-1)

# Phase 2: Identify region of interest
# (see Adaptive Refinement above)

# Phase 3: Fine scan in interesting region
# (targeted high-resolution scan)
```

### 2. Check Convergence

```python
# Test nev convergence on a single point
nev_values = [60, 100, 150, 200]
C_values = []

for nev_i in nev_values:
    params_test = ScanParameters(
        dQ_range=(15, 15, 1),  # Single point
        dE_range=(1.0, 1.0, 1),
        nev_initial=nev_i,
        nev_final=60,
        # ... other params
    )
    results = ParameterScanner(params_test).run_harmonic_scan(n_jobs=1)
    C_values.append(results.capture_coefficients[0, 0])

print("nev_i convergence:")
for nev, C in zip(nev_values, C_values):
    print(f"  nev={nev:3d}: C={C:.3e} cm³/s")

# If converged, use that nev for full scan
```

### 3. Physical Validation

```python
# Check for unphysical results
C = results.capture_coefficients

# Flag suspicious values
too_large = C > 1e-6
too_small = C < 1e-20
suspicious = too_large | too_small

if np.any(suspicious):
    n_suspicious = np.sum(suspicious)
    print(f"⚠️  {n_suspicious} suspicious values detected")

    # Investigate
    i_sus, j_sus = np.where(suspicious)
    for i, j in zip(i_sus[:3], j_sus[:3]):  # Print first 3
        dQ = results.dQ_grid[i]
        dE = results.dE_grid[j]
        C_val = C[i, j]
        print(f"  ΔQ={dQ:.1f}, ΔE={dE:.2f}: C={C_val:.2e}")
```

### 4. Save Intermediate Results

```python
# For long scans, save intermediate results
import time

start_time = time.time()
results = scanner.run_harmonic_scan(n_jobs=-1, show_progress=True)
elapsed = time.time() - start_time

print(f"Scan completed in {elapsed/60:.1f} minutes")

# Save immediately
results.save('scan_results.npz')
print("✓ Results saved")

# Also save metadata
import json
metadata = {
    'date': time.strftime('%Y-%m-%d %H:%M:%S'),
    'elapsed_time_minutes': elapsed / 60,
    'n_points': results.capture_coefficients.size,
    'parameters': {
        'dQ_range': list(params.dQ_range),
        'dE_range': list(params.dE_range),
        'temperature': params.temperature,
        'volume': params.volume,
    }
}
with open('scan_metadata.json', 'w') as f:
    json.dump(metadata, f, indent=2)
```

---

## Complete Example

```python
import numpy as np
from carriercapture.analysis import ParameterScanner, ScanParameters
from carriercapture.visualization import plot_scan_heatmap

# 1. Define scan parameters
params = ScanParameters(
    dQ_range=(0, 25, 25),      # 25 points: 0-25 amu^0.5·Å
    dE_range=(0, 2.5, 10),     # 10 points: 0-2.5 eV
    hbar_omega_i=0.008,        # 8 meV (typical for ZnO)
    hbar_omega_f=0.008,        # 8 meV
    temperature=300.0,         # Room temperature
    volume=1e-21,              # Typical supercell
    degeneracy=1,
    sigma=0.01,                # 10 meV broadening
    cutoff=0.25,               # 250 meV cutoff
    nev_initial=180,           # Neutral state: many phonons
    nev_final=60,              # Charged state: fewer phonons
)

print(f"Scan configuration:")
print(f"  Grid: {params.dQ_range[2]} × {params.dE_range[2]} = {params.dQ_range[2] * params.dE_range[2]} points")
print(f"  Temperature: {params.temperature} K")
print(f"  Phonon energy: {params.hbar_omega_i*1000:.1f} meV")

# 2. Create scanner and run
scanner = ParameterScanner(params, verbose=True)
results = scanner.run_harmonic_scan(n_jobs=-1, show_progress=True)

print(f"✓ Scan complete!")

# 3. Analyze results
C = results.capture_coefficients
valid = ~np.isnan(C)

print(f"\nResults summary:")
print(f"  Valid calculations: {np.sum(valid)} / {C.size}")
print(f"  Max C: {np.nanmax(C):.3e} cm³/s")
print(f"  Min C: {np.nanmin(C):.3e} cm³/s")
print(f"  Mean C: {np.nanmean(C):.3e} cm³/s")

# Find optimal parameters
i_max, j_max = np.unravel_index(np.nanargmax(C), C.shape)
print(f"\nMaximum capture:")
print(f"  ΔQ = {results.dQ_grid[i_max]:.2f} amu^0.5·Å")
print(f"  ΔE = {results.dE_grid[j_max]:.2f} eV")
print(f"  C = {C[i_max, j_max]:.3e} cm³/s")

# 4. Visualize
fig = plot_scan_heatmap(
    results,
    title=f"Carrier Capture Scan at {params.temperature}K",
    log_scale=True,
    colorscale='Viridis'
)
fig.write_html('scan_heatmap.html')
fig.show()

# 5. Save results
results.save('scan_results.npz')
print(f"\n✓ Results saved to scan_results.npz")
```

---

## See Also

- **[API Reference: ParameterScanner](../api/analysis.md)** - Complete API documentation
- **[CLI Reference: scan](../api/cli.md#scan---parameter-scan)** - Command-line usage
- **[Example Notebook: Parameter Scan](https://github.com/WMD-group/CarrierCapture.py/blob/main/examples/notebooks/03_parameter_scan.ipynb)** - Step-by-step example
