# CLI Reference

Command-line interface for CarrierCapture.

## Overview

CarrierCapture provides a comprehensive CLI built with [Click](https://click.palletsprojects.com/) for all common workflows:

- **`fit`** - Fit potential energy surfaces to data
- **`solve`** - Solve Schrödinger equation for phonon states
- **`capture`** - Calculate carrier capture coefficients
- **`scan`** - High-throughput parameter screening
- **`scan-plot`** - Visualize scan results
- **`viz`** - Launch interactive dashboard
- **`plot`** - Generate static plots

## Installation

The CLI is automatically available after installation:

```bash
pip install carriercapture
carriercapture --help
```

## Global Options

Available for all commands:

| Option | Description |
|--------|-------------|
| `--version` | Show version and exit |
| `-v, --verbose` | Increase verbosity (can repeat: `-v`, `-vv`, `-vvv`) |
| `--help` | Show help message and exit |

**Example:**
```bash
# Show version
carriercapture --version

# Run with verbose output
carriercapture fit data.dat -v

# Maximum verbosity
carriercapture capture config.yaml -vvv
```

---

## Commands

### `fit` - Fit Potential Energy Surface

Fit potential energy surface data to various functional forms.

**Usage:**
```bash
carriercapture fit DATA_FILE [OPTIONS]
```

**Arguments:**
- `DATA_FILE` - Path to potential data file (Q, E format)

**Options:**

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-f, --fit-type` | choice | `spline` | Fitting method: `spline`, `harmonic`, `morse`, `polynomial`, `morse_poly` |
| `-o, --output` | path | - | Output file (auto-detect format from extension) |
| `--order` | int | 4 | Spline order (for spline fitting) |
| `--smoothness` | float | 0.001 | Spline smoothness parameter |
| `--degree` | int | 4 | Polynomial degree (for polynomial fitting) |
| `--hw` | float | - | Phonon frequency in eV (for harmonic fitting) |
| `--Q0` | float | 0.0 | Equilibrium position (amu^0.5·Å) |
| `--E0` | float | 0.0 | Energy offset (eV) |
| `--plot` | flag | False | Generate plot of fitted potential |
| `--plot-output` | path | - | Save plot to file (requires `--plot`) |

**Examples:**

```bash
# Spline fit with custom parameters
carriercapture fit excited.dat -f spline --order 4 --smoothness 0.001 -o excited_fit.json

# Harmonic fit
carriercapture fit ground.dat -f harmonic --hw 0.03 --Q0 5.0 -o ground_fit.json

# Morse potential fit
carriercapture fit data.dat -f morse -o morse_fit.json

# Polynomial fit with plot
carriercapture fit data.dat -f polynomial --degree 6 --plot --plot-output fit.png

# Fit and save
carriercapture fit data.dat -f spline -o potential.json -v
```

**Supported Input Formats:**
- CSV (`.csv`, `.dat`) - Two-column format: Q (amu^0.5·Å), E (eV)
- Space or comma separated

**Supported Output Formats:**
- JSON (`.json`)
- YAML (`.yaml`, `.yml`)
- NPZ (`.npz`)
- DAT (`.dat`)

---

### `solve` - Solve Schrödinger Equation

Solve 1D Schrödinger equation to compute phonon eigenvalues and eigenvectors.

**Usage:**
```bash
carriercapture solve POTENTIAL_FILE [OPTIONS]
```

**Arguments:**
- `POTENTIAL_FILE` - Path to fitted potential (JSON, YAML, NPZ)

**Options:**

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-n, --nev` | int | 60 | Number of eigenvalues to compute |
| `-o, --output` | path | `*_solved.*` | Output file (default: add `_solved` suffix) |
| `--Q-range` | float×2 | - | Q grid range: `Q_min Q_max` (amu^0.5·Å) |
| `--npoints` | int | 3001 | Number of grid points for solving |
| `--plot` | flag | False | Plot eigenvalue spectrum |
| `--plot-output` | path | - | Save plot to file |
| `--plot-wavefunctions` | int | - | Number of wavefunctions to plot |

**Examples:**

```bash
# Solve with default parameters (60 states)
carriercapture solve potential.json

# Solve for more states with custom grid
carriercapture solve potential.json -n 180 --npoints 5000 -o solved.json

# Solve and plot eigenvalue spectrum
carriercapture solve potential.json --plot --plot-output spectrum.png

# Solve and visualize wavefunctions
carriercapture solve potential.json --plot-wavefunctions 10 --show

# Solve with custom Q range
carriercapture solve potential.json --Q-range -10 10 -n 180 -v
```

**Output:**
- Saves potential with computed `eigenvalues` and `eigenvectors` arrays
- Eigenvalues in ascending order (E₀, E₁, E₂, ...)
- Eigenvectors normalized on the Q grid

---

### `capture` - Calculate Capture Coefficient

Calculate carrier capture coefficients using multiphonon theory.

**Usage:**
```bash
carriercapture capture [CONFIG_FILE] [OPTIONS]
```

**Arguments:**
- `CONFIG_FILE` (optional) - YAML configuration file with all parameters

**Modes:**

1. **Config file mode**: Provide YAML file with all parameters
2. **Command-line mode**: Specify all options via CLI
3. **doped mode**: Integrate with doped package (`--doped` flag)

**Options:**

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--pot-i` | path | - | Initial state potential file |
| `--pot-f` | path | - | Final state potential file |
| `-W, --coupling` | float | - | Electron-phonon coupling (eV) |
| `-g, --degeneracy` | int | 1 | Degeneracy factor |
| `-V, --volume` | float | - | Supercell volume (cm³) |
| `--temp-range` | float×3 | 100 500 50 | Temperature range: `T_min T_max n_points` (K) |
| `--Q0` | float | - | Shift for coordinate operator (amu^0.5·Å) |
| `--cutoff` | float | 0.25 | Energy cutoff for overlaps (eV) |
| `--sigma` | float | 0.025 | Gaussian delta width (eV) |
| `-o, --output` | path | - | Output file (`.json`, `.yaml`, `.csv`, `.npz`) |
| `--plot` | flag | False | Generate Arrhenius plot |
| `--plot-output` | path | - | Save plot to file |

**doped Integration Options:**

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--doped` | path | - | Path to doped DefectEntry JSON.GZ file |
| `--charge-i` | int | - | Initial charge state |
| `--charge-f` | int | - | Final charge state |
| `--doped-path-i` | path | - | VASP calculation directory (initial state) |
| `--doped-path-f` | path | - | VASP calculation directory (final state) |
| `--n-images` | int | 10 | Number of interpolation images |
| `--auto-Q0` | flag | False | Auto-suggest Q0 from structure displacement |

**Examples:**

```bash
# From config file
carriercapture capture config.yaml -o results.json -v

# From command line
carriercapture capture \
  --pot-i excited.json --pot-f ground.json \
  -W 0.205 -V 1e-21 --temp-range 100 500 50 \
  --Q0 10.0 -o results.json

# With Arrhenius plot
carriercapture capture config.yaml \
  --plot --plot-output arrhenius.png -v

# doped integration mode
carriercapture capture --doped defect.json.gz \
  --charge-i 0 --charge-f +1 \
  --doped-path-i path_q0/ --doped-path-f path_q1/ \
  -W 0.205 -V 1e-21 --temp-range 100 500 50 \
  --auto-Q0 -o results.json -v

# Quick calculation with plotting
carriercapture capture config.yaml \
  --temp-range 200 400 20 \
  --plot --show -vv
```

**Config File Format (YAML):**

```yaml
potential_initial:
  file: excited_solved.json

potential_final:
  file: ground_solved.json

capture:
  W: 0.205              # eV
  degeneracy: 1
  volume: 1.0e-21       # cm³
  Q0: 10.0              # amu^0.5·Å
  cutoff: 0.25          # eV
  sigma: 0.025          # eV
  temperature:
    min: 100            # K
    max: 500            # K
    n_points: 50
```

---

### `scan` - Parameter Scan

High-throughput parameter screening over ΔQ and ΔE space.

**Usage:**
```bash
carriercapture scan [OPTIONS]
```

**Required Options:**

| Option | Type | Description |
|--------|------|-------------|
| `--dQ-min` | float | Minimum ΔQ value (amu^0.5·Å) |
| `--dQ-max` | float | Maximum ΔQ value (amu^0.5·Å) |
| `--dQ-points` | int | Number of ΔQ points |
| `--dE-min` | float | Minimum ΔE value (eV) |
| `--dE-max` | float | Maximum ΔE value (eV) |
| `--dE-points` | int | Number of ΔE points |
| `-o, --output` | path | Output file (`.npz` or `.h5`) |

**Optional Parameters:**

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--hbar-omega-i` | float | 0.008 | ℏω for initial state (eV) |
| `--hbar-omega-f` | float | 0.008 | ℏω for final state (eV) |
| `-T, --temperature` | float | 300.0 | Temperature (K) |
| `-V, --volume` | float | 1e-21 | Supercell volume (cm³) |
| `-g, --degeneracy` | int | 1 | Degeneracy factor |
| `--sigma` | float | 0.01 | Gaussian delta width (eV) |
| `--cutoff` | float | 0.25 | Energy cutoff for overlaps (eV) |
| `--nev-i` | int | 180 | Number of eigenvalues (initial state) |
| `--nev-f` | int | 60 | Number of eigenvalues (final state) |
| `-j, --n-jobs` | int | 1 | Parallel jobs (use `-1` for all cores) |
| `--no-progress` | flag | False | Disable progress bar |

**Examples:**

```bash
# Basic scan over ΔQ and ΔE
carriercapture scan \
  --dQ-min 0 --dQ-max 25 --dQ-points 25 \
  --dE-min 0 --dE-max 2.5 --dE-points 10 \
  -o scan_results.npz

# Parallel scan with 4 cores
carriercapture scan \
  --dQ-min 0 --dQ-max 25 --dQ-points 25 \
  --dE-min 0 --dE-max 2.5 --dE-points 10 \
  -j 4 -o results.npz -v

# Use all available cores
carriercapture scan \
  --dQ-min 0 --dQ-max 25 --dQ-points 50 \
  --dE-min 0 --dE-max 2.5 --dE-points 20 \
  -j -1 -o results.npz

# Custom phonon frequencies and temperature
carriercapture scan \
  --dQ-min 0 --dQ-max 25 --dQ-points 25 \
  --dE-min 0 --dE-max 2.5 --dE-points 10 \
  --hbar-omega-i 0.010 --hbar-omega-f 0.010 \
  -T 500 -o results_500K.npz -vv

# High-resolution scan
carriercapture scan \
  --dQ-min 0 --dQ-max 30 --dQ-points 100 \
  --dE-min 0 --dE-max 3.0 --dE-points 50 \
  -j -1 --no-progress -o high_res_scan.npz
```

**Output:**
- Saves `ScanResult` object with:
  - `dQ_grid`, `dE_grid` - Parameter grids
  - `capture_coefficients` - 2D array of C(T) values
  - Scan parameters and metadata

---

### `scan-plot` - Visualize Scan Results

Generate heatmaps or contour plots from parameter scan results.

**Usage:**
```bash
carriercapture scan-plot SCAN_FILE [OPTIONS]
```

**Arguments:**
- `SCAN_FILE` - Path to scan results file (`.npz` or `.h5`)

**Options:**

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--type` | choice | `heatmap` | Plot type: `heatmap`, `contour`, `both` |
| `--log-scale` | flag | False | Use log₁₀ scale for capture coefficients |
| `-o, --output` | path | - | Output file (`.html`, `.png`) |
| `--show` | flag | False | Display plot in browser |

**Examples:**

```bash
# Plot heatmap
carriercapture scan-plot results.npz --show

# Save with log scale
carriercapture scan-plot results.npz --log-scale -o heatmap.html

# Generate contour plot
carriercapture scan-plot results.npz --type contour --show

# Save both types
carriercapture scan-plot results.npz --type both --log-scale -o scan_viz.html

# Quick view
carriercapture scan-plot results.npz --show -v
```

---

### `viz` - Interactive Dashboard

Launch web-based interactive visualization dashboard.

**Usage:**
```bash
carriercapture viz [OPTIONS]
```

**Options:**

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--port` | int | 8050 | Port for Dash server |
| `--host` | str | 127.0.0.1 | Host address |
| `--debug` | flag | False | Run in debug mode with hot reloading |
| `--no-browser` | flag | False | Don't automatically open browser |
| `--data` | path | - | Load data file on startup |

**Examples:**

```bash
# Launch dashboard on default port (8050)
carriercapture viz

# Launch on custom port
carriercapture viz --port 8080

# Launch with data file preloaded
carriercapture viz --data potential.json

# Run in debug mode
carriercapture viz --debug

# Run on network (accessible from other machines)
carriercapture viz --host 0.0.0.0 --port 8050 -v

# Launch without opening browser
carriercapture viz --no-browser
```

**Dashboard Features:**
- Interactive potential energy surface plotting
- Real-time fitting with parameter adjustment
- Schrödinger equation solver with visualization
- Configuration coordinate diagram explorer
- Capture coefficient calculation
- Export results and figures

**Requirements:**
```bash
pip install carriercapture[viz]
```

---

### `plot` - Static Plots

Generate static publication-quality plots from potential data.

**Usage:**
```bash
carriercapture plot POTENTIAL_FILE [OPTIONS]
```

**Arguments:**
- `POTENTIAL_FILE` - Path to potential file (JSON, YAML, NPZ)

**Options:**

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--type` | choice | `potential` | Plot type: `potential`, `spectrum`, `both` |
| `--show-wf` | flag | False | Show wavefunctions on potential plot |
| `--max-wf` | int | 20 | Maximum number of wavefunctions to plot |
| `-o, --output` | path | - | Output file (`.html`, `.png`, `.pdf`, `.svg`) |
| `--width` | int | 900 | Figure width in pixels |
| `--height` | int | 600 | Figure height in pixels |
| `--show` | flag | False | Display plot in browser |

**Examples:**

```bash
# Plot potential energy surface
carriercapture plot potential.json --show

# Plot with wavefunctions
carriercapture plot potential.json --show-wf --max-wf 10 -o figure.html

# Plot eigenvalue spectrum
carriercapture plot potential.json --type spectrum -o spectrum.png

# Generate both plots
carriercapture plot potential.json --type both --show-wf --show

# High-resolution export
carriercapture plot potential.json \
  --show-wf --width 1200 --height 800 \
  -o figure.pdf

# Quick visualization
carriercapture plot potential.json --show-wf --show -v
```

---

## Common Workflows

### 1. Complete Workflow (Fit → Solve → Capture)

```bash
# Step 1: Fit potentials
carriercapture fit excited.dat -f spline -o excited.json -v
carriercapture fit ground.dat -f spline -o ground.json -v

# Step 2: Solve Schrödinger equation
carriercapture solve excited.json -n 180 -o excited_solved.json -v
carriercapture solve ground.json -n 60 -o ground_solved.json -v

# Step 3: Calculate capture coefficient
carriercapture capture \
  --pot-i excited_solved.json --pot-f ground_solved.json \
  -W 0.205 -V 1e-21 --Q0 10.0 \
  --temp-range 100 500 50 \
  -o capture_results.json --plot --plot-output arrhenius.png -vv
```

### 2. Config File Workflow

Create `config.yaml`:
```yaml
potential_initial:
  file: excited_solved.json
potential_final:
  file: ground_solved.json
capture:
  W: 0.205
  volume: 1.0e-21
  Q0: 10.0
  temperature:
    min: 100
    max: 500
    n_points: 50
```

Run:
```bash
carriercapture capture config.yaml -o results.json --plot -v
```

### 3. Parameter Screening Workflow

```bash
# Run scan
carriercapture scan \
  --dQ-min 0 --dQ-max 25 --dQ-points 25 \
  --dE-min 0 --dE-max 2.5 --dE-points 10 \
  -j -1 -o scan.npz -v

# Visualize results
carriercapture scan-plot scan.npz --log-scale --show -v
```

### 4. Interactive Exploration

```bash
# Launch dashboard for interactive exploration
carriercapture viz --port 8050
```

### 5. doped Integration Workflow

```bash
# Calculate capture from doped defect data
carriercapture capture --doped defect.json.gz \
  --charge-i 0 --charge-f +1 \
  --doped-path-i vasp_q0/ --doped-path-f vasp_q1/ \
  -W 0.205 -V 1e-21 --auto-Q0 \
  --temp-range 100 500 50 \
  -o doped_capture.json --plot -vv
```

---

## Tips and Tricks

### Verbosity Levels

```bash
# No verbosity (minimal output)
carriercapture fit data.dat -o fit.json

# Level 1 (-v): Progress and status
carriercapture fit data.dat -o fit.json -v

# Level 2 (-vv): Detailed information
carriercapture fit data.dat -o fit.json -vv

# Level 3 (-vvv): Debug output
carriercapture fit data.dat -o fit.json -vvv
```

### Parallel Processing

```bash
# Use specific number of cores
carriercapture scan ... -j 4

# Use all available cores
carriercapture scan ... -j -1

# Single-threaded (for debugging)
carriercapture scan ... -j 1
```

### File Format Auto-Detection

Output format is automatically detected from file extension:

```bash
carriercapture fit data.dat -o fit.json    # JSON
carriercapture fit data.dat -o fit.yaml    # YAML
carriercapture fit data.dat -o fit.npz     # NumPy compressed
carriercapture fit data.dat -o fit.dat     # Plain text
```

### Quick Visualization

```bash
# Fit and immediately visualize
carriercapture fit data.dat -f spline --plot --show

# Solve and plot wavefunctions
carriercapture solve potential.json --plot-wavefunctions 10 --show

# Calculate and plot capture
carriercapture capture config.yaml --plot --show
```

---

## Environment Variables

CarrierCapture respects these environment variables:

| Variable | Description |
|----------|-------------|
| `CARRIERCAPTURE_DATA_DIR` | Default data directory |
| `CARRIERCAPTURE_CACHE_DIR` | Cache directory for results |
| `OMP_NUM_THREADS` | OpenMP threads for BLAS/LAPACK |
| `MKL_NUM_THREADS` | MKL threads (if using Intel MKL) |

---

## Error Handling

### Common Errors

**1. Potential not fitted**
```bash
$ carriercapture solve potential.json
Error: Potential must be fitted before solving
```
→ Run `fit` command first

**2. Potential not solved**
```bash
$ carriercapture capture --pot-i excited.json --pot-f ground.json ...
Error: Initial potential must be solved (use 'solve' command)
```
→ Run `solve` command first

**3. Missing required parameters**
```bash
$ carriercapture capture --pot-i excited.json
Error: Must provide -W/--coupling (or config file)
```
→ Provide all required parameters

**4. Import errors (optional dependencies)**
```bash
$ carriercapture viz
Error: Dash dependencies not available. Install with: pip install carriercapture[viz]
```
→ Install optional dependencies

---

## See Also

- **[User Guide: CLI Usage](../user-guide/cli-usage.md)** - Comprehensive CLI guide with examples
- **[Getting Started: Quick Start](../getting-started/quick-start.md)** - 5-minute introduction
- **[API Reference: Core](core.md)** - Python API reference
