# CarrierCapture.py - Agent Development Guide

**Version**: 1.0.0
**Last Updated**: 2026-01-18
**Purpose**: Guide for AI agents and developers working with CarrierCapture.py

---

## Overview

CarrierCapture.py is a Python package for calculating **carrier capture rates** and **non-radiative recombination** in semiconductor defects using **multiphonon theory** (Huang-Rhys formalism). It's a complete rewrite of [CarrierCapture.jl](https://github.com/WMD-group/CarrierCapture.jl) with emphasis on clean code, performance, and interactive visualization.

### Key Capabilities

- **1D Schrödinger equation solver** for vibrational wavefunctions (ARPACK-based)
- **Potential energy surface fitting** (harmonic, Morse, spline)
- **Configuration coordinate diagrams** for defect charge state transitions
- **Capture coefficient calculations** via Fermi's golden rule
- **High-throughput parameter scanning** (parallel with joblib)
- **Interactive visualization** (Dash web dashboard)
- **Integration with doped** for DFT defect calculations

---

## Scientific Background

### The Physics

**Problem**: Calculate the rate at which free carriers (electrons/holes) are captured by defects in semiconductors through non-radiative recombination.

**Method**: Static coupling approximation in multiphonon emission theory

**Key Equation** (Capture Coefficient):

```
C(T) = (V · 2π/ℏ) · g · W² · Σᵢⱼ pᵢ |⟨χᵢ|Q-Q₀|χⱼ⟩|² δ(εᵢ - εⱼ)
```

Where:
- `V`: supercell volume (cm³)
- `g`: degeneracy factor (spin/orbital)
- `W`: electron-phonon coupling matrix element (eV)
- `pᵢ`: Boltzmann occupation of initial vibrational state `i`
- `χᵢ, χⱼ`: vibrational wavefunctions (eigenfunctions of 1D Schrödinger equation)
- `δ(εᵢ - εⱼ)`: energy-conserving delta function (Gaussian-broadened)
- `Q`: configuration coordinate (mass-weighted displacement, amu^0.5·Å)

### Workflow Summary

1. **Define potential energy surfaces**: Two parabolic curves (initial/final charge states) along configuration coordinate Q
2. **Solve Schrödinger equation**: Get vibrational eigenstates (χₙ, εₙ) on each surface
3. **Calculate overlap matrix**: ⟨χᵢ|Q|χⱼ⟩ between initial and final states
4. **Apply Fermi's golden rule**: Sum over thermally occupied states with energy conservation
5. **Result**: Capture coefficient C(T) in cm³/s

### References

- **Alkauskas et al. (2014)**: [Phys. Rev. B 90, 075202](https://doi.org/10.1103/PhysRevB.90.075202) - First-principles calculations
- **Huang & Rhys (1950)**: [Proc. R. Soc. Lond. A 204, 406](https://royalsocietypublishing.org/rspa/article/204/1078/406/8369/Theory-of-light-absorption-and-non-radiative) - Multiphonon theory

---

## Package Architecture

### Directory Structure

```
CarrierCapture.py/
├── src/carriercapture/
│   ├── __init__.py              # Package exports: Potential, ConfigCoordinate
│   ├── __version__.py           # Version string
│   ├── _constants.py            # Physical constants (AMU, HBAR_C, K_B)
│   ├── core/                    # Core physics engine
│   │   ├── potential.py         # Potential class (PES fitting, solving)
│   │   ├── schrodinger.py       # 1D Schrödinger solver (ARPACK)
│   │   ├── config_coord.py      # ConfigCoordinate (capture coefficients)
│   │   └── transfer_coord.py    # TransferCoordinate (experimental)
│   ├── analysis/
│   │   └── parameter_scan.py    # High-throughput scanning (joblib)
│   ├── io/
│   │   ├── readers.py           # Load potential data (JSON, DAT, CSV)
│   │   ├── writers.py           # Save potential data (JSON, NPZ, HDF5)
│   │   └── doped_interface.py   # Integration with doped package
│   ├── visualization/
│   │   ├── static.py            # Plotly plotting functions
│   │   ├── interactive.py       # Dash web dashboard
│   │   └── themes.py            # Color schemes and styles
│   └── cli/                     # Command-line interface
│       ├── main.py              # Click CLI entry point
│       └── commands/            # Subcommands (fit, solve, capture, scan, viz)
├── tests/                       # 88 tests (pytest)
├── examples/
│   ├── notebooks/               # Jupyter tutorials
│   └── data/                    # Example DFT data
├── benchmarks/                  # Validation against CarrierCapture.jl
├── docs/                        # MkDocs documentation
├── pyproject.toml               # Build config, dependencies
├── mkdocs.yml                   # Documentation config
└── README.md                    # User-facing documentation
```

### Core Classes

#### 1. **Potential** (`core/potential.py`)

Represents a potential energy surface (PES) along configuration coordinate Q.

**Key Attributes**:
```python
Q: ndarray                 # Configuration coordinate grid (amu^0.5·Å)
E: ndarray                 # Potential energy on grid (eV)
Q0: float                  # Equilibrium position (amu^0.5·Å)
E0: float                  # Minimum energy (eV)
fit_func: Callable         # Fitted function E(Q)
fit_type: str              # "harmonic", "morse", "spline", "polynomial"
eigenvalues: ndarray       # Eigenvalues εₙ (eV)
eigenvectors: ndarray      # Eigenfunctions χₙ(Q) [shape: (nev, len(Q))]
```

**Key Methods**:
```python
# Creation
Potential.from_harmonic(hw, Q0, E0)      # Analytical harmonic potential
Potential.from_file(filename)            # Load from DFT data

# Fitting
potential.fit(fit_type="spline", order=4, smoothness=0.001)

# Solving
potential.solve(nev=60)                   # Solve Schrödinger equation

# I/O
potential.to_json(filename)
Potential.from_json(filename)
```

#### 2. **ConfigCoordinate** (`core/config_coord.py`)

Manages two-state capture calculation.

**Key Attributes**:
```python
pot_i: Potential           # Initial state (excited)
pot_f: Potential           # Final state (ground)
W: float                   # Electron-phonon coupling (eV)
g: int                     # Degeneracy factor
overlap_matrix: ndarray    # ⟨χᵢ|Q|χⱼ⟩
capture_coefficient: ndarray  # C(T) in cm³/s
```

**Key Methods**:
```python
cc = ConfigCoordinate(pot_i, pot_f, W=0.068)

# Calculate overlap matrix
cc.calculate_overlap(Q0=5.0, sigma=0.025, cutoff=0.25)

# Calculate capture coefficient
cc.calculate_capture_coefficient(
    volume=1e-21,
    temperature=np.linspace(100, 500, 50)
)

# Results
cc.capture_coefficient  # Array of C(T) values
```

#### 3. **ScanResult** (`analysis/parameter_scan.py`)

High-throughput parameter scanning.

**Key Functions**:
```python
from carriercapture.analysis.parameter_scan import scan_parameters

results = scan_parameters(
    dQ_range=(0, 25, 25),      # (min, max, n_points)
    dE_range=(0, 2.5, 10),
    hw_i=0.008,
    hw_f=0.008,
    W=0.068,
    volume=1e-21,
    temperature=300.0,
    n_jobs=-1                   # Parallel execution
)

# Access results
results.capture_coefficients  # 2D array [dQ, dE]
results.barrier_heights       # 2D array [dQ, dE]
```

---

## Common Workflows

### 1. Basic Harmonic Oscillator (Python API)

```python
import numpy as np
from carriercapture.core.potential import Potential
from carriercapture.core.config_coord import ConfigCoordinate

# Create potentials
pot_i = Potential.from_harmonic(hw=0.008, Q0=0.0, E0=0.5)
pot_f = Potential.from_harmonic(hw=0.008, Q0=10.5, E0=0.0)

# Solve Schrödinger equation
pot_i.solve(nev=180)
pot_f.solve(nev=60)

# Calculate capture coefficient
cc = ConfigCoordinate(pot_i=pot_i, pot_f=pot_f, W=0.068)
cc.calculate_overlap(Q0=5.0)
cc.calculate_capture_coefficient(
    volume=1e-21,
    temperature=np.linspace(100, 500, 50)
)

print(f"C(300K) = {cc.capture_coefficient[20]:.3e} cm³/s")
```

### 2. DFT Data → Capture Coefficient (CLI)

```bash
# Fit potential from DFT data
carriercapture fit excited_state.dat -f spline -o 4 -s 0.001 -O excited.json

# Solve Schrödinger equation
carriercapture solve excited.json -n 180 -O excited_solved.json
carriercapture solve ground.json -n 60 -O ground_solved.json

# Calculate capture coefficient
carriercapture capture config.yaml -V 1e-21 --temp-range 100 500 50
```

**config.yaml**:
```yaml
initial: excited_solved.json
final: ground_solved.json
W: 0.068
volume: 1.0e-21
temperature: [100, 500, 50]
Q0: 5.0
sigma: 0.025
cutoff: 0.25
```

### 3. Interactive Visualization Dashboard

```bash
carriercapture viz --port 8050
```

**Current Defaults** (as of 2026-01-18):
- `hw`: 0.02 eV (20 meV phonon energy)
- `nev`: 500 eigenvalues
- Display: Wavefunctions shown by default
- **Note**: "Show Eigenvalues" option removed to declutter interface

**Workflow**:
1. Click "Generate Harmonic" with default hw=0.02
2. Click "Solve Schrödinger" with default nev=500
3. Wavefunctions automatically displayed (ψ₀, ψ₅₀, ψ₁₀₀, ... ψ₄₅₀)
4. Adjust "Wavefunction Scaling" slider to tune amplitude

### 4. Integration with doped

```python
from carriercapture.io.doped_interface import (
    load_defect_entry,
    create_potential_from_doped
)

# Load defect entry from doped calculation
defect = load_defect_entry("path/to/defect.json.gz")

# Create potentials for charge state transition (0 → +1)
pot_i = create_potential_from_doped(defect, charge_state=0)
pot_f = create_potential_from_doped(defect, charge_state=+1)

# Continue with standard workflow
pot_i.solve(nev=180)
pot_f.solve(nev=60)
# ... calculate capture coefficient
```

---

## Development Practices

### Code Style

- **Formatter**: Black (line length: 100)
- **Linter**: Ruff (E, F, I, N, W checks)
- **Type hints**: Comprehensive (mypy strict mode)
- **Docstrings**: NumPy style

**Example**:
```python
def solve(
    self,
    nev: int | None = None,
    maxiter: int = 10000
) -> None:
    """
    Solve 1D Schrödinger equation for vibrational states.

    Parameters
    ----------
    nev : int, optional
        Number of eigenvalues to compute
    maxiter : int, default=10000
        Maximum ARPACK iterations

    Raises
    ------
    ValueError
        If potential function is not fitted

    Examples
    --------
    >>> pot.fit(fit_type="spline")
    >>> pot.solve(nev=60)
    >>> len(pot.eigenvalues)
    60
    """
```

### Testing

**Framework**: pytest with pytest-cov

**Structure**:
```
tests/
├── test_potential.py         # Potential class tests
├── test_schrodinger.py       # Solver validation (analytical solutions)
├── test_config_coord.py      # Capture coefficient workflows
├── test_parameter_scan.py    # High-throughput scanning
├── test_visualization.py     # Plotting functions
├── test_io.py                # File I/O
└── test_cli.py               # Command-line interface
```

**Run tests**:
```bash
pytest tests/ -v --cov                           # All tests with coverage
pytest tests/test_potential.py -v                # Specific module
pytest tests/ -k "harmonic" -v                   # Tests matching keyword
pytest tests/ -x --pdb                           # Stop on first failure, debug
```

**Current Status** (2026-01-18):
- 88 tests passing
- Core modules: >90% coverage
- Python 3.9-3.12 supported

### Git Workflow

**Commit Message Format**:
```
Brief description of change (50 chars max)

Longer explanation if needed, wrapped at 72 characters.
Can include bullet points:
- Change 1
- Change 2

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>
```

**Example**:
```bash
git add src/carriercapture/visualization/interactive.py
git commit -m "Update visualizer defaults for better harmonic well display

Changes:
- Set hw default to 0.02 eV (20 meV) for clearer eigenvalue spacing
- Set nev default to 500 (shows ψ₀, ψ₅₀, ψ₁₀₀, ... ψ₄₅₀)
- Remove 'Show Eigenvalues' option to declutter interface

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>"
```

### Documentation

**System**: MkDocs with Material theme

**Build**:
```bash
mkdocs serve              # Local preview at http://127.0.0.1:8000
mkdocs build              # Build static site
mkdocs gh-deploy          # Deploy to GitHub Pages
```

**Structure**:
- `docs/index.md`: Landing page
- `docs/getting-started/`: Installation, quick start, basic concepts
- `docs/user-guide/`: Detailed usage guides
- `docs/api/`: API reference (auto-generated from docstrings)
- `docs/examples/`: Example gallery

---

## Important Constants and Units

**File**: `src/carriercapture/_constants.py`

```python
# Physical constants
AMU = 931.4940954e6          # eV/c² - atomic mass unit
HBAR_C = 0.19732697e-6       # eV·m - ℏc
HBAR = 6.582119514e-16       # eV·s - ℏ
K_B = 8.6173303e-5           # eV/K - Boltzmann constant

# Convergence criteria
OCC_CUTOFF = 1e-5            # Max occupation for partition function convergence
```

**Unit Conventions**:
- **Energy**: eV (electronvolts)
- **Length**: Å (angstroms) or cm
- **Mass**: amu (atomic mass units)
- **Temperature**: K (kelvin)
- **Configuration coordinate Q**: amu^0.5·Å (mass-weighted displacement)
- **Capture coefficient C**: cm³/s

---

## Typical Parameter Values

### Physical Systems

**Sn_Zn in ZnO** (standard example):
```python
hw = 0.008 eV           # 8 meV phonon
dQ = 10.5               # amu^0.5·Å shift
dE = 0.5 eV             # Energy difference
W = 0.068 eV            # Electron-phonon coupling
volume = 1e-21 cm³      # Supercell volume
temperature = 300 K     # Room temperature
nev_initial = 180       # Initial state eigenvalues
nev_final = 60          # Final state eigenvalues
```

**DX Centers in GaAs**:
```python
hw = 0.03 eV            # 30 meV (stronger mode)
dQ = 10.0
dE = 1.0 eV
W = 0.2 eV
```

### Visualization Defaults

**Current settings** (optimized for demo):
```python
hw = 0.02 eV            # 20 meV - good visual spacing
nev = 500               # Shows ~10 distributed wavefunctions
Q_range = (-25, 25)     # Wide enough for wavefunction decay
npoints = 5000          # Smooth potential curves
```

**Harmonic oscillator eigenvalues**:
```
E_n = ℏω * (n + 1/2)
```
- For hw=0.02 eV: E₀=0.01, E₁=0.03, E₂=0.05, ..., E₄₉₉=10.0 eV
- With nev=500, shows: ψ₀, ψ₅₀, ψ₁₀₀, ψ₁₅₀, ... ψ₄₅₀ (every 50th state)

---

## Recent Changes & Evolution

### 2026-01-18: Visualizer Improvements

**Commits**: `3c2641a`, `4981a77`, `3d81883`

**Changes**:
1. **Fix: Harmonic potential serialization**
   - Added `fit_params` to serialization
   - Recreate `fit_func` on deserialization (harmonic: analytical, others: interpolation)
   - **Impact**: "Generate Harmonic" → "Solve" now works in dashboard

2. **Fix: Wavefunction distribution**
   - Changed from showing first 5 states (clustered at bottom)
   - Now shows every Nth state (up to 10 total, distributed across energy range)
   - **Impact**: Wavefunctions fill the harmonic well properly

3. **Update: Dashboard defaults**
   - hw: 0.008 → 0.02 eV (clearer spacing)
   - nev: 60 → 500 (wider energy coverage)
   - nev max: 500 → 1000 (allow higher precision)
   - Removed "Show Eigenvalues" checkbox (declutter)
   - Default display: Show wavefunctions (not eigenvalues)

### 2026-01-18: Documentation & Dash API

**Commits**: `dd994ba`, `33f6880`

**Changes**:
- Removed empty Theory/Tutorials/Development sections from docs
- Added Huang & Rhys (1950) paper URL as clickable link
- Fixed deprecated Dash API: `app.run_server()` → `app.run()`
- Deployed documentation to GitHub Pages

### Benchmark Validation

**File**: `benchmarks/benchmark_sn_zn.py`

**Results** (vs CarrierCapture.jl):
- Initial eigenvalues: 0.005% difference ✓
- Final eigenvalues: 0.02% difference ✓
- Capture coefficient (300K): 1.5% difference ✓

**Conclusion**: Python matches Julia within ~1-2% (floating-point precision)

---

## Common Gotchas & Tips

### 1. Potential Must Be Fitted Before Solving

**Error**: `ValueError: Must fit potential before solving`

**Cause**: Calling `potential.solve()` when `fit_func` is None

**Solutions**:
```python
# Option A: Use from_harmonic (auto-fitted)
pot = Potential.from_harmonic(hw=0.008, Q0=0.0, E0=0.5)
pot.solve(nev=60)  # ✓ Works

# Option B: Manually fit data
pot = Potential()
pot.Q_data = data[:, 0]
pot.E_data = data[:, 1]
pot.fit(fit_type="spline", order=4, smoothness=0.001)  # Required!
pot.solve(nev=60)  # ✓ Works
```

### 2. Eigenvalue Convergence

**Rule of thumb**: `nev > 4 * kB * T / hw`

At T=300K (kB*T ≈ 0.026 eV):
- hw=0.008 eV → need nev > 13 (use 60+ for safety)
- hw=0.020 eV → need nev > 5 (use 30+ for safety)

**Check convergence**:
```python
from carriercapture._constants import OCC_CUTOFF

# Calculate occupation of highest state
p_max = cc.pot_i.boltzmann_occupation(temperature)[-1]
if p_max > OCC_CUTOFF:
    print(f"Warning: Need more eigenvalues! p_max = {p_max:.2e}")
```

### 3. Units in Configuration Coordinate

Q is **mass-weighted displacement**: `Q = √m * ΔR`

**Example**: Moving an atom with mass 50 amu by 0.5 Å:
```python
m = 50  # amu
delta_R = 0.5  # Å
Q = np.sqrt(m) * delta_R  # ≈ 3.54 amu^0.5·Å
```

### 4. Serialization Caveat

**Problem**: `fit_func` (Python function) cannot be JSON-serialized

**Solution**: We serialize `fit_params` and `E` values, then recreate:
- Harmonic: Reconstruct analytical function from parameters
- Others: Use cubic interpolation from stored `E` values

**Impact**: Slight loss of precision for non-harmonic fits after deserialization

### 5. Dashboard State Management

The Dash dashboard uses `dcc.Store` to preserve state between callbacks:
```python
# Potential data stored as JSON-serializable dict
"potential-store-initial": serialize_potential(pot)

# Scan results stored separately
"scan-results-store": serialize_scan_results(results)
```

**Tip**: When debugging callbacks, check `ctx.triggered_id` to see which button was clicked.

### 6. Parallel Execution

Parameter scanning uses `joblib` for parallelization:
```python
from carriercapture.analysis.parameter_scan import scan_parameters

results = scan_parameters(
    ...,
    n_jobs=-1   # Use all CPU cores
)
```

**Warning**: Progress bars may behave oddly with `n_jobs > 1` (joblib limitation).

---

## CLI Command Reference

### Main Commands

```bash
carriercapture fit <input> [OPTIONS]          # Fit potential energy surface
carriercapture solve <input> [OPTIONS]        # Solve Schrödinger equation
carriercapture capture <config> [OPTIONS]     # Calculate capture coefficient
carriercapture scan [OPTIONS]                 # Parameter scanning
carriercapture scan-plot <input> [OPTIONS]    # Visualize scan results
carriercapture viz [OPTIONS]                  # Launch Dash dashboard
carriercapture plot <input> [OPTIONS]         # Generate static plots
```

### Key Options

**fit**:
```bash
-f, --fit-type [spline|harmonic|morse|polynomial]
-o, --order INTEGER              # Spline/polynomial order
-s, --smoothness FLOAT           # Spline smoothing parameter
-O, --output PATH                # Save fitted potential
```

**solve**:
```bash
-n, --nev INTEGER                # Number of eigenvalues
-O, --output PATH                # Save solved potential
```

**capture**:
```bash
-V, --volume FLOAT               # Supercell volume (cm³)
--temp-range MIN MAX NPOINTS     # Temperature range (K)
-o, --output PATH                # Save results (NPZ or HDF5)
```

**viz**:
```bash
--port INTEGER                   # Server port (default: 8050)
--host TEXT                      # Host address (default: 127.0.0.1)
--debug                          # Enable debug mode
```

---

## API Quick Reference

### Core Imports

```python
from carriercapture import Potential, ConfigCoordinate
from carriercapture.core.schrodinger import solve_schrodinger_1d
from carriercapture.analysis.parameter_scan import scan_parameters
from carriercapture.visualization.static import (
    plot_potential,
    plot_config_coordinate,
    plot_scan_heatmap
)
```

### Potential Creation

```python
# Harmonic (analytical)
pot = Potential.from_harmonic(hw=0.008, Q0=0.0, E0=0.5)

# From file
pot = Potential.from_file("data.dat")
pot.fit(fit_type="spline", order=4, smoothness=0.001)

# From doped
from carriercapture.io.doped_interface import create_potential_from_doped
pot = create_potential_from_doped(defect_entry, charge_state=0)
```

### Solving

```python
pot.solve(nev=180)

# Access results
pot.eigenvalues        # Array of εₙ (eV)
pot.eigenvectors       # Array of χₙ(Q), shape: (nev, len(Q))
```

### Capture Calculation

```python
cc = ConfigCoordinate(pot_i, pot_f, W=0.068, g=1)
cc.calculate_overlap(Q0=5.0, sigma=0.025, cutoff=0.25)
cc.calculate_capture_coefficient(
    volume=1e-21,
    temperature=np.linspace(100, 500, 50)
)

# Results
cc.capture_coefficient  # Array of C(T) in cm³/s
cc.overlap_matrix       # Matrix of ⟨χᵢ|Q|χⱼ⟩
```

### Visualization

```python
from carriercapture.visualization.static import plot_potential

fig = plot_potential(pot, show_eigenvalues=True, show_wavefunctions=True)
fig.show()
```

---

## Example Notebooks

Located in `examples/notebooks/`:

1. **`01_harmonic_sn_zn.ipynb`** - Complete workflow with harmonic potentials (Sn_Zn in ZnO)
2. **`03_parameter_scan.ipynb`** - High-throughput screening example

**To run**:
```bash
pip install carriercapture[notebook]
cd examples/notebooks
jupyter notebook
```

---

## Troubleshooting

### ARPACK Convergence Issues

**Error**: `RuntimeError: ARPACK failed to converge`

**Solutions**:
1. Increase `maxiter`: `pot.solve(nev=60, maxiter=20000)`
2. Reduce `nev`: Try fewer eigenvalues
3. Check potential quality: Ensure smooth, well-behaved fit
4. Expand Q_range: Make sure potential is well-defined

### Import Errors

**Error**: `ModuleNotFoundError: No module named 'carriercapture'`

**Solution**: Install in editable mode:
```bash
pip install -e ".[dev]"
```

### Dash Dashboard Won't Start

**Error**: `app.run_server has been replaced by app.run`

**Solution**: Update to latest version (fixed as of commit `33f6880`)

### Slow Performance

**Tips**:
1. Reduce `npoints` in potential (5000 → 3000)
2. Use parallel execution: `n_jobs=-1`
3. Use NumPy/SciPy vectorized operations
4. Profile with `cProfile` to identify bottlenecks

---

## Dependencies

**Core** (always installed):
```
numpy>=1.21
scipy>=1.7
pandas>=1.3
pyyaml>=5.4
click>=8.0
rich>=10.0
plotly>=5.0
dash>=2.0
pymatgen>=2022.0
joblib>=1.0
```

**Development** (`pip install -e ".[dev]"`):
```
pytest, pytest-cov, pytest-xdist
mypy
black, ruff
pre-commit
```

**Documentation** (`pip install -e ".[docs]"`):
```
mkdocs, mkdocs-material
mkdocstrings[python]
mkdocs-git-revision-date-localized-plugin
```

**Notebooks** (`pip install -e ".[notebook]"`):
```
jupyter, ipywidgets
```

**doped Integration** (`pip install -e ".[doped]"`):
```
doped>=2.0, monty>=2024.0
```

---

## Resources

### Documentation
- **GitHub**: https://github.com/WMD-group/CarrierCapture.py
- **Docs**: https://WMD-group.github.io/CarrierCapture.py/
- **Original Julia**: https://github.com/WMD-group/CarrierCapture.jl

### Key Papers
- [Alkauskas et al. (2014) - Phys. Rev. B 90, 075202](https://doi.org/10.1103/PhysRevB.90.075202)
- [Huang & Rhys (1950) - Proc. R. Soc. Lond. A 204, 406](https://royalsocietypublishing.org/rspa/article/204/1078/406/8369/Theory-of-light-absorption-and-non-radiative)

### Related Tools
- **doped**: Defect calculations with pymatgen/VASP
- **pymatgen**: Materials analysis framework
- **phonopy**: Phonon calculations

---

## Contributing

### Before Submitting PR

1. **Run tests**: `pytest tests/ -v --cov`
2. **Format code**: `black src/ tests/`
3. **Lint**: `ruff check src/ tests/`
4. **Type check**: `mypy src/`
5. **Update docs**: If adding features, update relevant .md files
6. **Add tests**: New features need test coverage

### PR Template

```markdown
## Description
Brief description of changes

## Type of Change
- [ ] Bug fix
- [ ] New feature
- [ ] Breaking change
- [ ] Documentation update

## Testing
- [ ] Tests added/updated
- [ ] All tests passing
- [ ] Code formatted (black)
- [ ] Linted (ruff)

## Checklist
- [ ] Updated CHANGELOG.md
- [ ] Updated documentation
- [ ] Added example if appropriate
```

---

## Version History

- **1.0.0** (2026-01): Complete rewrite from Julia, feature parity achieved
- Benchmark validation: Within 1-2% of CarrierCapture.jl
- Dashboard improvements: Fixed serialization, wavefunction distribution, defaults
- Documentation deployed to GitHub Pages

---

## Notes for AI Agents

### When Modifying Code

1. **Always read files before editing**: Never propose changes without reading the file
2. **Respect existing patterns**: Follow NumPy docstring style, type hints, Black formatting
3. **Run tests after changes**: Ensure nothing breaks
4. **Update comments**: Keep code and comments in sync
5. **Use specialized tools**: Use Edit/Read/Write tools, not bash cat/sed/awk

### When Adding Features

1. **Check existing examples**: Look at notebooks, tests for patterns
2. **Consider physical validity**: Units, typical values, energy conservation
3. **Add tests**: New features must have test coverage
4. **Update documentation**: Add to relevant user-guide or API docs
5. **Benchmark if relevant**: Compare against Julia or analytical solutions

### When Debugging

1. **Read error messages carefully**: ARPACK failures, convergence issues
2. **Check parameter ranges**: Are eigenvalues reasonable? Is nev sufficient?
3. **Verify units**: Configuration coordinate is mass-weighted (amu^0.5·Å)
4. **Use pytest -x --pdb**: Stop on first failure and debug interactively
5. **Profile performance**: Use cProfile or line_profiler for bottlenecks

### Important Files to Understand

- `core/potential.py`: Central class, ~1000 lines, handles PES and solving
- `core/config_coord.py`: Capture coefficient calculation, ~600 lines
- `visualization/interactive.py`: Dash dashboard, ~1200 lines
- `tests/test_config_coord.py`: Complete workflow tests, line 213-250 key example

---

**Last Updated**: 2026-01-18
**Maintained By**: CarrierCapture Team + Claude (Anthropic)
**Questions?**: See GitHub Issues or Discussions
