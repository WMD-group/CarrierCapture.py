# Jupyter Notebook Examples

Interactive tutorials and examples in Jupyter notebook format.

## Overview

CarrierCapture includes Jupyter notebooks that demonstrate complete workflows from data loading to results visualization. These notebooks are:

- **Interactive** - Run code cells, modify parameters, see results immediately
- **Educational** - Detailed explanations with theory and practice
- **Reproducible** - All data included, ready to run
- **Extensible** - Easy to adapt for your own calculations

---

## Available Notebooks

### 📓 01: Harmonic Oscillator (Sn in ZnO)

**File**: `examples/notebooks/01_harmonic_sn_zn.ipynb`

**Status**: ✅ Available

**Topics Covered**:
- Creating harmonic potentials
- Solving Schrödinger equation
- Calculating Franck-Condon overlaps
- Computing capture coefficient
- Temperature-dependent analysis
- Arrhenius plots

**Learning Objectives**:
- Understand basic workflow
- Master Python API
- Interpret results physically
- Create publication plots

**Data**: Harmonic approximation for Sn substituting Zn in ZnO

**Estimated Time**: 30 minutes

#### Quick Start

```bash
# Navigate to examples
cd examples/notebooks

# Launch Jupyter
jupyter notebook 01_harmonic_sn_zn.ipynb
```

Or run cells online:
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/WMD-group/CarrierCapture.py/blob/main/examples/notebooks/01_harmonic_sn_zn.ipynb)

#### Key Code Snippets

**Creating potentials:**
```python
from carriercapture.core import Potential

pot_i = Potential.from_harmonic(hw=0.008, Q0=0.0, E0=0.5)
pot_f = Potential.from_harmonic(hw=0.008, Q0=10.5, E0=0.0)
```

**Solving:**
```python
pot_i.solve(nev=180)
pot_f.solve(nev=60)
```

**Capture calculation:**
```python
from carriercapture.core import ConfigCoordinate

cc = ConfigCoordinate(pot_i, pot_f, W=0.068, degeneracy=1)
cc.calculate_overlap(Q0=5.0, cutoff=0.25, sigma=0.025)
cc.calculate_capture_coefficient(volume=1e-21, temperature=temperature)
```

---

### 📓 02: Anharmonic Potential (DX Center)

**File**: `examples/notebooks/02_anharmonic_dx_center.ipynb`

**Status**: ✅ Available

**Topics Covered**:
- Generating/loading Q-E data for anharmonic potentials
- Spline and Morse fitting for anharmonic potentials
- Quality assessment of fits (R², RMSE)
- Anharmonic effects on energy level spacing
- Comparison with harmonic approximation
- Impact on carrier capture coefficients

**Learning Objectives**:
- Work with anharmonic potential data
- Master fitting techniques (spline, Morse, harmonic)
- Understand how anharmonicity affects capture rates
- Validate fit quality quantitatively

**Data**: Synthetic DX-center-like Morse potential data

**Estimated Time**: 30 minutes

#### Quick Start

```bash
jupyter notebook 02_anharmonic_dx_center.ipynb
```

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/WMD-group/CarrierCapture.py/blob/main/examples/notebooks/02_anharmonic_dx_center.ipynb)

#### Key Code Snippets

```python
# Generate synthetic anharmonic data
Q_data, E_data = read_csv_data('dx_center_data.csv')

# Create and fit potential
pot = Potential(Q_data=Q_data, E_data=E_data)
pot.fit(fit_type='spline', order=4, smoothness=0.001)

# Assess fit quality
Q_fit = np.linspace(Q_data.min(), Q_data.max(), 500)
E_fit = pot(Q_fit)

# Compare with harmonic
pot_harmonic = Potential(Q_data=Q_data, E_data=E_data)
pot_harmonic.fit(fit_type='harmonic')

# Analyze anharmonic effects
```

---

### 📓 03: Parameter Scanning

**File**: `examples/notebooks/03_parameter_scan.ipynb`

**Status**: ✅ Available

**Topics Covered**:
- Setting up parameter scans
- Defining (ΔQ, ΔE) grids
- Parallel execution
- Result analysis and visualization
- Identifying optimal parameters
- Materials screening workflow

**Learning Objectives**:
- Perform high-throughput screening
- Optimize computational efficiency
- Analyze 2D parameter space
- Extract design principles

**Data**: Systematic scan over defect parameter space

**Estimated Time**: 30 minutes

#### Quick Start

```bash
jupyter notebook 03_parameter_scan.ipynb
```

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/WMD-group/CarrierCapture.py/blob/main/examples/notebooks/03_parameter_scan.ipynb)

#### Key Code Snippets

**Defining scan:**
```python
from carriercapture.analysis import ParameterScanner, ScanParameters

params = ScanParameters(
    dQ_range=(0, 25, 25),
    dE_range=(0, 2.5, 10),
    hbar_omega_i=0.008,
    temperature=300.0,
    volume=1e-21
)
```

**Running scan:**
```python
scanner = ParameterScanner(params, verbose=True)
results = scanner.run_harmonic_scan(n_jobs=-1, show_progress=True)
```

**Visualizing:**
```python
from carriercapture.visualization import plot_scan_heatmap

fig = plot_scan_heatmap(results, log_scale=True)
fig.show()
```

---

### 📓 04: Interactive Visualization

**File**: `examples/notebooks/04_interactive_viz.ipynb`

**Status**: ✅ Available

**Topics Covered**:
- Launching Dash dashboard from notebook and CLI
- Dashboard features tour (fitting, capture, scanning, comparison)
- Programmatic control of visualizations
- Exporting figures (HTML, PNG, PDF)
- Exporting data (JSON, NPZ)
- Customizing themes and styles

**Learning Objectives**:
- Launch and navigate the interactive dashboard
- Use static plotting functions programmatically
- Export publication-quality figures
- Customize visualization themes

**Data**: Example potentials created in notebook

**Estimated Time**: 20 minutes

#### Quick Start

```bash
jupyter notebook 04_interactive_viz.ipynb
```

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/WMD-group/CarrierCapture.py/blob/main/examples/notebooks/04_interactive_viz.ipynb)

#### Key Code Snippets

```python
# Launch dashboard
from carriercapture.visualization.interactive import create_app, run_server

app = create_app()
app.run(port=8050, jupyter_mode="tab")

# Or from CLI: carriercapture viz --port 8050
```

---

## Running Notebooks

### Local Installation

```bash
# Install CarrierCapture with notebook dependencies
pip install carriercapture[dev]

# Navigate to examples
cd examples/notebooks

# Launch Jupyter
jupyter notebook

# Or use JupyterLab
jupyter lab
```

### Online (Google Colab)

1. Click the "Open in Colab" badge on any notebook
2. Notebooks will open in Google Colab (free, no installation)
3. First cell installs CarrierCapture automatically
4. Run all cells sequentially

**Note**: Colab uses Python 3.10+ and all dependencies are pre-installed.

### VS Code

```bash
# Install VS Code Python extension
# Open notebook file (*.ipynb)
# Select Python kernel
# Run cells with Shift+Enter
```

---

## Notebook Structure

All notebooks follow a consistent structure:

### 1. Setup
```python
# Imports
import numpy as np
from carriercapture.core import Potential, ConfigCoordinate
from carriercapture.visualization import plot_potential, plot_capture_coefficient

# Configure matplotlib (if used)
%matplotlib inline
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = (10, 6)
```

### 2. Data Loading
```python
# Load or create data
# Example: from file
Q_data, E_data = read_csv_data('potential_data.csv')

# Example: analytical
pot = Potential.from_harmonic(hw=0.008, Q0=0.0, E0=0.0)
```

### 3. Calculation
```python
# Fit potential
pot.fit(fit_type='spline', order=4, smoothness=0.001)

# Solve Schrödinger equation
pot.solve(nev=180)

# Calculate capture
cc = ConfigCoordinate(pot_i, pot_f, W=0.205)
cc.calculate_overlap(Q0=10.0)
cc.calculate_capture_coefficient(volume=1e-21, temperature=temperature)
```

### 4. Analysis
```python
# Extract results
C_300K = cc.capture_coefficient[temperature == 300][0]
print(f"C(300K) = {C_300K:.3e} cm³/s")

# Statistical analysis
print(f"Max C: {cc.capture_coefficient.max():.3e} cm³/s")
print(f"Min C: {cc.capture_coefficient.min():.3e} cm³/s")
```

### 5. Visualization
```python
# Create plots
fig1 = plot_potential(pot, show_wavefunctions=True)
fig2 = plot_capture_coefficient(cc, arrhenius=True)

# Display
fig1.show()
fig2.show()

# Save
fig1.write_html('potential.html')
fig2.write_image('capture.png')
```

### 6. Discussion
- Interpretation of results
- Physical insights
- Comparison with literature
- Limitations and caveats

---

## Tips for Using Notebooks

### 1. Run Cells Sequentially

Notebooks are designed to run top-to-bottom:
```python
# Cell 1: Imports
# Cell 2: Data loading
# Cell 3: Calculation
# Cell 4: Results
```

**Don't skip cells** - later cells depend on earlier ones.

### 2. Modify Parameters

Encouraged! Change parameters to see effects:

```python
# Original
pot = Potential.from_harmonic(hw=0.008, Q0=0.0, E0=0.0)

# Try different phonon energy
pot = Potential.from_harmonic(hw=0.010, Q0=0.0, E0=0.0)  # 10 meV instead of 8 meV

# See how results change
```

### 3. Save Your Work

```python
# Save modified notebook: File → Save

# Export results
from carriercapture.io import save_potential
save_potential(pot, 'my_potential.json')

# Export figures
fig.write_html('my_figure.html')
```

### 4. Add Your Own Cells

Insert cells to:
- Explore data
- Try different parameters
- Create additional plots
- Add notes

```python
# Add markdown cell for notes
# Add code cell for experiments
```

### 5. Use Keyboard Shortcuts

| Action | Shortcut |
|--------|----------|
| Run cell | Shift+Enter |
| Run cell (stay) | Ctrl+Enter |
| Insert cell above | A |
| Insert cell below | B |
| Delete cell | D, D |
| Change to markdown | M |
| Change to code | Y |

---

## Troubleshooting

### Import Errors

```python
# If you get: ModuleNotFoundError: No module named 'carriercapture'

# Solution 1: Install package
!pip install carriercapture

# Solution 2: Install from source
!pip install git+https://github.com/WMD-group/CarrierCapture.py.git

# Solution 3: Local installation
# cd to repo root, then:
!pip install -e .
```

### Plotting Issues

```python
# If plots don't show:

# For Jupyter Notebook
%matplotlib inline

# For Plotly (interactive)
import plotly.io as pio
pio.renderers.default = 'notebook'  # or 'colab' for Google Colab

# For JupyterLab
# Install extension:
# jupyter labextension install jupyterlab-plotly
```

### Memory Issues (Large Scans)

```python
# If notebook crashes during parameter scan:

# Solution 1: Reduce grid size
params = ScanParameters(
    dQ_range=(0, 25, 15),  # Reduce from 25 to 15 points
    dE_range=(0, 2.5, 8),  # Reduce from 10 to 8 points
    # ...
)

# Solution 2: Use fewer parallel jobs
results = scanner.run_harmonic_scan(n_jobs=4)  # Instead of -1

# Solution 3: Save intermediate results
results = scanner.run_harmonic_scan(n_jobs=-1)
results.save('intermediate_results.npz')
```

---

## Contributing Notebooks

Have a useful example? Contribute it!

### Notebook Guidelines

1. **Self-contained** - Include all necessary imports and data
2. **Well-documented** - Markdown cells explaining each step
3. **Tested** - Run all cells before submitting
4. **Clean output** - Clear output before committing (optional)
5. **Licensed** - MIT license (same as package)

### Submission Process

```bash
# Fork repository
# Create branch
git checkout -b add-notebook-name

# Add notebook
cp your_notebook.ipynb examples/notebooks/05_your_topic.ipynb

# Add data if needed
cp your_data.csv examples/data/

# Test notebook
jupyter notebook examples/notebooks/05_your_topic.ipynb
# Run all cells

# Commit and push
git add examples/notebooks/05_your_topic.ipynb
git commit -m "Add notebook: Your Topic"
git push origin add-notebook-name

# Open pull request on GitHub
```

### Notebook Template

Use this template for new notebooks:

```markdown
# Notebook Title

**Author**: Your Name
**Date**: YYYY-MM-DD
**Topics**: Topic 1, Topic 2, Topic 3

## Overview

Brief description of what this notebook demonstrates.

## Learning Objectives

- Objective 1
- Objective 2
- Objective 3

## Setup

\```python
# Imports
import numpy as np
from carriercapture.core import Potential, ConfigCoordinate
\```

## Data

Description of data used.

## Calculation

Step-by-step calculation with explanations.

## Results

Analysis and interpretation.

## Conclusion

Summary and key takeaways.

## References

- Reference 1
- Reference 2
```

---

## See Also

- **[Getting Started: Quick Start](../getting-started/quick-start.md)** - 5-minute introduction
- **[User Guide](../user-guide/potentials.md)** - Detailed usage guides
- **[Examples: Gallery](gallery.md)** - Example visualizations
- **[API Reference](../api/core.md)** - Complete API documentation
