# Example Gallery

Showcase of visualizations and calculations you can create with CarrierCapture.

## Overview

This gallery demonstrates the types of plots, diagrams, and analyses possible with CarrierCapture.py. All examples use real data and are reproducible with the code snippets provided.

---

## Potential Energy Surfaces

### Harmonic Oscillator

Simple harmonic potential with wavefunctions:

```python
from carriercapture.core import Potential
from carriercapture.visualization import plot_potential

# Create harmonic potential
pot = Potential.from_harmonic(hw=0.008, Q0=0.0, E0=0.0, Q_range=(-20, 20))
pot.solve(nev=60)

# Plot with wavefunctions
fig = plot_potential(
    pot,
    title="Harmonic Oscillator (ℏω = 8 meV)",
    show_wavefunctions=True,
    n_states=10,
    show_eigenvalues=True
)
fig.show()
```

**Features**:
- Equally spaced energy levels (quantum harmonic oscillator)
- Gaussian-like wavefunctions
- Ground state localized at potential minimum

---

### Anharmonic Potential (Spline Fit)

Fitted potential from DFT data:

```python
import numpy as np
from carriercapture.core import Potential

# Load Q-E data from DFT
Q_data = np.array([0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20])
E_data = np.array([0.5, 0.4, 0.25, 0.12, 0.05, 0.01, 0.0, 0.02, 0.05, 0.09, 0.14])

# Create and fit potential
pot = Potential(Q_data=Q_data, E_data=E_data)
pot.fit(fit_type='spline', order=4, smoothness=0.001)
pot.solve(nev=60)

# Plot
from carriercapture.visualization import plot_potential
fig = plot_potential(
    pot,
    title="Anharmonic Potential (DFT + Spline Fit)",
    show_wavefunctions=True,
    n_states=12
)
fig.show()
```

**Features**:
- Smooth spline interpolation through DFT points
- Unequally spaced energy levels (anharmonicity)
- Wavefunctions adapt to potential shape

---

## Configuration Coordinate Diagrams

### Sn in ZnO (q=0 → q=+1)

Classic CC diagram for carrier capture:

```python
from carriercapture.io import load_potential_from_file
from carriercapture.core import Potential
from carriercapture.visualization import plot_configuration_coordinate

# Load solved potentials
pot_i = Potential.from_dict(load_potential_from_file('Sn_Zn_excited.json'))
pot_f = Potential.from_dict(load_potential_from_file('Sn_Zn_ground.json'))

# Create CC diagram
fig = plot_configuration_coordinate(
    pot_i,
    pot_f,
    Q0=10.0,
    show_crossing=True,
    show_wavefunctions=True,
    n_states=8,
    title="Sn_Zn in ZnO: Carrier Capture"
)
fig.show()
```

**Features**:
- Two potential energy surfaces (initial and final states)
- Crossing point between potentials
- Franck-Condon transitions (vertical arrows)
- Configuration coordinate offset (ΔQ)

---

### DX Center in AlGaAs

Deep defect with large structural relaxation:

```python
# DX center: large ΔQ ~ 20 amu^0.5·Å
fig = plot_configuration_coordinate(
    pot_dx_initial,
    pot_dx_final,
    Q0=20.0,
    show_crossing=True,
    show_wavefunctions=True,
    title="DX Center: Large Lattice Relaxation"
)
fig.show()
```

**Features**:
- Large configuration coordinate offset
- Deep crossing point (high barrier)
- Wide wavefunction spread

---

## Capture Coefficient Analysis

### Arrhenius Plot

Temperature-dependent capture coefficient:

```python
from carriercapture.core import ConfigCoordinate
from carriercapture.visualization import plot_capture_coefficient
import numpy as np

# Calculate capture coefficient
cc = ConfigCoordinate(pot_i, pot_f, W=0.205, degeneracy=1)
cc.calculate_overlap(Q0=10.0, cutoff=0.25, sigma=0.025)

temperature = np.linspace(100, 500, 50)
cc.calculate_capture_coefficient(volume=1e-21, temperature=temperature)

# Arrhenius plot
fig = plot_capture_coefficient(
    cc,
    arrhenius=True,
    show_activation_energy=True,
    title="Capture Coefficient: Sn_Zn in ZnO"
)
fig.show()
```

**Features**:
- Log scale (log₁₀ C vs 1000/T)
- Linear region → activation energy
- Plateau at low temperature (tunneling)
- Temperature labels on top axis

---

### Multiple Transitions Comparison

Compare capture for different charge state transitions:

```python
import plotly.graph_objects as go

# Calculate for multiple transitions
transitions = {
    'q=0→+1': cc_01,
    'q=+1→+2': cc_12,
    'q=-1→0': cc_m10,
}

fig = go.Figure()

for label, cc in transitions.items():
    fig.add_trace(go.Scatter(
        x=1000/cc.temperature,
        y=np.log10(cc.capture_coefficient),
        mode='lines+markers',
        name=label,
        line=dict(width=2)
    ))

fig.update_layout(
    title="Capture Coefficient: Multi-Transition Comparison",
    xaxis_title="1000/T (K⁻¹)",
    yaxis_title="log₁₀(C) [cm³/s]",
    template='plotly_white'
)
fig.show()
```

**Features**:
- Multiple curves on same plot
- Direct comparison of capture rates
- Identify fastest/slowest transitions

---

## Parameter Scanning

### 2D Heatmap (ΔQ vs ΔE)

High-throughput screening:

```python
from carriercapture.analysis import ParameterScanner, ScanParameters
from carriercapture.visualization import plot_scan_heatmap

# Run parameter scan
params = ScanParameters(
    dQ_range=(0, 25, 25),
    dE_range=(0, 2.5, 10),
    hbar_omega_i=0.008,
    hbar_omega_f=0.008,
    temperature=300.0,
    volume=1e-21
)

scanner = ParameterScanner(params)
results = scanner.run_harmonic_scan(n_jobs=-1)

# Plot heatmap
fig = plot_scan_heatmap(
    results,
    log_scale=True,
    title="Materials Screening: C(300K) Landscape",
    colorscale='Viridis'
)
fig.show()
```

**Features**:
- 2D grid: ΔQ (structural) vs ΔE (electronic)
- Log scale color for capture coefficient
- Identify "sweet spots" for fast capture
- Guide defect engineering

---

### Contour Plot

Iso-capture contours:

```python
fig = plot_scan_heatmap(
    results,
    plot_type='contour',
    log_scale=True,
    title="Iso-Capture Contours"
)
fig.show()
```

**Features**:
- Contour lines for constant capture
- Easy to identify optimal parameter regions
- Interpolated smooth transitions

---

## Overlap Matrix Analysis

### Franck-Condon Factors

Visualize wavefunction overlaps:

```python
from carriercapture.visualization import plot_overlap_matrix

# Calculate overlaps
cc.calculate_overlap(Q0=10.0, cutoff=0.25, sigma=0.025)

# Plot overlap matrix
fig = plot_overlap_matrix(
    cc,
    log_scale=True,
    title="Franck-Condon Overlap Matrix",
    colorscale='RdBu'
)
fig.show()
```

**Features**:
- Heatmap: initial states (y-axis) vs final states (x-axis)
- Log scale shows wide range of overlaps
- Diagonal pattern for small ΔQ
- Off-diagonal for large ΔQ (multiphonon)

---

## Eigenvalue Spectra

### Energy Level Diagrams

Compare phonon level spacing:

```python
from carriercapture.visualization import plot_eigenvalue_spectrum

fig = plot_eigenvalue_spectrum(
    pot,
    max_states=20,
    show_spacing=True,
    title="Vibrational Energy Levels"
)
fig.show()
```

**Features**:
- Horizontal lines = energy levels
- Spacing annotations
- Check harmonic vs anharmonic character

---

## Multi-Panel Figures

### Complete Analysis Figure

Publication-ready composite:

```python
from plotly.subplots import make_subplots
import plotly.graph_objects as go

# Create 2×2 subplot
fig = make_subplots(
    rows=2, cols=2,
    subplot_titles=(
        'Configuration Coordinate Diagram',
        'Capture Coefficient',
        'Overlap Matrix',
        'Parameter Scan'
    ),
    specs=[
        [{'type': 'xy'}, {'type': 'xy'}],
        [{'type': 'heatmap'}, {'type': 'heatmap'}]
    ]
)

# Panel 1: CC diagram
fig.add_trace(
    go.Scatter(x=pot_i.Q, y=pot_i.E, mode='lines', name='Initial'),
    row=1, col=1
)
fig.add_trace(
    go.Scatter(x=pot_f.Q, y=pot_f.E, mode='lines', name='Final'),
    row=1, col=1
)

# Panel 2: Capture coefficient
fig.add_trace(
    go.Scatter(
        x=1000/cc.temperature,
        y=np.log10(cc.capture_coefficient),
        mode='lines+markers',
        name='C(T)'
    ),
    row=1, col=2
)

# Panel 3: Overlap matrix
fig.add_trace(
    go.Heatmap(
        z=np.log10(np.abs(cc.overlap_matrix) + 1e-30),
        colorscale='Viridis',
        showscale=False
    ),
    row=2, col=1
)

# Panel 4: Parameter scan
fig.add_trace(
    go.Heatmap(
        x=scan_results.dE_grid,
        y=scan_results.dQ_grid,
        z=np.log10(scan_results.capture_coefficients),
        colorscale='Viridis'
    ),
    row=2, col=2
)

# Update layout
fig.update_layout(
    height=900,
    width=1200,
    title_text="Complete Carrier Capture Analysis",
    showlegend=False
)

fig.update_xaxes(title_text="Q (amu^0.5·Å)", row=1, col=1)
fig.update_xaxes(title_text="1000/T (K⁻¹)", row=1, col=2)
fig.update_xaxes(title_text="Final state", row=2, col=1)
fig.update_xaxes(title_text="ΔE (eV)", row=2, col=2)

fig.update_yaxes(title_text="E (eV)", row=1, col=1)
fig.update_yaxes(title_text="log₁₀(C)", row=1, col=2)
fig.update_yaxes(title_text="Initial state", row=2, col=1)
fig.update_yaxes(title_text="ΔQ (amu^0.5·Å)", row=2, col=2)

fig.show()

# Save for publication
fig.write_image('figure_complete_analysis.pdf', scale=2)
```

**Features**:
- Four panels in one figure
- Consistent styling
- Publication-ready export
- Tells complete story

---

## Interactive Dashboard

### Web-Based Exploration

Real-time parameter tuning:

```bash
# Launch dashboard
carriercapture viz --port 8050

# Or from Python
from carriercapture.visualization import run_server
run_server(port=8050)
```

**Features**:
- Upload Q-E data
- Adjust fitting parameters (order, smoothness)
- Solve Schrödinger equation
- Calculate capture interactively
- Export results and figures
- All in web browser

---

## Code Examples

### Custom Styling

Apply consistent theme:

```python
# Define custom colors
COLORS = {
    'primary': '#1f77b4',
    'secondary': '#ff7f0e',
    'tertiary': '#2ca02c',
    'background': '#ffffff',
    'grid': '#e5e5e5'
}

def apply_custom_theme(fig):
    """Apply custom theme to figure."""
    fig.update_layout(
        font=dict(family='Arial, sans-serif', size=12),
        plot_bgcolor=COLORS['background'],
        paper_bgcolor=COLORS['background'],
        xaxis=dict(gridcolor=COLORS['grid']),
        yaxis=dict(gridcolor=COLORS['grid']),
        template='plotly_white'
    )
    return fig

# Use
fig = plot_potential(pot)
fig = apply_custom_theme(fig)
fig.show()
```

### Batch Processing

Generate all figures for a calculation:

```python
def create_figure_set(pot_i, pot_f, cc, output_dir='figures'):
    """Generate complete figure set."""
    import os
    os.makedirs(output_dir, exist_ok=True)

    # Figure 1: Potentials
    fig1 = plot_potential(pot_i, show_wavefunctions=True, title="Initial State")
    fig1.write_html(f'{output_dir}/01_potential_initial.html')

    # Figure 2: CC diagram
    fig2 = plot_configuration_coordinate(pot_i, pot_f, Q0=10.0, show_crossing=True)
    fig2.write_html(f'{output_dir}/02_cc_diagram.html')

    # Figure 3: Capture
    fig3 = plot_capture_coefficient(cc, arrhenius=True)
    fig3.write_html(f'{output_dir}/03_capture_arrhenius.html')

    # Figure 4: Overlaps
    fig4 = plot_overlap_matrix(cc, log_scale=True)
    fig4.write_html(f'{output_dir}/04_overlap_matrix.html')

    print(f"✓ Created 4 figures in {output_dir}/")

# Use
create_figure_set(pot_initial, pot_final, cc)
```

---

## Tips for Great Visualizations

### 1. Choose Appropriate Plot Type

- **Line plots**: Temperature dependence, 1D scans
- **Heatmaps**: 2D parameter scans, overlap matrices
- **Scatter plots**: Discrete data points, comparisons
- **Contour plots**: Iso-value lines, boundaries

### 2. Use Log Scales Wisely

```python
# Capture coefficients span many orders of magnitude
# → Always use log scale
fig = plot_capture_coefficient(cc, arrhenius=True)  # Already log scale

# Overlap matrices have wide dynamic range
# → Log scale for better visibility
fig = plot_overlap_matrix(cc, log_scale=True)

# Parameter scans
# → Log scale to see full range
fig = plot_scan_heatmap(results, log_scale=True)
```

### 3. Add Annotations

```python
# Annotate key points
fig.add_annotation(
    x=x_point, y=y_point,
    text="Maximum capture",
    showarrow=True,
    arrowhead=2,
    arrowcolor='red'
)

# Add text box
fig.add_annotation(
    xref='paper', yref='paper',
    x=0.95, y=0.95,
    text="T = 300K<br>W = 0.205 eV",
    showarrow=False,
    bgcolor='white',
    bordercolor='black',
    borderwidth=1
)
```

### 4. Consistent Color Schemes

- **Viridis**: Good default, colorblind-friendly
- **RdBu**: Diverging data (positive/negative)
- **Plasma**: High contrast
- **Greens/Blues**: Single-hue gradients

### 5. Export High Resolution

```python
# For publications
fig.write_image('figure.pdf', scale=2)  # Vector, 2× resolution
fig.write_image('figure.png', width=1600, height=1200, scale=3)  # Raster, high DPI

# For presentations
fig.write_image('slide.png', width=1920, height=1080)  # 1080p

# For web
fig.write_html('interactive.html', include_plotlyjs='cdn')  # Smaller file size
```

---

## See Also

- **[User Guide: Visualization](../user-guide/visualization.md)** - Comprehensive plotting guide
- **[API Reference: Visualization](../api/visualization.md)** - Complete API
- **[Examples: Notebooks](notebooks.md)** - Jupyter notebook examples
- **[Example Notebook: Sn in ZnO](https://github.com/WMD-group/CarrierCapture.py/blob/main/examples/notebooks/01_harmonic_sn_zn.ipynb)** - Complete example
