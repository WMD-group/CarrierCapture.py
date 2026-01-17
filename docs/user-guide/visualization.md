# Visualization

Publication-quality plots and interactive dashboards for exploring carrier capture calculations.

## Overview

CarrierCapture provides two complementary visualization approaches:

1. **Static Plots** (Plotly) - Publication-quality figures for papers and presentations
2. **Interactive Dashboard** (Dash) - Web-based exploration and parameter tuning

All visualizations are built on Plotly for consistent styling and interactivity.

---

## Static Plots

### Installation

Static plotting requires only the base package:

```bash
pip install carriercapture  # Plotly included by default
```

### Quick Start

```python
from carriercapture.visualization import plot_potential
from carriercapture.core import Potential
from carriercapture.io import load_potential_from_file

# Load solved potential
pot = Potential.from_dict(load_potential_from_file('potential_solved.json'))

# Plot with wavefunctions
fig = plot_potential(
    pot,
    title="Potential Energy Surface",
    show_wavefunctions=True,
    n_states=10
)

# Display in browser
fig.show()

# Save to file
fig.write_html('potential.html')
fig.write_image('potential.png', width=800, height=600)
```

---

## Available Plot Types

### 1. Potential Energy Surface

```python
from carriercapture.visualization import plot_potential

# Basic plot
fig = plot_potential(pot)
fig.show()

# With wavefunctions
fig = plot_potential(
    pot,
    show_wavefunctions=True,
    n_states=10,              # Number of states to show
    show_eigenvalues=True,    # Show horizontal lines at eigenvalues
    wavefunction_scale=0.3,   # Scaling factor for visibility
    title="Ground State Potential"
)
fig.show()

# Customize appearance
fig = plot_potential(
    pot,
    show_wavefunctions=True,
    line_color='blue',
    wavefunction_colors='rainbow',  # or specific color
    show_grid=True
)
fig.update_layout(
    width=900,
    height=600,
    font=dict(size=14)
)
fig.show()
```

**Use cases:**
- Visualize fitted potential quality
- Display eigenvalue spectrum
- Show wavefunction character
- Compare different potentials

### 2. Configuration Coordinate Diagram

```python
from carriercapture.visualization import plot_configuration_coordinate

# Basic CC diagram
fig = plot_configuration_coordinate(
    pot_initial,
    pot_final,
    Q0=10.0,                  # Configuration coordinate shift
    show_crossing=True,       # Show potential crossing
    show_wavefunctions=True,  # Show phonon wavefunctions
    n_states=10,              # States per potential
    title="Configuration Coordinate Diagram"
)
fig.show()

# Highlight specific transitions
fig = plot_configuration_coordinate(
    pot_initial,
    pot_final,
    Q0=10.0,
    show_crossing=True,
    show_wavefunctions=True,
    highlight_transitions=[(0, 0), (1, 1), (2, 2)],  # (initial, final) pairs
    transition_color='red',
    transition_width=2
)
fig.show()
```

**Use cases:**
- Illustrate carrier capture mechanism
- Show Franck-Condon transitions
- Display potential crossing point
- Explain multiphonon emission

### 3. Eigenvalue Spectrum

```python
from carriercapture.visualization import plot_eigenvalue_spectrum

# Simple spectrum
fig = plot_eigenvalue_spectrum(
    pot,
    max_states=20,
    show_spacing=True,  # Annotate level spacing
    title="Vibrational Energy Levels"
)
fig.show()

# Compare multiple potentials
fig = go.Figure()
fig.add_trace(plot_eigenvalue_spectrum(pot_initial, name="Initial"))
fig.add_trace(plot_eigenvalue_spectrum(pot_final, name="Final"))
fig.update_layout(title="Energy Level Comparison")
fig.show()
```

**Use cases:**
- Check phonon level spacing
- Verify harmonic vs. anharmonic character
- Compare different charge states

### 4. Capture Coefficient vs Temperature

```python
from carriercapture.visualization import plot_capture_coefficient

# Standard Arrhenius plot
fig = plot_capture_coefficient(
    cc,  # ConfigCoordinate object after calculation
    arrhenius=True,              # Log scale, 1000/T x-axis
    show_activation_energy=False, # Fit and show E_a
    title="Capture Coefficient vs Temperature"
)
fig.show()

# Linear temperature scale
fig = plot_capture_coefficient(
    cc,
    arrhenius=False,  # T on x-axis
    log_scale=True,   # Log C on y-axis
    title="C(T) Linear Scale"
)
fig.show()

# With experimental data overlay
fig = plot_capture_coefficient(cc, arrhenius=True)

# Add experimental data
import plotly.graph_objects as go
T_exp = [200, 250, 300, 350, 400]
C_exp = [1e-12, 5e-12, 2e-11, 7e-11, 2e-10]
fig.add_trace(go.Scatter(
    x=[1000/T for T in T_exp],
    y=[np.log10(C) for C in C_exp],
    mode='markers',
    name='Experiment',
    marker=dict(size=10, symbol='diamond')
))
fig.show()
```

**Use cases:**
- Present capture coefficient results
- Extract activation energy
- Compare with experiments
- Analyze temperature dependence

### 5. Overlap Matrix Heatmap

```python
from carriercapture.visualization import plot_overlap_matrix

# Log scale (recommended for overlaps)
fig = plot_overlap_matrix(
    cc,  # ConfigCoordinate object after overlap calculation
    log_scale=True,
    title="Franck-Condon Overlap Matrix",
    colorscale='RdBu'
)
fig.show()

# Linear scale
fig = plot_overlap_matrix(
    cc,
    log_scale=False,
    show_values=True,  # Annotate values (for small matrices)
    threshold=1e-6     # Only show values above threshold
)
fig.show()
```

**Use cases:**
- Identify dominant transitions
- Check energy-conserving delta function
- Diagnose overlap calculation issues

### 6. Parameter Scan Heatmap

```python
from carriercapture.visualization import plot_scan_heatmap

# Basic heatmap
fig = plot_scan_heatmap(
    scan_results,  # ScanResult object
    log_scale=True,
    title="Parameter Scan: Capture Coefficient",
    colorscale='Viridis'
)
fig.show()

# Contour plot
fig = plot_scan_heatmap(
    scan_results,
    plot_type='contour',  # 'heatmap', 'contour', or 'both'
    log_scale=True,
    contour_levels=10
)
fig.show()

# Annotated with maximum
fig = plot_scan_heatmap(scan_results, log_scale=True)

# Find and annotate maximum
C = scan_results.capture_coefficients
i_max, j_max = np.unravel_index(np.nanargmax(C), C.shape)
dQ_max = scan_results.dQ_grid[i_max]
dE_max = scan_results.dE_grid[j_max]

fig.add_annotation(
    x=dE_max, y=dQ_max,
    text=f"Max: C={C[i_max,j_max]:.2e}",
    showarrow=True,
    arrowhead=2
)
fig.show()
```

**Use cases:**
- Visualize parameter scan results
- Identify optimal parameters
- Show design principles
- Materials screening

---

## Customization

### Themes and Styling

```python
from carriercapture.visualization import COLORS, get_default_layout

# Access color scheme
colors = COLORS
primary_color = colors['primary']
secondary_color = colors['secondary']

# Get default layout
layout = get_default_layout(
    title="My Custom Plot",
    xaxis_title="Configuration Coordinate (amu^0.5·Å)",
    yaxis_title="Energy (eV)"
)

# Create custom figure
import plotly.graph_objects as go
fig = go.Figure(layout=layout)

# Add your traces
fig.add_trace(go.Scatter(
    x=pot.Q,
    y=pot.E,
    mode='lines',
    line=dict(color=primary_color, width=2)
))

fig.show()
```

### Publication Style

```python
from carriercapture.visualization import apply_publication_style

# Create plot
fig = plot_potential(pot, show_wavefunctions=True)

# Apply publication styling
fig = apply_publication_style(
    fig,
    width=800,
    height=600,
    font_size=14,
    line_width=2.5,
    show_legend=True
)

# Export high-resolution
fig.write_image('figure_1.pdf', scale=2)  # 2× resolution
fig.write_image('figure_1.png', scale=3)  # 3× resolution
```

### Custom Color Scales

```python
# Custom colorscale for heatmaps
custom_colorscale = [
    [0.0, 'rgb(255,255,255)'],   # White
    [0.5, 'rgb(100,149,237)'],   # Cornflower blue
    [1.0, 'rgb(25,25,112)']      # Midnight blue
]

fig = plot_scan_heatmap(
    scan_results,
    log_scale=True,
    colorscale=custom_colorscale
)
fig.show()

# Named Plotly colorscales:
# - 'Viridis', 'Plasma', 'Inferno', 'Magma'
# - 'RdBu', 'RdYlBu', 'Spectral'
# - 'Blues', 'Greens', 'Reds'
```

### Layout Customization

```python
fig = plot_potential(pot)

# Update layout
fig.update_layout(
    title=dict(
        text="Potential Energy Surface",
        font=dict(size=18, family='Arial, sans-serif')
    ),
    xaxis=dict(
        title="Q (amu<sup>0.5</sup>·Å)",
        title_font=dict(size=14),
        tickfont=dict(size=12),
        showgrid=True,
        gridcolor='lightgray'
    ),
    yaxis=dict(
        title="E (eV)",
        title_font=dict(size=14),
        range=[0, 2.0]  # Set axis range
    ),
    width=900,
    height=600,
    template='plotly_white',  # or 'plotly', 'plotly_dark', 'seaborn', etc.
    showlegend=True,
    legend=dict(
        x=0.02, y=0.98,
        bgcolor='rgba(255,255,255,0.8)',
        bordercolor='black',
        borderwidth=1
    )
)

fig.show()
```

---

## Exporting Figures

### HTML (Interactive)

```python
# Save interactive HTML
fig.write_html(
    'figure.html',
    include_plotlyjs='cdn',  # Use CDN for smaller file size
    # include_plotlyjs=True,  # Embed Plotly.js (larger, offline-ready)
)

# With custom configuration
fig.write_html(
    'figure.html',
    config={
        'displayModeBar': True,
        'displaylogo': False,
        'modeBarButtonsToRemove': ['pan2d', 'lasso2d']
    }
)
```

### Static Images

Requires `kaleido` package:

```bash
pip install kaleido
```

```python
# PNG (raster)
fig.write_image('figure.png', width=1200, height=800, scale=2)

# PDF (vector)
fig.write_image('figure.pdf', width=8, height=6)  # inches

# SVG (vector)
fig.write_image('figure.svg')

# JPEG
fig.write_image('figure.jpg', width=1600, height=1200)
```

### Multiple Formats

```python
# Save all formats
for fmt in ['html', 'png', 'pdf', 'svg']:
    if fmt == 'html':
        fig.write_html(f'figure.{fmt}')
    else:
        fig.write_image(f'figure.{fmt}', width=1200, height=800, scale=2)

print("✓ Saved to figure.html, figure.png, figure.pdf, figure.svg")
```

---

## Interactive Dashboard

### Installation

Dashboard requires additional dependencies:

```bash
pip install carriercapture[viz]
# Installs: dash, dash-bootstrap-components
```

### Launching the Dashboard

**From Python:**

```python
from carriercapture.visualization import create_app, run_server

# Create Dash app
app = create_app()

# Run server
run_server(
    port=8050,
    debug=False,
    host='127.0.0.1'
)

# Browser opens automatically to http://127.0.0.1:8050
```

**From Command Line:**

```bash
# Launch dashboard
carriercapture viz

# Custom port
carriercapture viz --port 8080

# With data preloaded
carriercapture viz --data potential.json

# Debug mode (for development)
carriercapture viz --debug

# Don't open browser automatically
carriercapture viz --no-browser
```

### Dashboard Features

The interactive dashboard provides:

1. **Data Loading**
   - Load potential data from files
   - Import Q-E data from CSV
   - Create test potentials

2. **Potential Fitting**
   - Interactive parameter adjustment
   - Real-time fit visualization
   - Multiple fitting methods
   - Fit quality metrics

3. **Schrödinger Solver**
   - Adjust number of eigenvalues
   - Visualize wavefunctions
   - Inspect eigenvalue spectrum

4. **Configuration Coordinate**
   - Two-potential CC diagram
   - Adjust Q₀ interactively
   - Show transitions

5. **Capture Calculation**
   - Set all parameters via UI
   - Calculate C(T) on the fly
   - Arrhenius plot
   - Export results

6. **Parameter Scanning**
   - Define scan ranges
   - Monitor progress
   - Visualize heatmap
   - Export scan results

### Dashboard Workflow

```python
# 1. Launch dashboard
from carriercapture.visualization import run_server
run_server(port=8050)

# 2. In browser (http://127.0.0.1:8050):
#    - Upload Q-E data or load example
#    - Fit potential (adjust parameters)
#    - Solve for eigenvalues
#    - Repeat for second potential
#    - Set up capture calculation
#    - View results
#    - Export figures and data

# 3. Results saved to downloads folder
```

### Programmatic Dashboard

```python
from carriercapture.visualization import create_app
from carriercapture.core import Potential

# Create app with pre-loaded data
pot_i = Potential.from_dict(load_potential_from_file('excited.json'))
pot_f = Potential.from_dict(load_potential_from_file('ground.json'))

app = create_app(
    pot_initial=pot_i,
    pot_final=pot_f,
    default_params={
        'W': 0.205,
        'Q0': 10.0,
        'volume': 1e-21,
        'degeneracy': 1
    }
)

# Run
app.run_server(debug=False, port=8050)
```

---

## Gallery Examples

### Example 1: CC Diagram with Transitions

```python
from carriercapture.visualization import plot_configuration_coordinate
import plotly.graph_objects as go

# Create CC diagram
fig = plot_configuration_coordinate(
    pot_initial,
    pot_final,
    Q0=10.0,
    show_crossing=True,
    show_wavefunctions=True,
    n_states=8
)

# Highlight specific transition (e.g., 0→3)
Q = pot_initial.Q
psi_i = pot_initial.eigenvectors[0, :]
psi_f = pot_final.eigenvectors[3, :]

# Add transition arrow
fig.add_annotation(
    x=5, y=pot_initial.eigenvalues[0],
    ax=5 + 10.0, ay=pot_final.eigenvalues[3],
    xref='x', yref='y',
    axref='x', ayref='y',
    showarrow=True,
    arrowhead=2,
    arrowsize=1,
    arrowwidth=2,
    arrowcolor='red'
)

fig.show()
```

### Example 2: Multi-Panel Figure

```python
from plotly.subplots import make_subplots
import plotly.graph_objects as go

# Create 2×2 subplot layout
fig = make_subplots(
    rows=2, cols=2,
    subplot_titles=(
        'Initial State', 'Final State',
        'Overlap Matrix', 'Capture Coefficient'
    )
)

# Panel 1: Initial potential
fig.add_trace(
    go.Scatter(x=pot_i.Q, y=pot_i.E, mode='lines', name='Initial'),
    row=1, col=1
)

# Panel 2: Final potential
fig.add_trace(
    go.Scatter(x=pot_f.Q, y=pot_f.E, mode='lines', name='Final'),
    row=1, col=2
)

# Panel 3: Overlap matrix
fig.add_trace(
    go.Heatmap(z=np.log10(np.abs(cc.overlap_matrix)+1e-30),
               colorscale='Viridis'),
    row=2, col=1
)

# Panel 4: Capture coefficient
fig.add_trace(
    go.Scatter(x=1000/cc.temperature, y=np.log10(cc.capture_coefficient),
               mode='markers+lines', name='C(T)'),
    row=2, col=2
)

# Update layout
fig.update_layout(height=800, width=1200, showlegend=False)
fig.update_xaxes(title_text="Q (amu^0.5·Å)", row=1, col=1)
fig.update_xaxes(title_text="Q (amu^0.5·Å)", row=1, col=2)
fig.update_xaxes(title_text="Final state", row=2, col=1)
fig.update_xaxes(title_text="1000/T (K^-1)", row=2, col=2)
fig.update_yaxes(title_text="E (eV)", row=1, col=1)
fig.update_yaxes(title_text="E (eV)", row=1, col=2)
fig.update_yaxes(title_text="Initial state", row=2, col=1)
fig.update_yaxes(title_text="log₁₀(C)", row=2, col=2)

fig.show()
```

### Example 3: Animated Temperature Scan

```python
import numpy as np
import plotly.graph_objects as go

# Calculate C(T) for multiple parameter sets
temperatures = np.linspace(100, 500, 50)
Q0_values = [8, 10, 12, 14]
results = {}

for Q0 in Q0_values:
    cc_temp = ConfigCoordinate(pot_i, pot_f, W=0.205, degeneracy=1)
    cc_temp.calculate_overlap(Q0=Q0, cutoff=0.25, sigma=0.025)
    cc_temp.calculate_capture_coefficient(volume=1e-21, temperature=temperatures)
    results[Q0] = cc_temp.capture_coefficient

# Create animated figure
fig = go.Figure()

# Add trace for each Q0
for Q0, C in results.items():
    fig.add_trace(go.Scatter(
        x=1000/temperatures,
        y=np.log10(C),
        mode='lines+markers',
        name=f'Q₀={Q0} amu^0.5·Å',
        line=dict(width=2)
    ))

fig.update_layout(
    title="Capture Coefficient: Q₀ Dependence",
    xaxis_title="1000/T (K⁻¹)",
    yaxis_title="log₁₀(C) [cm³/s]",
    hovermode='x unified'
)

fig.show()
```

---

## Best Practices

### 1. Consistent Styling

```python
# Define consistent theme for all figures
PLOT_CONFIG = {
    'width': 900,
    'height': 600,
    'font_size': 14,
    'line_width': 2,
    'template': 'plotly_white'
}

def style_figure(fig):
    """Apply consistent styling."""
    fig.update_layout(
        width=PLOT_CONFIG['width'],
        height=PLOT_CONFIG['height'],
        font=dict(size=PLOT_CONFIG['font_size']),
        template=PLOT_CONFIG['template']
    )
    return fig

# Use for all plots
fig1 = plot_potential(pot)
fig1 = style_figure(fig1)

fig2 = plot_capture_coefficient(cc)
fig2 = style_figure(fig2)
```

### 2. Publication-Ready Exports

```python
# High-resolution export for publication
def export_publication_figure(fig, basename):
    """Export figure in multiple formats for publication."""
    # Vector formats (preferred for publications)
    fig.write_image(f'{basename}.pdf', width=1200, height=800)
    fig.write_image(f'{basename}.svg')

    # Raster format (for previews)
    fig.write_image(f'{basename}.png', width=1600, height=1200, scale=2)

    # Interactive version
    fig.write_html(f'{basename}.html')

    print(f"✓ Exported {basename} in PDF, SVG, PNG, HTML")

# Use
fig = plot_configuration_coordinate(pot_i, pot_f, Q0=10.0)
fig = apply_publication_style(fig)
export_publication_figure(fig, 'figure_2_cc_diagram')
```

### 3. Batch Visualization

```python
# Create all figures for a calculation
def create_all_figures(pot_i, pot_f, cc, output_dir='figures'):
    """Generate complete figure set."""
    import os
    os.makedirs(output_dir, exist_ok=True)

    # 1. Initial potential
    fig1 = plot_potential(pot_i, show_wavefunctions=True, title="Initial State")
    fig1.write_html(f'{output_dir}/01_potential_initial.html')

    # 2. Final potential
    fig2 = plot_potential(pot_f, show_wavefunctions=True, title="Final State")
    fig2.write_html(f'{output_dir}/02_potential_final.html')

    # 3. CC diagram
    fig3 = plot_configuration_coordinate(pot_i, pot_f, Q0=cc.Q0, show_crossing=True)
    fig3.write_html(f'{output_dir}/03_cc_diagram.html')

    # 4. Overlap matrix
    fig4 = plot_overlap_matrix(cc, log_scale=True)
    fig4.write_html(f'{output_dir}/04_overlap_matrix.html')

    # 5. Capture coefficient
    fig5 = plot_capture_coefficient(cc, arrhenius=True)
    fig5.write_html(f'{output_dir}/05_capture_arrhenius.html')

    print(f"✓ Created 5 figures in {output_dir}/")

# Use
create_all_figures(pot_initial, pot_final, cc, output_dir='my_calculation_figures')
```

---

## Troubleshooting

### Dashboard won't start

```python
# Check if port is already in use
import socket

def check_port(port):
    """Check if port is available."""
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    result = sock.connect_ex(('127.0.0.1', port))
    sock.close()
    return result != 0  # True if available

if not check_port(8050):
    print("Port 8050 is busy. Try a different port:")
    print("  carriercapture viz --port 8051")
else:
    print("Port 8050 is available")
```

### Image export fails

```python
# Install kaleido if missing
try:
    import kaleido
    print("✓ Kaleido installed")
except ImportError:
    print("⚠️  Kaleido not found. Install with:")
    print("  pip install kaleido")

# Alternative: Use orca (deprecated but sometimes works)
# pip install plotly-orca
```

### Figure not displaying in Jupyter

```python
# Use Plotly's Jupyter renderer
import plotly.io as pio
pio.renderers.default = 'notebook'  # or 'jupyterlab', 'colab'

fig.show()
```

---

## See Also

- **[API Reference: Visualization](../api/visualization.md)** - Complete API documentation
- **[CLI Reference: viz](../api/cli.md#viz---interactive-dashboard)** - Dashboard CLI options
- **[CLI Reference: plot](../api/cli.md#plot---static-plots)** - Static plotting CLI
- **[Tutorial 4: Dashboard](../tutorials/04-interactive-dashboard.md)** - Interactive tutorial
- **[Examples: Gallery](../examples/gallery.md)** - Gallery of example plots
