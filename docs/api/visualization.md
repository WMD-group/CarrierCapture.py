# Visualization API

Plotting functions and interactive dashboards for visualizing results.

## Static Plots

Publication-quality plots using Plotly.

::: carriercapture.visualization.static
    options:
      show_root_heading: true
      show_source: false
      heading_level: 3

## Interactive Dashboard

Web-based interactive dashboard using Dash.

!!! info "Optional Dependency"
    The dashboard requires `pip install carriercapture[viz]`

::: carriercapture.visualization.interactive
    options:
      show_root_heading: true
      show_source: false
      heading_level: 3

## Themes

Styling and color schemes for plots.

::: carriercapture.visualization.themes
    options:
      show_root_heading: true
      show_source: false
      heading_level: 3

## Usage Examples

### Plot Potential Energy Surface

```python
from carriercapture.visualization import plot_potential
from carriercapture.core import Potential

# Create and solve potential
pot = Potential.from_harmonic(hw=0.008, Q0=0.0, E0=0.0)
pot.solve(nev=60)

# Plot with wavefunctions
fig = plot_potential(
    pot,
    title="Harmonic Potential",
    show_wavefunctions=True,
    n_states=10,  # Show first 10 states
    show_eigenvalues=True
)

# Display
fig.show()

# Save
fig.write_html('potential.html')
fig.write_image('potential.png', width=800, height=600)
```

### Plot Capture Coefficient

```python
from carriercapture.visualization import plot_capture_coefficient

# After calculating capture coefficient
cc = ConfigCoordinate(...)
cc.calculate_capture_coefficient(...)

# Arrhenius plot
fig = plot_capture_coefficient(
    cc,
    title="Capture Coefficient vs Temperature",
    arrhenius=True,  # Log scale, 1000/T x-axis
    show_activation_energy=True
)

fig.show()
```

### Configuration Coordinate Diagram

```python
from carriercapture.visualization import plot_configuration_coordinate

# Plot CC diagram with both potentials
fig = plot_configuration_coordinate(
    pot_initial,
    pot_final,
    show_crossing=True,
    show_wavefunctions=True
)

fig.show()
```

### Eigenvalue Spectrum

```python
from carriercapture.visualization import plot_eigenvalue_spectrum

fig = plot_eigenvalue_spectrum(
    pot,
    title="Vibrational Energy Levels",
    max_states=20
)

fig.show()
```

### Overlap Matrix Heatmap

```python
from carriercapture.visualization import plot_overlap_matrix

# After calculating overlaps
cc = ConfigCoordinate(...)
cc.calculate_overlap(...)

fig = plot_overlap_matrix(
    cc,
    title="Franck-Condon Overlaps",
    log_scale=True  # Use log scale for small values
)

fig.show()
```

### Parameter Scan 2D Heatmap

```python
from carriercapture.visualization import plot_scan_heatmap
from carriercapture.analysis import ParameterScanner

# After running scan
results = scanner.run_harmonic_scan(...)

fig = plot_scan_heatmap(
    results,
    title="Capture Coefficient Scan",
    log_scale=True,
    colorscale='Viridis'
)

fig.show()
```

### Interactive Dashboard

```python
from carriercapture.visualization import create_app

# Create Dash app
app = create_app(
    pot_initial=pot_i,
    pot_final=pot_f,
    cc=cc
)

# Run server
app.run_server(debug=False, port=8050)
```

Or from command line:

```bash
carriercapture viz --port 8050
```

## Plot Customization

### Themes

```python
from carriercapture.visualization import COLORS, get_default_layout

# Access color scheme
colors = COLORS
primary = colors['primary']
secondary = colors['secondary']

# Get default layout
layout = get_default_layout(
    title="My Plot",
    xaxis_title="Configuration Coordinate (amu^0.5·Å)",
    yaxis_title="Energy (eV)"
)

# Create custom figure
fig = go.Figure(layout=layout)
fig.add_trace(...)
```

### Publication Style

```python
from carriercapture.visualization import apply_publication_style

fig = plot_potential(pot)

# Apply publication styling
fig = apply_publication_style(
    fig,
    width=800,
    height=600,
    font_size=14,
    line_width=2
)

# Export high-res
fig.write_image('figure.pdf', scale=2)
```

## Plot Types

| Function | Purpose | Key Options |
|----------|---------|-------------|
| `plot_potential` | PES with wavefunctions | `show_wavefunctions`, `n_states` |
| `plot_capture_coefficient` | C(T) vs temperature | `arrhenius`, `show_activation_energy` |
| `plot_configuration_coordinate` | CC diagram | `show_crossing`, `show_wavefunctions` |
| `plot_eigenvalue_spectrum` | Energy levels | `max_states`, `show_spacing` |
| `plot_overlap_matrix` | Franck-Condon factors | `log_scale` |
| `plot_scan_heatmap` | 2D parameter scan | `log_scale`, `colorscale` |

## Interactive Features

All Plotly plots support:

- **Zoom**: Click and drag to zoom
- **Pan**: Shift + drag to pan
- **Hover**: Hover for data values
- **Export**: Camera icon to save as PNG
- **Reset**: Home icon to reset axes

## See Also

- **[User Guide: Visualization](../user-guide/visualization.md)** - Comprehensive visualization guide
- **[Tutorial 4: Dashboard](../tutorials/04-interactive-dashboard.md)** - Interactive tutorial
- **[Examples: Gallery](../examples/gallery.md)** - Gallery of example plots
