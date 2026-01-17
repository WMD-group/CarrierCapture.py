# CarrierCapture.py

**Modern Python package for computing carrier capture rates and non-radiative recombination in semiconductors using multiphonon theory.**

[![Tests](https://github.com/WMD-group/CarrierCapture.py/actions/workflows/tests.yml/badge.svg)](https://github.com/WMD-group/CarrierCapture.py/actions/workflows/tests.yml)
[![Python Version](https://img.shields.io/badge/python-3.11%20%7C%203.12-blue)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

---

## ✨ Features

### 🧮 Complete Multiphonon Theory Implementation
- 1D Schrödinger equation solver (ARPACK-based)
- Harmonic, Morse, and spline potential fitting
- Configuration coordinate diagrams
- Capture coefficient calculations

### ⚡ High-Performance Computing
- Parallel parameter scanning with `joblib`
- Optimized with NumPy/SciPy (within 10-20% of Julia speed)
- Support for HDF5 and NPZ result storage

### 🎨 Rich Visualization
- Publication-quality plots with Plotly
- Interactive Dash dashboard (web-based)
- Real-time parameter exploration
- Arrhenius plots, CC diagrams, 2D heatmaps

### 🔗 Integration with [doped](https://github.com/SMTG-Bham/doped)
- Direct compatibility with doped defect calculation workflows
- Load configuration coordinate data from doped outputs
- Automated mass-weighted displacement calculations
- Seamless end-to-end: DFT → defect analysis → capture rates

---

## 🚀 Quick Start

### Installation

```bash
pip install carriercapture
```

For development or from source:

```bash
git clone https://github.com/WMD-group/CarrierCapture.py.git
cd CarrierCapture.py
pip install -e ".[dev]"
```

### Python API

```python
import numpy as np
from carriercapture.core.potential import Potential
from carriercapture.core.config_coord import ConfigCoordinate

# Create harmonic potentials
pot_initial = Potential.from_harmonic(hw=0.008, Q0=0.0, E0=0.5)
pot_final = Potential.from_harmonic(hw=0.008, Q0=10.5, E0=0.0)

# Solve Schrödinger equation
pot_initial.solve(nev=180)
pot_final.solve(nev=60)

# Calculate capture coefficient
cc = ConfigCoordinate(pot_i=pot_initial, pot_f=pot_final, W=0.068)
cc.calculate_overlap(Q0=5.0)
cc.calculate_capture_coefficient(
    volume=1e-21,
    temperature=np.linspace(100, 500, 50)
)

print(f"Capture coefficient at 300K: {cc.capture_coefficient[20]:.3e} cm³/s")
```

### Command-Line Interface

```bash
# Fit potential from DFT data
carriercapture fit data/excited.dat -f spline -o 4 -s 0.001 -O excited.json

# Solve Schrödinger equation
carriercapture solve excited.json -n 180 -O excited_solved.json

# Calculate capture coefficient
carriercapture capture config.yaml -V 1e-21 --temp-range 100 500 50

# High-throughput parameter scan
carriercapture scan --dQ-min 0 --dQ-max 25 --dQ-points 25 \
                    --dE-min 0 --dE-max 2.5 --dE-points 10 \
                    -j -1 -o scan_results.npz

# Launch interactive dashboard
carriercapture viz --port 8050
```

---

## 📚 Documentation Sections

### [Getting Started](getting-started/installation.md)
New to CarrierCapture? Start here!

- **[Installation](getting-started/installation.md)** - Setup instructions
- **[Quick Start](getting-started/quick-start.md)** - 5-minute introduction
- **[Basic Concepts](getting-started/basic-concepts.md)** - Key concepts explained
- **[First Calculation](getting-started/first-calculation.md)** - Step-by-step tutorial

### [User Guide](user-guide/potentials.md)
Learn how to use CarrierCapture effectively.

- **[Potentials](user-guide/potentials.md)** - Creating and fitting potentials
- **[Capture Coefficients](user-guide/capture-coefficients.md)** - Calculating capture rates
- **[Parameter Scanning](user-guide/parameter-scanning.md)** - High-throughput screening
- **[Visualization](user-guide/visualization.md)** - Plots and dashboards
- **[doped Integration](user-guide/doped-integration.md)** - Using with doped
- **[CLI Usage](user-guide/cli-usage.md)** - Command-line interface
- **[Configuration](user-guide/configuration.md)** - YAML configuration files

### [Theory](theory/multiphonon-theory.md)
Understand the science behind CarrierCapture.

- **[Multiphonon Theory](theory/multiphonon-theory.md)** - Huang-Rhys formalism
- **[Configuration Coordinates](theory/configuration-coordinates.md)** - CC diagrams
- **[Equations](theory/equations.md)** - Mathematical formulation
- **[References](theory/references.md)** - Literature and citations

### [Tutorials](tutorials/index.md)
Hands-on examples with Jupyter notebooks.

- **[1. Harmonic Oscillator](tutorials/01-harmonic-oscillator.md)** - Basic workflow
- **[2. DX Center](tutorials/02-dx-center.md)** - Real-world anharmonic example
- **[3. Parameter Scan](tutorials/03-parameter-scan.md)** - High-throughput screening
- **[4. Interactive Dashboard](tutorials/04-interactive-dashboard.md)** - Web visualization

### [API Reference](api/index.md)
Complete API documentation.

- **[Core](api/core.md)** - Potential, ConfigCoordinate, Schrödinger solver
- **[Analysis](api/analysis.md)** - Parameter scanning
- **[I/O](api/io.md)** - File I/O and doped interface
- **[Visualization](api/visualization.md)** - Plotting functions
- **[CLI](api/cli.md)** - Command-line interface

### [Development](development/contributing.md)
Contributing to CarrierCapture.

- **[Contributing](development/contributing.md)** - How to contribute
- **[Testing](development/testing.md)** - Running tests
- **[Architecture](development/architecture.md)** - Code structure

---

## 🔬 Scientific Background

CarrierCapture implements the **static coupling approximation** for non-radiative carrier capture via multiphonon emission (Huang-Rhys theory).

### Key Equation

The capture coefficient is calculated as:

$$C(T) = \frac{V \cdot 2\pi}{\hbar} \cdot g \cdot W^2 \cdot \sum_{i,j} p_i |\langle\chi_i|Q-Q_0|\chi_j\rangle|^2 \delta(\varepsilon_i - \varepsilon_j)$$

Where:

- $V$: supercell volume
- $g$: degeneracy factor
- $W$: electron-phonon coupling matrix element
- $p_i$: thermal occupation of initial state $i$
- $\chi_i, \chi_j$: vibrational wavefunctions
- $\delta$: energy-conserving delta function (Gaussian broadened)

### References

1. **Alkauskas et al.** (2014) - [Phys. Rev. B **90**, 075202](https://doi.org/10.1103/PhysRevB.90.075202)
2. **Huang & Rhys** (1950) - Proc. R. Soc. Lond. A **204**, 406

---

## 🎯 Project Status

| Component | Status |
|-----------|--------|
| Core Engine | ✅ Complete |
| CLI | ✅ Complete |
| Visualization | ✅ Complete |
| Parameter Scanning | ✅ Complete |
| doped Integration | ✅ Complete |
| Documentation | ✅ Complete |
| Test Coverage | ✅ 88 tests |
| PyPI Release | 🔄 Planned |

---

## 🤝 Contributing

Contributions are welcome! See our [Contributing Guide](development/contributing.md) for details.

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/WMD-group/CarrierCapture.py/blob/main/LICENSE) file for details.

---

## 🙏 Acknowledgments

- Original Julia implementation: [WMD-group/CarrierCapture.jl](https://github.com/WMD-group/CarrierCapture.jl)
- Built with: NumPy, SciPy, Plotly, Dash, Click
- Compatible with: [doped](https://github.com/SMTG-Bham/doped) (defect calculations)
