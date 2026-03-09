# CarrierCapture.py Examples

This directory contains tutorial notebooks and example data for CarrierCapture.py.

## Tutorial Notebooks

### 01. Basic Harmonic Workflow (`notebooks/01_harmonic_sn_zn.ipynb`)
Introduction to CarrierCapture with a simple harmonic oscillator model. Covers:
- Creating harmonic potentials
- Solving the Schrödinger equation
- Calculating capture coefficients
- Visualizing results

**Difficulty**: Beginner
**Time**: ~10 minutes

### 02. Anharmonic DX Center (`notebooks/02_anharmonic_dx_center.ipynb`)
Real-world example using the DX center in GaAs. Covers:
- Loading experimental/DFT data
- Fitting anharmonic potentials (spline)
- Configuration coordinate diagrams
- Comparison with literature values

**Difficulty**: Intermediate
**Time**: ~20 minutes

### 03. Parameter Scanning (`notebooks/03_parameter_scan.ipynb`)
High-throughput materials screening. Covers:
- Setting up parameter sweeps (ΔQ, ΔE, ℏω)
- Parallel execution strategies
- Visualizing 2D heatmaps
- Identifying optimal materials

**Difficulty**: Intermediate
**Time**: ~15 minutes

### 04. Interactive Dashboard (`notebooks/04_interactive_viz.ipynb`)
Using the web-based visualization dashboard. Covers:
- Launching the Dash server
- Interactive fitting and parameter exploration
- Real-time capture coefficient calculation
- Exporting publication-quality figures

**Difficulty**: Beginner
**Time**: ~15 minutes

## Example Data

The `data/` directory contains sample input files:
- `dx_center_*.dat` - DX center in GaAs potential energy data
- `sn_zn_*.dat` - Sn_Zn defect in ZnO
- Example configuration files for CLI usage

## Running the Notebooks

### Installation
```bash
# Install CarrierCapture with notebook dependencies
pip install -e ".[dev]"

# Or install jupyter separately
pip install jupyter matplotlib
```

### Launch Jupyter
```bash
cd examples/notebooks
jupyter notebook
```

### Run notebooks in order
The notebooks are designed to be followed sequentially, but each is self-contained.

## Command-Line Examples

### Quick Start: Harmonic Oscillator
```bash
# Generate harmonic potential and calculate capture coefficient
python -c "
from carriercapture.core.potential import Potential
from carriercapture.core.config_coord import ConfigCoordinate
import numpy as np

# Create two harmonic potentials
pot_i = Potential.from_harmonic(hw=0.008, Q0=0.0, E0=0.5)
pot_f = Potential.from_harmonic(hw=0.008, Q0=10.0, E0=0.0)

# Solve Schrödinger equation
pot_i.solve(nev=180)
pot_f.solve(nev=60)

# Calculate capture coefficient
cc = ConfigCoordinate(pot_i=pot_i, pot_f=pot_f, W=0.068)
cc.calculate_overlap(Q0=5.0, cutoff=0.25, sigma=0.01)
cc.calculate_capture_coefficient(volume=1e-21, temperature=np.linspace(100, 500, 50))

print(f'Capture coefficient at 300K: {cc.capture_coefficient[20]:.3e} cm³/s')
"
```

### Full CLI Workflow
```bash
# 1. Fit potential from data
carriercapture fit data/dx_center_excited.dat -f spline -o 4 -s 0.001 -O excited.json

# 2. Solve Schrödinger equation
carriercapture solve excited.json -n 180 -O excited_solved.json

# 3. Calculate capture coefficient
carriercapture capture config.yaml -V 1e-21 --temp-range 100 500 50 -O results.json

# 4. Visualize results
carriercapture plot results.json --show

# 5. Parameter scan
carriercapture scan --dQ-min 0 --dQ-max 25 --dQ-points 25 \
                    --dE-min 0 --dE-max 2.5 --dE-points 10 \
                    -j 4 -o scan_results.npz

# 6. Plot scan results
carriercapture scan-plot scan_results.npz --log-scale --show

# 7. Launch interactive dashboard
carriercapture viz --port 8050
```

## Validation Against Julia

The `validation/` directory contains test cases that compare Python results against the original Julia implementation:
- Harmonic oscillator eigenvalues
- DX center capture coefficients
- Sn_Zn example from literature

All validation tests should pass within specified tolerances:
- Eigenvalues: `rtol=1e-4`
- Capture coefficients: `rtol=1e-2`

## Additional Resources

- **Documentation**: https://carriercapture.readthedocs.io
- **Julia version**: https://github.com/WMD-group/CarrierCapture.jl
- **Paper**: Alkauskas et al., Phys. Rev. B (2014) - Multiphonon theory

## Contributing

Found a bug or have an example to share? Please open an issue or PR at:
https://github.com/WMD-group/CarrierCapture.py

## Citation

If you use CarrierCapture in your research, please cite:

```bibtex
@article{alkauskas2014,
  title={First-principles calculations of luminescence spectrum line shapes for defects in semiconductors},
  author={Alkauskas, Audrius and Yan, Qimin and Van de Walle, Chris G},
  journal={Physical Review B},
  volume={90},
  number={7},
  pages={075202},
  year={2014},
  publisher={APS}
}
```
