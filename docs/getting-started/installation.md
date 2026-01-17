# Installation

This guide will help you install CarrierCapture.py on your system.

## System Requirements

- Python 3.9 or later (Python 3.11-3.12 recommended)
- pip package manager
- Git (for development installation)

### Platform Support

CarrierCapture.py works on:

- **Linux** (tested on Ubuntu 20.04+)
- **macOS** (tested on macOS 12+)
- **Windows** (via WSL or native Python)

## Quick Installation (PyPI)

!!! info "Coming Soon"
    CarrierCapture.py will be available on PyPI. For now, install from source.

```bash
# Future release
pip install carriercapture
```

## Installation from Source

### 1. Clone the Repository

```bash
git clone https://github.com/WMD-group/CarrierCapture.py.git
cd CarrierCapture.py
```

### 2. Create Virtual Environment (Recommended)

Using `venv`:

```bash
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
```

Using `conda`:

```bash
conda create -n carriercapture python=3.11
conda activate carriercapture
```

### 3. Install Package

For users:

```bash
pip install -e .
```

For developers:

```bash
pip install -e ".[dev]"
```

## Optional Dependencies

### Interactive Visualization

For the Dash dashboard and advanced plotting:

```bash
pip install carriercapture[viz]
```

This includes:

- Plotly
- Dash
- Additional visualization tools

### doped Integration

For defect calculation workflows with [doped](https://github.com/SMTG-Bham/doped):

```bash
pip install carriercapture[doped]
```

This includes:

- doped package
- monty for serialization
- pymatgen (already included in base install)

### Documentation Building

For building documentation locally:

```bash
pip install carriercapture[docs]
```

### Jupyter Notebooks

For running tutorial notebooks:

```bash
pip install carriercapture[notebook]
```

### All Optional Dependencies

Install everything:

```bash
pip install -e ".[dev,docs,notebook,doped,viz]"
```

## Verify Installation

After installation, verify that CarrierCapture works:

### Check CLI

```bash
carriercapture --help
```

You should see the help message with available commands.

### Check Python API

```bash
python -c "import carriercapture; print(carriercapture.__version__)"
```

### Run Tests

```bash
pytest tests/ -v
```

## Troubleshooting

### Import Errors

If you get import errors, ensure your virtual environment is activated:

```bash
which python  # Should point to your venv
pip list | grep carriercapture  # Should show installed version
```

### SciPy/NumPy Installation Issues

On some systems, you may need to install system dependencies first:

**Ubuntu/Debian:**

```bash
sudo apt-get install python3-dev gfortran libopenblas-dev liblapack-dev
```

**macOS:**

```bash
brew install gcc openblas lapack
```

### Permission Denied

If you get permission errors, use `pip install --user` or create a virtual environment.

### Windows-Specific Issues

On Windows, we recommend using Windows Subsystem for Linux (WSL) for the best experience. Alternatively, ensure you have:

- Visual C++ Build Tools installed
- Latest version of pip: `python -m pip install --upgrade pip`

## Updating CarrierCapture

### From Source

```bash
cd CarrierCapture.py
git pull
pip install -e . --upgrade
```

### From PyPI (Future)

```bash
pip install --upgrade carriercapture
```

## Uninstallation

To remove CarrierCapture:

```bash
pip uninstall carriercapture
```

## Development Installation

For contributing to CarrierCapture, install in editable mode with all development tools:

```bash
git clone https://github.com/WMD-group/CarrierCapture.py.git
cd CarrierCapture.py
pip install -e ".[dev]"
pre-commit install  # Optional: set up git hooks
```

This installs:

- pytest for testing
- black, ruff for code formatting
- mypy for type checking
- All other development dependencies

## Next Steps

Now that you have CarrierCapture installed, proceed to:

- **[Quick Start](quick-start.md)** - Run your first calculation in 5 minutes
- **[Basic Concepts](basic-concepts.md)** - Understand key concepts
- **[User Guide](../user-guide/potentials.md)** - Learn how to use CarrierCapture

## Getting Help

If you encounter issues during installation:

1. Check the [GitHub Issues](https://github.com/WMD-group/CarrierCapture.py/issues)
2. Search for your error message
3. Open a new issue with:
   - Your operating system
   - Python version (`python --version`)
   - Error message
   - Steps to reproduce
