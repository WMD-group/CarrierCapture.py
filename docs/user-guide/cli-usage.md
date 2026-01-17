# CLI Usage Guide

Practical guide to using CarrierCapture from the command line.

## Overview

The CarrierCapture CLI provides a complete workflow from raw data to capture coefficients:

```bash
carriercapture <command> [options]
```

**Available commands:**
- `fit` - Fit potential energy surfaces
- `solve` - Solve Schrödinger equation
- `capture` - Calculate capture coefficients
- `scan` - High-throughput parameter screening
- `scan-plot` - Visualize scan results
- `viz` - Launch interactive dashboard
- `plot` - Generate static plots

---

## Getting Help

```bash
# General help
carriercapture --help

# Command-specific help
carriercapture fit --help
carriercapture solve --help
carriercapture capture --help

# Show version
carriercapture --version
```

---

## Common Workflows

### Workflow 1: Complete Analysis (Fit → Solve → Capture)

Starting from Q-E data files to capture coefficient:

```bash
# Step 1: Fit excited state potential
carriercapture fit excited_data.dat \
  -f spline --order 4 --smoothness 0.001 \
  -o excited.json -v

# Step 2: Fit ground state potential
carriercapture fit ground_data.dat \
  -f spline --order 4 --smoothness 0.001 \
  -o ground.json -v

# Step 3: Solve for excited state (more eigenvalues)
carriercapture solve excited.json \
  -n 180 -o excited_solved.json -v

# Step 4: Solve for ground state
carriercapture solve ground.json \
  -n 60 -o ground_solved.json -v

# Step 5: Calculate capture coefficient
carriercapture capture \
  --pot-i excited_solved.json \
  --pot-f ground_solved.json \
  -W 0.205 -V 1e-21 --Q0 10.0 \
  --temp-range 100 500 50 \
  -o capture_results.json \
  --plot --plot-output arrhenius.png \
  -vv
```

**Expected output:**
```
Loading initial potential: excited_solved.json
Loading final potential: ground_solved.json

Initial potential: 180 states
Final potential: 60 states

Capture parameters:
  W (coupling): 0.205 eV
  g (degeneracy): 1
  V (volume): 1.00e-21 cm³
  Q0: 10.0 amu^0.5·Å
  Energy cutoff: 0.25 eV
  Delta width: 0.025 eV

Calculating wavefunction overlaps...
✓ Computed 2847/10800 non-zero overlaps (26.4%)

Calculating capture coefficient...
  Temperature: 100 - 500 K (50 points)
✓ Calculation completed successfully

Results:
  C(T=100K) = 1.234e-12 cm³/s
  C(T=300K) = 5.678e-11 cm³/s
  C(T=500K) = 2.345e-10 cm³/s

Saving results to: capture_results.json
✓ Saved successfully

✓ Plot saved to: arrhenius.png
```

### Workflow 2: Using Config File

More convenient for repeated calculations:

**Create `config.yaml`:**
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

**Run:**
```bash
carriercapture capture config.yaml \
  -o results.json --plot -v
```

### Workflow 3: Parameter Screening

Screen materials across (ΔQ, ΔE) space:

```bash
# Run parameter scan
carriercapture scan \
  --dQ-min 0 --dQ-max 25 --dQ-points 25 \
  --dE-min 0 --dE-max 2.5 --dE-points 10 \
  --hbar-omega-i 0.008 --hbar-omega-f 0.008 \
  -T 300 -V 1e-21 -g 1 \
  -j -1 \
  -o scan_results.npz -v

# Visualize results
carriercapture scan-plot scan_results.npz \
  --log-scale --show

# Save visualization
carriercapture scan-plot scan_results.npz \
  --log-scale -o heatmap.html
```

**Output:**
```
Setting up parameter scan...
  ΔQ: 25 points from 0 to 25 amu^0.5·Å
  ΔE: 10 points from 0 to 2.5 eV
  ℏω_i = 0.0080 eV, ℏω_f = 0.0080 eV
  Temperature: 300.0 K
  Volume: 1.00e-21 cm³
  Parallel jobs: 8

Starting scan...
100%|████████████| 250/250 [02:15<00:00,  1.85it/s]

Saving results to: scan_results.npz
✓ Scan complete!

Results summary:
  Grid size: 25 × 10 = 250 points
  Successful calculations: 250/250
  Capture coefficient range:
    Min: 1.234e-16 cm³/s
    Max: 5.678e-10 cm³/s
    Mean: 3.456e-12 cm³/s
```

### Workflow 4: Quick Visualization

Visualize fitted and solved potentials:

```bash
# Plot potential with wavefunctions
carriercapture plot potential_solved.json \
  --show-wf --max-wf 10 --show

# Plot eigenvalue spectrum
carriercapture plot potential_solved.json \
  --type spectrum --show

# Generate both plots and save
carriercapture plot potential_solved.json \
  --type both --show-wf \
  -o potential_figure.html
```

### Workflow 5: Interactive Exploration

Launch dashboard for interactive analysis:

```bash
# Basic launch
carriercapture viz

# Custom port
carriercapture viz --port 8080

# With data preloaded
carriercapture viz --data potential.json

# Debug mode (for development)
carriercapture viz --debug --no-browser
```

---

## Practical Tips

### 1. Verbosity Levels

Control output detail:

```bash
# Quiet (minimal output)
carriercapture fit data.dat -o fit.json

# Standard (-v): Show progress
carriercapture fit data.dat -o fit.json -v

# Verbose (-vv): Detailed information
carriercapture fit data.dat -o fit.json -vv

# Debug (-vvv): Maximum detail
carriercapture fit data.dat -o fit.json -vvv
```

### 2. Chaining Commands with Shell

```bash
# Chain commands with &&
carriercapture fit excited.dat -f spline -o excited.json -v && \
carriercapture solve excited.json -n 180 -o excited_solved.json -v && \
echo "✓ Excited state ready"

# Process multiple files
for file in data/*.dat; do
    basename=$(basename "$file" .dat)
    carriercapture fit "$file" -f spline -o "fitted/${basename}.json" -v
done
```

### 3. Batch Processing Script

**`process_all.sh`:**
```bash
#!/bin/bash
# Process all defect data

set -e  # Exit on error

DEFECTS=("Sn_Zn" "Sn_O" "V_Zn")

for defect in "${DEFECTS[@]}"; do
    echo "Processing $defect..."

    # Fit
    carriercapture fit data/${defect}_excited.dat \
        -f spline -o ${defect}_excited.json -v
    carriercapture fit data/${defect}_ground.dat \
        -f spline -o ${defect}_ground.json -v

    # Solve
    carriercapture solve ${defect}_excited.json -n 180 \
        -o ${defect}_excited_solved.json -v
    carriercapture solve ${defect}_ground.json -n 60 \
        -o ${defect}_ground_solved.json -v

    # Capture
    carriercapture capture \
        --pot-i ${defect}_excited_solved.json \
        --pot-f ${defect}_ground_solved.json \
        -W 0.205 -V 1e-21 --Q0 10.0 \
        -o results/${defect}_capture.json \
        --plot --plot-output results/${defect}_arrhenius.png \
        -v

    echo "✓ $defect complete"
    echo
done

echo "✓ All defects processed"
```

Run:
```bash
chmod +x process_all.sh
./process_all.sh
```

### 4. Output Redirection

```bash
# Save output to log file
carriercapture capture config.yaml -vv > capture.log 2>&1

# Separate stdout and stderr
carriercapture capture config.yaml -v > output.txt 2> errors.txt

# Append to log
carriercapture fit data.dat -o fit.json -v >> workflow.log 2>&1
```

### 5. Parallel Scans

```bash
# Use all CPU cores
carriercapture scan \
  --dQ-min 0 --dQ-max 25 --dQ-points 50 \
  --dE-min 0 --dE-max 2.5 --dE-points 20 \
  -j -1 \
  -o scan_high_res.npz

# Specify core count
carriercapture scan ... -j 8 -o scan.npz

# Monitor with top/htop in another terminal
# to see CPU usage
```

---

## File Format Guide

### Input Files

**Q-E Data (CSV/DAT):**
```
# Configuration coordinate (amu^0.5·Å), Energy (eV)
0.0, 0.0
2.5, 0.05
5.0, 0.18
7.5, 0.35
10.0, 0.48
...
```

Can be space or comma separated:
```
0.0 0.0
2.5 0.05
5.0 0.18
```

**Config File (YAML):**
```yaml
potential_initial:
  file: excited_solved.json

potential_final:
  file: ground_solved.json

capture:
  W: 0.205
  degeneracy: 1
  volume: 1.0e-21
  Q0: 10.0
  cutoff: 0.25
  sigma: 0.025
  temperature:
    min: 100
    max: 500
    n_points: 50
```

### Output Files

**JSON (Human-readable):**
```bash
carriercapture fit data.dat -o potential.json
```

Contains: Q_data, E_data, fit parameters, eigenvalues, eigenvectors (if solved)

**NPZ (NumPy, compact):**
```bash
carriercapture solve potential.json -o solved.npz
```

Faster to load, smaller file size for large arrays

**HDF5 (Large datasets):**
```bash
carriercapture scan ... -o scan.h5
```

Efficient for very large scans

---

## Troubleshooting

### Command not found

```bash
$ carriercapture --version
bash: carriercapture: command not found
```

**Fix:**
```bash
# Check if installed
pip show carriercapture

# Reinstall if needed
pip install --upgrade carriercapture

# Check PATH
which python
which pip

# If using conda
conda list | grep carriercapture
```

### Permission denied

```bash
$ carriercapture fit data.dat -o fit.json
Error: Permission denied: fit.json
```

**Fix:**
```bash
# Check write permissions
ls -l $(dirname fit.json)

# Use different output location
carriercapture fit data.dat -o ~/output/fit.json

# Or fix permissions
chmod u+w .
```

### File not found

```bash
$ carriercapture fit my_data.dat -o fit.json
Error: [Errno 2] No such file or directory: 'my_data.dat'
```

**Fix:**
```bash
# Check file exists
ls -l my_data.dat

# Use absolute path
carriercapture fit /full/path/to/my_data.dat -o fit.json

# Check current directory
pwd
ls
```

### Port already in use (viz)

```bash
$ carriercapture viz
Error: Address already in use
```

**Fix:**
```bash
# Use different port
carriercapture viz --port 8051

# Or kill process using port 8050
lsof -ti:8050 | xargs kill -9  # macOS/Linux
```

### Out of memory (large scans)

```bash
$ carriercapture scan ... --dQ-points 100 --dE-points 100 -j -1
MemoryError: Unable to allocate array
```

**Fix:**
```bash
# Reduce grid resolution
carriercapture scan ... --dQ-points 50 --dE-points 50 -j -1

# Use fewer parallel jobs
carriercapture scan ... -j 4

# Or split into multiple scans
carriercapture scan --dQ-min 0 --dQ-max 12.5 ... -o scan_1.npz
carriercapture scan --dQ-min 12.5 --dQ-max 25 ... -o scan_2.npz
```

---

## Environment Variables

CarrierCapture respects these environment variables:

```bash
# Set number of threads for BLAS/LAPACK (used by NumPy/SciPy)
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4

# Default data directory
export CARRIERCAPTURE_DATA_DIR=$HOME/carriercapture_data

# Run command
carriercapture fit data.dat -o $CARRIERCAPTURE_DATA_DIR/fit.json
```

Add to `~/.bashrc` or `~/.zshrc`:
```bash
# CarrierCapture settings
export OMP_NUM_THREADS=8
export CARRIERCAPTURE_DATA_DIR=$HOME/carriercapture_data
```

---

## Advanced Examples

### Example 1: doped Integration

```bash
# Calculate capture from doped VASP calculations
carriercapture capture \
  --doped defects/Sn_Zn_0.json.gz \
  --charge-i 0 --charge-f +1 \
  --doped-path-i vasp/path_q0/ \
  --doped-path-f vasp/path_q1/ \
  -W 0.205 -V 1e-21 \
  --temp-range 100 500 50 \
  --auto-Q0 \
  -o Sn_Zn_capture.json \
  --plot --plot-output Sn_Zn_arrhenius.png \
  -vv
```

### Example 2: Multi-Temperature Scan

```bash
# Scan at different temperatures
for T in 200 300 400 500; do
    echo "Scanning at T=${T}K..."
    carriercapture scan \
        --dQ-min 0 --dQ-max 25 --dQ-points 25 \
        --dE-min 0 --dE-max 2.5 --dE-points 10 \
        -T $T -V 1e-21 -g 1 \
        -j -1 -o scan_T${T}K.npz -v

    # Visualize
    carriercapture scan-plot scan_T${T}K.npz \
        --log-scale -o heatmap_T${T}K.html
done

echo "✓ All temperatures scanned"
```

### Example 3: Convergence Testing

```bash
# Test nev convergence
for nev in 40 60 80 100; do
    echo "Testing nev=${nev}..."
    carriercapture solve potential.json \
        -n $nev -o solved_nev${nev}.json

    carriercapture capture \
        --pot-i solved_nev${nev}.json \
        --pot-f ground_solved.json \
        -W 0.205 -V 1e-21 --Q0 10.0 \
        --temp-range 300 300 1 \
        -o capture_nev${nev}.json

    # Extract C(300K) from JSON
    python -c "import json; \
        data = json.load(open('capture_nev${nev}.json')); \
        print(f'nev={$nev}: C(300K)={data[\"capture_coefficient\"][0]:.3e} cm³/s')"
done
```

---

## Shell Integration

### Bash Completion

Create `~/.carriercapture-completion.bash`:
```bash
# CarrierCapture bash completion

_carriercapture() {
    local cur prev opts
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"

    # Main commands
    if [ $COMP_CWORD -eq 1 ]; then
        opts="fit solve capture scan scan-plot viz plot --help --version"
        COMPREPLY=( $(compgen -W "${opts}" -- ${cur}) )
        return 0
    fi

    # File completion for remaining arguments
    COMPREPLY=( $(compgen -f -- ${cur}) )
}

complete -F _carriercapture carriercapture
```

Add to `~/.bashrc`:
```bash
source ~/.carriercapture-completion.bash
```

### Aliases

Add convenient aliases to `~/.bashrc`:
```bash
# CarrierCapture aliases
alias ccfit='carriercapture fit'
alias ccsolve='carriercapture solve'
alias cccapture='carriercapture capture'
alias ccscan='carriercapture scan'
alias ccviz='carriercapture viz'
alias ccplot='carriercapture plot'
```

Use:
```bash
ccfit data.dat -f spline -o fit.json -v
ccsolve fit.json -n 180 -v
ccviz
```

---

## Integration with Other Tools

### With Jupyter Notebooks

```bash
# Run CLI command from notebook
!carriercapture fit data.dat -o fit.json -v

# Or use %%bash magic
%%bash
carriercapture capture config.yaml -o results.json -v
```

### With Makefiles

**`Makefile`:**
```makefile
.PHONY: all fit solve capture clean

all: capture

fit: excited.json ground.json

%.json: data/%.dat
	carriercapture fit $< -f spline -o $@ -v

solve: excited_solved.json ground_solved.json

%_solved.json: %.json
	carriercapture solve $< -n 180 -o $@ -v

capture: capture_results.json

capture_results.json: excited_solved.json ground_solved.json config.yaml
	carriercapture capture config.yaml -o $@ --plot -v

clean:
	rm -f *.json *.png *.html
```

Use:
```bash
make all  # Run complete workflow
make clean  # Clean up
```

---

## See Also

- **[API Reference: CLI](../api/cli.md)** - Complete CLI command reference
- **[User Guide: Configuration](configuration.md)** - YAML configuration files
- **[Getting Started: Quick Start](../getting-started/quick-start.md)** - Quick introduction
- **[Examples](../examples/notebooks.md)** - Example calculations
