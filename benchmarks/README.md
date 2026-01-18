# CarrierCapture.jl vs CarrierCapture.py Benchmark

This directory contains benchmarks to validate that CarrierCapture.py produces identical results to the original Julia implementation.

## Quick Start

Run the full benchmark with a single command:

```bash
./benchmarks/run_benchmark.sh
```

This will:
1. Run the Julia reference calculation
2. Run the Python benchmark and compare results
3. Generate a detailed report

## Requirements

- **Julia** with CarrierCapture.jl installed
- **Python 3.11+** with CarrierCapture.py installed

## Files

- `run_benchmark.sh` - Main runner script (runs both Julia and Python)
- `run_julia_reference.jl` - Julia reference calculation
- `benchmark_sn_zn.py` - Python benchmark + comparison logic
- `reference_data/` - Julia output (JSON format)
- `results/` - Benchmark reports (JSON format)

## Test Case

**Example**: Sn_Zn in ZnO (harmonic approximation)

**Parameters**:
- Phonon energy: ℏω = 8 meV
- Configuration coordinate shift: ΔQ = 10.5 amu^0.5·Å
- Energy offset: ΔE = 0.5 eV
- Electron-phonon coupling: W = 0.068 eV

**Compared Quantities**:
- Initial state eigenvalues (first 20)
- Final state eigenvalues (first 20)
- Capture coefficient at 300K

**Tolerances**:
- Eigenvalues: relative tolerance 1e-4
- Capture coefficient: relative tolerance 1e-2

## Manual Usage

Run steps separately:

```bash
# Step 1: Generate Julia reference data
julia benchmarks/run_julia_reference.jl

# Step 2: Run Python benchmark
python benchmarks/benchmark_sn_zn.py
```

## Output

**Julia Reference** (`reference_data/sn_zn_julia_reference.json`):
```json
{
  "parameters": {...},
  "eigenvalues_initial": [0.504, 0.512, 0.520, ...],
  "eigenvalues_final": [0.004, 0.012, 0.020, ...],
  "capture_coefficient_300K": 1.23e-12
}
```

**Benchmark Report** (`results/benchmark_report.json`):
```json
{
  "test_case": "Sn_Zn in ZnO (Harmonic)",
  "parameters": {...},
  "comparisons": {
    "eigenvalues_initial": {
      "passed": true,
      "max_relative_difference": 1.5e-6,
      "tolerance": 1e-4
    },
    ...
  },
  "overall_passed": true
}
```

## Expected Results

If everything works correctly, you should see:

```
============================================================
CarrierCapture.jl vs CarrierCapture.py Benchmark
============================================================

Test Case: Sn_Zn in ZnO (Harmonic Approximation)

...

============================================================
Comparison Results
============================================================

1. Initial Eigenvalues (first 20 states):
   Max relative diff:  1.23e-06
   Tolerance:          1.00e-04
   Status:             ✓ PASS

2. Final Eigenvalues (first 20 states):
   Max relative diff:  2.34e-06
   Tolerance:          1.00e-04
   Status:             ✓ PASS

3. Capture Coefficient (300K):
   Python:      1.2345e-12 cm³/s
   Julia:       1.2346e-12 cm³/s
   Relative diff: 8.10e-05
   Tolerance:     1.00e-02
   Status:        ✓ PASS

============================================================
Overall: ✓ ALL TESTS PASSED
============================================================
```

## Troubleshooting

### Julia script fails

The Julia script may need adjustments for the actual CarrierCapture.jl API. Check:
- Function names (e.g., `potential()` vs `Potential()`)
- Parameter names and order
- Module structure (`using CarrierCapture` vs submodules)

Consult [CarrierCapture.jl documentation](https://github.com/WMD-group/CarrierCapture.jl).

### Python benchmark fails to find reference data

Make sure you run the Julia reference first:
```bash
julia benchmarks/run_julia_reference.jl
```

### Tests fail (exceed tolerance)

Small differences are expected due to:
- Floating-point rounding differences between languages
- Compiler optimizations
- BLAS/LAPACK library versions

If differences are > 0.1%, investigate:
1. Check eigenvalue magnitudes are reasonable (~0.004-0.5 eV)
2. Verify same parameters used in both implementations
3. Check grid size and numerical integration settings

## Adding to README

After successful benchmark, add results to main README.md:

```markdown
## 🔬 Validation Against Julia

Validated against CarrierCapture.jl for the Sn_Zn in ZnO example:

| Observable | Max Relative Diff | Tolerance | Status |
|------------|-------------------|-----------|--------|
| Initial eigenvalues | < 1e-5 | 1e-4 | ✓ PASS |
| Final eigenvalues | < 1e-5 | 1e-4 | ✓ PASS |
| Capture coefficient (300K) | < 1e-3 | 1e-2 | ✓ PASS |

See `benchmarks/` for benchmark code.
```
