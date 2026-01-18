#!/usr/bin/env python3
"""
CarrierCapture.py Benchmark vs Julia Reference
===============================================

Runs the Sn_Zn in ZnO example using CarrierCapture.py and compares
results against Julia reference data to validate numerical accuracy.

Test Case: Sn substituting Zn in ZnO
Parameters: From examples/notebooks/01_harmonic_sn_zn.ipynb
"""

import json
import sys
from pathlib import Path

import numpy as np

# Add src to path for local development
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from carriercapture.core.potential import Potential
from carriercapture.core.config_coord import ConfigCoordinate


def compare_arrays(python_vals, julia_vals, name, rtol):
    """
    Compare arrays and return detailed comparison.

    Parameters
    ----------
    python_vals : array-like
        Python results
    julia_vals : array-like
        Julia reference results
    name : str
        Name of the comparison
    rtol : float
        Relative tolerance

    Returns
    -------
    dict
        Comparison results including pass/fail status
    """
    python_vals = np.array(python_vals)
    julia_vals = np.array(julia_vals)

    abs_diff = np.abs(python_vals - julia_vals)
    rel_diff = abs_diff / np.abs(julia_vals)

    max_rel_diff = np.max(rel_diff)
    passed = bool(max_rel_diff < rtol)

    return {
        "name": name,
        "passed": passed,
        "max_relative_difference": float(max_rel_diff),
        "tolerance": rtol,
        "max_absolute_difference": float(np.max(abs_diff)),
        "mean_relative_difference": float(np.mean(rel_diff))
    }


def main():
    print("=" * 60)
    print("CarrierCapture.jl vs CarrierCapture.py Benchmark")
    print("=" * 60)
    print("\nTest Case: Sn_Zn in ZnO (Harmonic Approximation)")

    # Load Julia reference data
    ref_path = Path(__file__).parent / "reference_data" / "sn_zn_julia_reference.json"

    if not ref_path.exists():
        print(f"\n✗ ERROR: Julia reference data not found at {ref_path}")
        print("\nPlease run Julia reference first:")
        print("  julia benchmarks/run_julia_reference.jl")
        sys.exit(1)

    print(f"\nLoading Julia reference data from:")
    print(f"  {ref_path}")

    with open(ref_path) as f:
        julia_results = json.load(f)

    # Extract parameters
    params = julia_results["parameters"]
    print("\nParameters:")
    print(f"  ℏω = {params['hw']} eV")
    print(f"  ΔQ = {params['dQ']} amu^0.5·Å")
    print(f"  ΔE = {params['dE']} eV")
    print(f"  W = {params['W']} eV")
    print(f"  Volume = {params['volume']} cm³")
    print(f"  Temperature = {params['temperature']} K")
    print(f"  Grid points = {params['npoints']}")

    # Run Python calculations
    print("\nStep 1: Creating harmonic potentials...")

    # Initial state (excited): Q0=0.0, E0=0.5 eV
    pot_initial = Potential.from_harmonic(
        hw=params['hw'],
        Q0=0.0,
        E0=params['dE'],
        Q_range=(params['Q_range'][0], params['Q_range'][1]),
        npoints=params['npoints']
    )

    # Final state (ground): Q0=10.5, E0=0.0 eV
    pot_final = Potential.from_harmonic(
        hw=params['hw'],
        Q0=params['dQ'],
        E0=0.0,
        Q_range=(params['Q_range'][0], params['Q_range'][1]),
        npoints=params['npoints']
    )

    print("  Initial state: E0=0.5 eV, Q0=0.0")
    print("  Final state: E0=0.0 eV, Q0=10.5")

    # Solve Schrödinger equation
    print("\nStep 2: Solving Schrödinger equation...")
    pot_initial.solve(nev=params['nev_initial'])
    pot_final.solve(nev=params['nev_final'])

    print(f"  Initial state: Found {len(pot_initial.eigenvalues)} eigenvalues")
    print(f"    E₀ = {pot_initial.eigenvalues[0]:.6f} eV")
    print(f"    E₁ = {pot_initial.eigenvalues[1]:.6f} eV")
    print(f"    E₂ = {pot_initial.eigenvalues[2]:.6f} eV")

    print(f"  Final state: Found {len(pot_final.eigenvalues)} eigenvalues")
    print(f"    E₀ = {pot_final.eigenvalues[0]:.6f} eV")
    print(f"    E₁ = {pot_final.eigenvalues[1]:.6f} eV")
    print(f"    E₂ = {pot_final.eigenvalues[2]:.6f} eV")

    # Calculate capture coefficient
    print("\nStep 3: Calculating capture coefficient...")

    # Use same crossing point as Julia if available
    Q0_crossing = params.get('Q0_crossing', 5.0)

    cc = ConfigCoordinate(
        pot_i=pot_initial,
        pot_f=pot_final,
        W=params['W']
    )

    cc.calculate_overlap(Q0=Q0_crossing, sigma=0.025)
    cc.calculate_capture_coefficient(
        volume=params['volume'],
        temperature=np.array([params['temperature']])
    )

    C_300K_python = cc.capture_coefficient[0]

    print(f"  Overlap matrix: {cc.overlap_matrix.shape}")
    print(f"  Q0 (crossing) = {Q0_crossing} amu^0.5·Å")
    print(f"  C(300K) = {C_300K_python:.6e} cm³/s")

    # Compare results
    print("\n" + "=" * 60)
    print("Comparison Results")
    print("=" * 60)

    # Compare initial eigenvalues
    n_compare = min(20, len(pot_initial.eigenvalues), len(julia_results["eigenvalues_initial"]))
    eig_initial_comp = compare_arrays(
        pot_initial.eigenvalues[:n_compare],
        julia_results["eigenvalues_initial"][:n_compare],
        "Initial eigenvalues",
        rtol=1e-4
    )

    print(f"\n1. Initial Eigenvalues (first {n_compare} states):")
    print(f"   Max relative diff:  {eig_initial_comp['max_relative_difference']:.2e}")
    print(f"   Mean relative diff: {eig_initial_comp['mean_relative_difference']:.2e}")
    print(f"   Max absolute diff:  {eig_initial_comp['max_absolute_difference']:.2e} eV")
    print(f"   Tolerance:          {eig_initial_comp['tolerance']:.2e}")
    print(f"   Status:             {'✓ PASS' if eig_initial_comp['passed'] else '✗ FAIL'}")

    # Compare final eigenvalues
    n_compare_f = min(20, len(pot_final.eigenvalues), len(julia_results["eigenvalues_final"]))
    eig_final_comp = compare_arrays(
        pot_final.eigenvalues[:n_compare_f],
        julia_results["eigenvalues_final"][:n_compare_f],
        "Final eigenvalues",
        rtol=1e-4
    )

    print(f"\n2. Final Eigenvalues (first {n_compare_f} states):")
    print(f"   Max relative diff:  {eig_final_comp['max_relative_difference']:.2e}")
    print(f"   Mean relative diff: {eig_final_comp['mean_relative_difference']:.2e}")
    print(f"   Max absolute diff:  {eig_final_comp['max_absolute_difference']:.2e} eV")
    print(f"   Tolerance:          {eig_final_comp['tolerance']:.2e}")
    print(f"   Status:             {'✓ PASS' if eig_final_comp['passed'] else '✗ FAIL'}")

    # Compare capture coefficient
    C_300K_julia = julia_results["capture_coefficient_300K"]
    capture_rel_diff = abs(C_300K_python - C_300K_julia) / abs(C_300K_julia)
    capture_passed = bool(capture_rel_diff < 1e-2)

    capture_comp = {
        "name": "Capture coefficient at 300K",
        "python_value": float(C_300K_python),
        "julia_value": C_300K_julia,
        "relative_difference": float(capture_rel_diff),
        "absolute_difference": float(abs(C_300K_python - C_300K_julia)),
        "tolerance": 1e-2,
        "passed": capture_passed
    }

    print(f"\n3. Capture Coefficient (300K):")
    print(f"   Python:      {capture_comp['python_value']:.6e} cm³/s")
    print(f"   Julia:       {capture_comp['julia_value']:.6e} cm³/s")
    print(f"   Relative diff: {capture_comp['relative_difference']:.2e}")
    print(f"   Absolute diff: {capture_comp['absolute_difference']:.2e} cm³/s")
    print(f"   Tolerance:     {capture_comp['tolerance']:.2e}")
    print(f"   Status:        {'✓ PASS' if capture_comp['passed'] else '✗ FAIL'}")

    # Overall status
    overall_passed = all([
        eig_initial_comp["passed"],
        eig_final_comp["passed"],
        capture_comp["passed"]
    ])

    # Save detailed report
    report = {
        "test_case": "Sn_Zn in ZnO (Harmonic)",
        "parameters": params,
        "comparisons": {
            "eigenvalues_initial": eig_initial_comp,
            "eigenvalues_final": eig_final_comp,
            "capture_coefficient": capture_comp
        },
        "overall_passed": overall_passed
    }

    report_path = Path(__file__).parent / "results" / "benchmark_report.json"
    report_path.parent.mkdir(parents=True, exist_ok=True)

    with open(report_path, "w") as f:
        json.dump(report, f, indent=2)

    print(f"\nDetailed report saved to:")
    print(f"  {report_path}")

    # Final summary
    print("\n" + "=" * 60)
    if overall_passed:
        print("Overall: ✓ ALL TESTS PASSED")
        print("\nConclusion: Python implementation matches Julia results")
        print("within numerical precision!")
    else:
        print("Overall: ✗ SOME TESTS FAILED")
        print("\nSome comparisons exceeded tolerance thresholds.")
        print("Check the detailed report for more information.")
    print("=" * 60 + "\n")

    # Exit with appropriate code
    sys.exit(0 if overall_passed else 1)


if __name__ == "__main__":
    main()
