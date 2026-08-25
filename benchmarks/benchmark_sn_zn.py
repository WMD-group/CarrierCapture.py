#!/usr/bin/env python3
"""
CarrierCapture.py Benchmark vs Julia Reference
===============================================

Runs the Sn_Zn in ZnO example using CarrierCapture.py and compares
results against CarrierCapture.jl reference data in three tiers:

1. Algorithmic equivalence (binding): Python rebuilt with CarrierCapture.jl's
   numerical conventions must reproduce the Julia reference to near machine
   precision. This proves both codes implement the same physics.
2. Native accuracy (binding): Python's native eigenvalues must match the
   analytic harmonic result E_n = E0 + hw*(n + 1/2).
3. Native vs Julia (informational): the native results differ from Julia by
   ~1.5% at this grid, entirely due to two CarrierCapture.jl conventions:
   its finite-difference kinetic term uses grid spacing dq = (Q_max-Q_min)/N
   while its grid range(Qi, Qf, length=N) actually has spacing
   (Q_max-Q_min)/(N-1), and it integrates overlaps with the rectangle rule.
   Python uses the true grid spacing and the trapezoid rule, and is closer
   to the analytic eigenvalues. Both codes converge to the same answer as
   npoints -> infinity.
"""

import json
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

# Add src to path for local development
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from carriercapture.core.potential import Potential
from carriercapture._constants import AMU, HBAR_C, HBAR, K_B


def solve_julia_convention(pot, nev, npoints, q_range):
    """
    Solve with CarrierCapture.jl's solve1D_Quantum grid convention.

    CarrierCapture.jl builds its finite-difference kinetic term with
    dq = (Q_max - Q_min) / N, while its grid range(Qi, Qf, length=N) has
    true spacing (Q_max - Q_min) / (N - 1). Reproducing that off-by-one
    here lets us verify algorithmic equivalence to machine precision.
    """
    Q = np.linspace(*q_range, npoints)
    dq = (Q[-1] - Q[0]) / npoints
    kinetic = (HBAR_C * 1e10) ** 2 / AMU / (2 * dq**2)
    H = sp.diags(
        [np.full(npoints - 1, -kinetic), 2 * kinetic + pot(Q), np.full(npoints - 1, -kinetic)],
        [-1, 0, 1],
    )
    vals, vecs = eigsh(H.tocsc(), k=nev, which="SA", tol=0)
    order = np.argsort(vals)
    vals, vecs = vals[order], vecs[:, order]
    vecs /= np.sqrt(dq * np.sum(vecs**2, axis=0))
    return vals, vecs, Q, dq


def capture_julia_convention(params):
    """C(300K) with Julia's grid convention and rectangle-rule overlaps."""
    pot_i = Potential.from_harmonic(
        hw=params["hw"], Q0=0.0, E0=params["dE"],
        Q_range=tuple(params["Q_range"]), npoints=params["npoints"],
    )
    pot_f = Potential.from_harmonic(
        hw=params["hw"], Q0=params["dQ"], E0=0.0,
        Q_range=tuple(params["Q_range"]), npoints=params["npoints"],
    )
    vals_i, vecs_i, Q, dq = solve_julia_convention(
        pot_i, params["nev_initial"], params["npoints"], params["Q_range"]
    )
    vals_f, vecs_f, _, _ = solve_julia_convention(
        pot_f, params["nev_final"], params["npoints"], params["Q_range"]
    )

    # Rectangle-rule overlaps S_ij = dq * sum(psi_i * (Q - Q0) * psi_j)
    # (errstate: subnormal wavefunction tails trip spurious BLAS warnings)
    operator = Q - params["Q0_crossing"]
    with np.errstate(all="ignore"):
        S = dq * (vecs_i * operator[:, None]).T @ vecs_f
    dE = vals_i[:, None] - vals_f[None, :]
    delta = np.exp(-(dE**2) / (2 * params["sigma"] ** 2)) / (params["sigma"] * np.sqrt(2 * np.pi))
    mask = np.abs(dE) < params["cutoff"]
    S, delta = np.where(mask, S, 0.0), np.where(mask, delta, 0.0)

    beta = 1.0 / (K_B * params["temperature"])
    occupation = np.exp(-beta * vals_i)
    occupation /= occupation.sum()
    prefactor = params["volume"] * 2 * np.pi / HBAR * params["W"] ** 2
    C = prefactor * np.sum(occupation[:, None] * S**2 * delta)
    return vals_i, vals_f, C


def capture_native(params):
    """C(300K) with CarrierCapture.py's native (more accurate) numerics."""
    pot_i = Potential.from_harmonic(
        hw=params["hw"], Q0=0.0, E0=params["dE"],
        Q_range=tuple(params["Q_range"]), npoints=params["npoints"],
    )
    pot_f = Potential.from_harmonic(
        hw=params["hw"], Q0=params["dQ"], E0=0.0,
        Q_range=tuple(params["Q_range"]), npoints=params["npoints"],
    )
    pot_i.solve(nev=params["nev_initial"])
    pot_f.solve(nev=params["nev_final"])

    from carriercapture.core.config_coord import ConfigCoordinate

    cc = ConfigCoordinate(pot_i=pot_i, pot_f=pot_f, W=params["W"])
    cc.calculate_overlap(Q0=params["Q0_crossing"], cutoff=params["cutoff"], sigma=params["sigma"])
    cc.calculate_capture_coefficient(
        volume=params["volume"], temperature=np.array([params["temperature"]])
    )
    return pot_i.eigenvalues, pot_f.eigenvalues, cc.capture_coefficient[0]


def compare(python_vals, reference_vals, name, rtol, n=None):
    """Relative comparison of arrays or scalars."""
    python_vals = np.atleast_1d(np.asarray(python_vals, dtype=float))
    reference_vals = np.atleast_1d(np.asarray(reference_vals, dtype=float))
    if n is not None:
        python_vals, reference_vals = python_vals[:n], reference_vals[:n]
    rel = np.abs(python_vals - reference_vals) / np.abs(reference_vals)
    return {
        "name": name,
        "passed": bool(np.max(rel) < rtol),
        "max_relative_difference": float(np.max(rel)),
        "tolerance": rtol,
    }


def main():
    print("=" * 60)
    print("CarrierCapture.jl vs CarrierCapture.py Benchmark")
    print("=" * 60)
    print("\nTest Case: Sn_Zn in ZnO (Harmonic Approximation)")

    ref_path = Path(__file__).parent / "reference_data" / "sn_zn_julia_reference.json"
    if not ref_path.exists():
        print(f"\nERROR: Julia reference data not found at {ref_path}")
        print("Run: julia benchmarks/run_julia_reference.jl")
        sys.exit(1)

    with open(ref_path) as f:
        julia = json.load(f)
    params = julia["parameters"]
    n_states = 20

    print(f"\nParameters: hw={params['hw']} eV, dQ={params['dQ']} amu^0.5*Ang, "
          f"dE={params['dE']} eV, W={params['W']} eV/(amu^0.5*Ang), "
          f"npoints={params['npoints']}")

    # --- Tier 1: algorithmic equivalence (Julia conventions emulated) ---
    print("\nTier 1: Algorithmic equivalence (Julia conventions emulated)")
    em_i, em_f, C_emulated = capture_julia_convention(params)
    tier1 = [
        compare(em_i, julia["eigenvalues_initial"], "Initial eigenvalues (emulated)", 1e-9, n_states),
        compare(em_f, julia["eigenvalues_final"], "Final eigenvalues (emulated)", 1e-9, n_states),
        compare(C_emulated, julia["capture_coefficient_300K"], "C(300K) (emulated)", 1e-6),
    ]
    for c in tier1:
        print(f"  {c['name']}: max rel diff {c['max_relative_difference']:.2e} "
              f"(tol {c['tolerance']:.0e}) {'PASS' if c['passed'] else 'FAIL'}")

    # --- Tier 2: native accuracy vs analytic harmonic eigenvalues ---
    print("\nTier 2: Native accuracy vs analytic E_n = E0 + hw*(n + 1/2)")
    nat_i, nat_f, C_native = capture_native(params)
    n_arr = np.arange(n_states)
    analytic_i = params["dE"] + params["hw"] * (n_arr + 0.5)
    analytic_f = params["hw"] * (n_arr + 0.5)
    tier2 = [
        compare(nat_i, analytic_i, "Initial eigenvalues vs analytic", 5e-4, n_states),
        compare(nat_f, analytic_f, "Final eigenvalues vs analytic", 5e-4, n_states),
    ]
    for c in tier2:
        print(f"  {c['name']}: max rel diff {c['max_relative_difference']:.2e} "
              f"(tol {c['tolerance']:.0e}) {'PASS' if c['passed'] else 'FAIL'}")

    # --- Tier 3: native vs Julia (informational) ---
    print("\nTier 3: Native C(300K) vs Julia (informational)")
    C_julia = julia["capture_coefficient_300K"]
    tier3 = compare(C_native, C_julia, "C(300K) native vs Julia", 2e-2)
    tier3["python_value"] = float(C_native)
    tier3["julia_value"] = C_julia
    tier3["explanation"] = (
        "The gap is entirely due to CarrierCapture.jl's numerical conventions "
        "(finite-difference spacing dq = dQ/N vs the true grid spacing dQ/(N-1), "
        "and rectangle-rule overlap integration), verified by Tier 1. Python's "
        "native numerics are closer to the analytic eigenvalues (Tier 2). Both "
        "codes converge to the same answer with increasing npoints."
    )
    print(f"  Python: {C_native:.6e} cm^3/s | Julia: {C_julia:.6e} cm^3/s | "
          f"rel diff {tier3['max_relative_difference']:.2e} "
          f"(tol {tier3['tolerance']:.0e}) {'PASS' if tier3['passed'] else 'FAIL'}")

    overall = all(c["passed"] for c in tier1 + tier2 + [tier3])

    report = {
        "test_case": "Sn_Zn in ZnO (Harmonic)",
        "parameters": params,
        "tier1_algorithmic_equivalence": tier1,
        "tier2_native_vs_analytic": tier2,
        "tier3_native_vs_julia": tier3,
        "overall_passed": overall,
    }
    report_path = Path(__file__).parent / "results" / "benchmark_report.json"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with open(report_path, "w") as f:
        json.dump(report, f, indent=2)
    print(f"\nReport saved to {report_path}")

    print("\n" + "=" * 60)
    print("Overall: " + ("ALL TESTS PASSED" if overall else "SOME TESTS FAILED"))
    print("=" * 60 + "\n")
    sys.exit(0 if overall else 1)


if __name__ == "__main__":
    main()
