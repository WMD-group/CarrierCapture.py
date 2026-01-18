#!/usr/bin/env julia
"""
CarrierCapture.jl Reference Calculation
========================================

Runs the Sn_Zn in ZnO harmonic example using CarrierCapture.jl
and saves results to JSON for comparison with Python implementation.

Test Case: Sn substituting Zn in ZnO
Parameters: From examples/notebooks/01_harmonic_sn_zn.ipynb
"""

# Activate CarrierCapture.jl environment
using Pkg
carriercapture_path = joinpath(@__DIR__, "..", "..", "CarrierCapture.jl")
Pkg.activate(carriercapture_path)
Pkg.instantiate()  # Install dependencies if needed

# Add JSON package if not already installed
if !haskey(Pkg.project().dependencies, "JSON")
    Pkg.add("JSON")
end

using CarrierCapture
using JSON

println("="^60)
println("CarrierCapture.jl Reference Calculation")
println("="^60)
println("\nTest Case: Sn_Zn in ZnO (Harmonic Approximation)")

# Parameters (matching Python example exactly)
const hw = 0.008  # eV - phonon energy (8 meV)
const Q_min = -20.0  # amu^0.5·Å
const Q_max = 20.0   # amu^0.5·Å
const npoints = 5000
const nev_initial = 180
const nev_final = 60
const W = 0.068  # eV - electron-phonon coupling
const volume = 1e-21  # cm³
const temperature = 300.0  # K
const cut_off = 0.25  # eV
const σ = 0.025  # eV - Gaussian broadening

println("\nParameters:")
println("  ℏω = $hw eV")
println("  ΔQ = 10.5 amu^0.5·Å")
println("  ΔE = 0.5 eV")
println("  W = $W eV")
println("  Volume = $volume cm³")
println("  Temperature = $temperature K")
println("  Grid points = $npoints")

# Create Q range
Q = range(Q_min, stop=Q_max, length=npoints)

# Create potentials
println("\nStep 1: Creating harmonic potentials...")

# Initial state (excited): Q0=0.0, E0=0.5 eV
pot_i = Potential()
pot_i.name = "initial"
pot_i.Q0 = 0.0
pot_i.E0 = 0.5
pot_i.nev = nev_initial
pot_i.func_type = "harmonic"
pot_i.params = Dict("hw" => hw)
pot_i.Q = collect(Q)  # Convert range to array
fit_pot!(pot_i)

# Final state (ground): Q0=10.5, E0=0.0 eV
pot_f = Potential()
pot_f.name = "final"
pot_f.Q0 = 10.5
pot_f.E0 = 0.0
pot_f.nev = nev_final
pot_f.func_type = "harmonic"
pot_f.params = Dict("hw" => hw)
pot_f.Q = collect(Q)  # Convert range to array
fit_pot!(pot_f)

println("  Initial state: E0=0.5 eV, Q0=0.0")
println("  Final state: E0=0.0 eV, Q0=10.5")

# Solve Schrödinger equation
println("\nStep 2: Solving Schrödinger equation...")
solve_pot!(pot_i)
solve_pot!(pot_f)

println("  Initial state: Found $(length(pot_i.ϵ)) eigenvalues")
println("    E₀ = $(pot_i.ϵ[1]) eV")
println("    E₁ = $(pot_i.ϵ[2]) eV")
println("    E₂ = $(pot_i.ϵ[3]) eV")

println("  Final state: Found $(length(pot_f.ϵ)) eigenvalues")
println("    E₀ = $(pot_f.ϵ[1]) eV")
println("    E₁ = $(pot_f.ϵ[2]) eV")
println("    E₂ = $(pot_f.ϵ[3]) eV")

# Calculate capture coefficient
println("\nStep 3: Calculating capture coefficient...")

# Create configuration coordinate
cc = conf_coord(pot_i, pot_f)
cc.W = W
cc.g = 1  # degeneracy

# Find crossing point
Q₀ = 5.0  # Approximate midpoint between Q0_initial (0) and Q0_final (10.5)

# Calculate overlap matrix
calc_overlap!(cc; cut_off = cut_off, σ = σ, Q₀ = Q₀)

# Calculate capture coefficient at 300K
temp_array = [temperature]
calc_capt_coeff!(cc, volume, temp_array)

C_300K = cc.capt_coeff[1]

println("  Overlap matrix: $(size(cc.overlap_matrix))")
println("  Q₀ (crossing) = $Q₀ amu^0.5·Å")
println("  C(300K) = $(C_300K) cm³/s")

# Save results to JSON
println("\nStep 4: Saving results to JSON...")

results = Dict(
    "parameters" => Dict(
        "hw" => hw,
        "dQ" => 10.5,
        "dE" => 0.5,
        "W" => W,
        "nev_initial" => nev_initial,
        "nev_final" => nev_final,
        "temperature" => temperature,
        "volume" => volume,
        "npoints" => npoints,
        "Q_range" => [Q_min, Q_max],
        "Q0_crossing" => Q₀,
        "cutoff" => cut_off,
        "sigma" => σ
    ),
    "eigenvalues_initial" => pot_i.ϵ[1:min(20, length(pot_i.ϵ))],
    "eigenvalues_final" => pot_f.ϵ[1:min(20, length(pot_f.ϵ))],
    "capture_coefficient_300K" => C_300K,
    "julia_version" => string(VERSION),
    "carriercapture_version" => "CarrierCapture.jl (local)"
)

output_path = joinpath(@__DIR__, "reference_data", "sn_zn_julia_reference.json")
mkpath(dirname(output_path))
open(output_path, "w") do f
    JSON.print(f, results, 4)
end

println("  Saved to: $output_path")

println("\n" * "="^60)
println("Julia reference calculation complete!")
println("="^60)
