# Basic Concepts

Understanding the key concepts behind CarrierCapture will help you use it effectively.

## Overview

CarrierCapture calculates **carrier capture rates** for defects in semiconductors using **multiphonon theory** (Huang-Rhys formalism). The key idea: charge carriers can be captured non-radiatively via emission of multiple phonons.

## Key Concepts

### 1. Potential Energy Surfaces (PES)

A **potential energy surface** describes how the total energy of a system varies with atomic configuration.

**Mathematical form:**

$$E(Q) = \text{energy as function of configuration coordinate } Q$$

Where $Q$ is the **mass-weighted displacement** coordinate (units: amu$^{0.5}$·Å).

**Example - Harmonic potential:**

$$E(Q) = E_0 + \frac{1}{2}k(Q - Q_0)^2$$

Where:

- $E_0$: minimum energy
- $Q_0$: equilibrium position
- $k$: force constant (related to phonon frequency)

In CarrierCapture:

```python
from carriercapture.core import Potential

# Create harmonic potential
pot = Potential.from_harmonic(
    hw=0.008,  # ℏω phonon energy (eV)
    Q0=5.0,    # Equilibrium at Q = 5
    E0=0.0     # Ground state energy
)
```

### 2. Configuration Coordinate Diagrams

A **configuration coordinate (CC) diagram** plots the potential energy surfaces for different electronic states on the same axes.

**Physical picture:**

- **Horizontal axis**: Configuration coordinate $Q$ (atomic positions)
- **Vertical axis**: Total energy $E$
- **Curves**: PES for different charge states

**Example:**

```
E (eV) │
       │     Initial state (excited)
   1.5 │    ╱‾‾‾╲
       │   ╱     ╲
   1.0 │  ╱       ╲
       │ ╱         ╲
   0.5 │╱           ╲___
       │              Final state (ground)
   0.0 │              ╲___╱
       └────────────────────────────── Q
            Q₀        Q₁
```

**Key features:**

1. **Displacement**: $\Delta Q = Q_1 - Q_0$ (shift between minima)
2. **Energy difference**: $\Delta E = E_{\text{initial}}(Q_1) - E_{\text{final}}(Q_1)$
3. **Crossing point**: Where curves intersect (barrier for capture)

### 3. Vibrational Wavefunctions

At each electronic state, the atoms vibrate around equilibrium. These vibrations are **quantized** with discrete energy levels.

**Schrödinger equation:**

$$\hat{H}\psi_n(Q) = \varepsilon_n \psi_n(Q)$$

Where:

- $\psi_n(Q)$: wavefunction for state $n$
- $\varepsilon_n$: energy eigenvalue
- $n = 0, 1, 2, ...$: quantum number

**Harmonic oscillator solutions:**

$$\varepsilon_n = E_0 + \hbar\omega(n + \tfrac{1}{2})$$

In CarrierCapture:

```python
# Solve for wavefunctions
pot.solve(nev=60)  # Get first 60 states

# Access results
energies = pot.eigenvalues   # Shape: (60,)
wavefunctions = pot.eigenvectors  # Shape: (60, N_Q)
```

### 4. Franck-Condon Principle

**Statement**: Electronic transitions are vertical on a CC diagram (atoms don't move during transition).

**Implication**: After a sudden electronic transition, atoms are in a non-equilibrium configuration and must relax via:

1. **Radiative**: Photon emission
2. **Non-radiative**: Phonon emission (what CarrierCapture calculates)

### 5. Carrier Capture Coefficient

The **capture coefficient** $C(T)$ quantifies how fast charge carriers are captured (units: cm³/s).

**Physical meaning:**

$$\text{Capture rate} = C(T) \times n$$

Where $n$ is the carrier concentration.

**Temperature dependence:**

- **Low T**: Thermally activated (exponential increase with T)
- **High T**: Multiphonon cascade (can decrease with T)

**Typical values:**

| Capture Type | $C$ (cm³/s) | Example |
|--------------|-------------|---------|
| Very slow | 10⁻²⁰ | Deep traps at low T |
| Slow | 10⁻¹⁵ | Shallow defects |
| Moderate | 10⁻¹² | Typical defects at 300K |
| Fast | 10⁻⁸ | Resonant capture |

### 6. Electron-Phonon Coupling

The **coupling strength** $W$ describes how strongly electronic and vibrational degrees of freedom interact.

**Physical meaning:**

- **Large $W$**: Strong coupling → fast non-radiative capture
- **Small $W$**: Weak coupling → slow capture (more radiative)

**Typical range**: $W \sim 0.01$ to $0.5$ eV

**In multiphonon theory:**

$$C \propto W^2$$

So capture rate scales quadratically with coupling!

### 7. Overlap Integrals

The **overlap** between initial and final wavefunctions determines transition probability.

**Mathematical form:**

$$S_{ij} = \langle \psi_i | \hat{O} | \psi_j \rangle = \int \psi_i^*(Q) \cdot \hat{O} \cdot \psi_j(Q) \, dQ$$

Where $\hat{O}$ is typically $(Q - Q_0)$ for dipole coupling.

**Physical interpretation:**

- **Large overlap**: Easy transition
- **Small overlap**: Forbidden/unlikely transition

**Selection rules**: Dictate which $i \to j$ transitions are allowed.

## Putting It Together: The Full Calculation

### Step-by-Step Process

1. **Define potentials** for initial and final electronic states
2. **Solve Schrödinger equation** to get vibrational states
3. **Calculate overlaps** between vibrational wavefunctions
4. **Apply energy conservation** (delta function)
5. **Thermally weight** contributions from different initial states
6. **Sum over all pathways** to get total capture coefficient

### Mathematical Formula

The capture coefficient at temperature $T$ is:

$$C(T) = \frac{V \cdot 2\pi}{\hbar} \cdot g \cdot W^2 \cdot \sum_{i,j} p_i |S_{ij}|^2 \delta(\varepsilon_i - \varepsilon_j)$$

Where:

- $V$: supercell volume
- $g$: degeneracy factor (often = 1)
- $W$: electron-phonon coupling
- $p_i = \frac{e^{-\beta\varepsilon_i}}{Z}$: thermal occupation (Boltzmann)
- $S_{ij}$: overlap integral
- $\delta(\Delta E)$: energy-conserving delta function (Gaussian)

### In Code

```python
from carriercapture.core import Potential, ConfigCoordinate

# 1. Define potentials
pot_i = Potential.from_harmonic(hw=0.008, Q0=0.0, E0=0.5)
pot_f = Potential.from_harmonic(hw=0.008, Q0=10.0, E0=0.0)

# 2. Solve for states
pot_i.solve(nev=180)
pot_f.solve(nev=60)

# 3-6. Calculate capture coefficient
cc = ConfigCoordinate(pot_i=pot_i, pot_f=pot_f, W=0.068)
cc.calculate_overlap(Q0=5.0)
cc.calculate_capture_coefficient(volume=1e-21, temperature=[300])

# Result
print(f"C(300K) = {cc.capture_coefficient[0]:.3e} cm³/s")
```

## Common Terminology

| Term | Symbol | Units | Meaning |
|------|--------|-------|---------|
| Configuration coordinate | $Q$ | amu$^{0.5}$·Å | Mass-weighted displacement |
| Phonon energy | $\hbar\omega$ | eV | Vibrational quantum |
| Displacement | $\Delta Q$ | amu$^{0.5}$·Å | Shift between states |
| Energy difference | $\Delta E$ | eV | Vertical separation |
| Coupling | $W$ | eV | Electron-phonon interaction |
| Capture coefficient | $C$ | cm³/s | Capture rate constant |
| Huang-Rhys factor | $S$ | - | Effective phonon number |
| Reorganization energy | $\lambda$ | eV | Energy to relax |

## Physical Intuition

### Why Multiphonon?

For large $\Delta E \gg \hbar\omega$, a single phonon can't bridge the gap. Instead:

1. **Cascade emission**: Many phonons emitted sequentially
2. **Effective phonon number**: $n_{\text{eff}} \sim \Delta E / \hbar\omega$

Example: $\Delta E = 1$ eV, $\hbar\omega = 0.01$ eV $\Rightarrow$ need ~100 phonons!

### Temperature Effects

**Low Temperature:**

- Only ground state occupied ($p_0 \approx 1$)
- Needs thermal activation to reach crossing point
- $C \propto e^{-E_a/k_BT}$

**High Temperature:**

- Many initial states populated
- More pathways available
- But overlap factors decrease (orthogonality)
- Complex temperature dependence

### Material Dependence

**Good light emitters** (LEDs):

- Small $\Delta Q$ → small overlap
- Weak coupling $W$ → slow non-radiative
- $\Rightarrow$ Radiative process wins!

**Efficient carrier capture** (detectors, solar cells):

- Large $\Delta Q$ → large overlap
- Strong coupling $W$ → fast non-radiative
- $\Rightarrow$ Carriers captured quickly!

## Assumptions and Limitations

### Static Coupling Approximation

**Assumption**: Potential surfaces are independent of electronic state.

**Limitation**: Breaks down for strongly correlated electrons.

### Harmonic Approximation

**Assumption**: Potentials are parabolic near minimum.

**Limitation**: Anharmonic effects matter for large displacements.

CarrierCapture can use **spline** or **Morse** potentials to go beyond harmonic!

### Single Configuration Coordinate

**Assumption**: One effective coordinate $Q$ captures all relaxation.

**Limitation**: Real systems have many vibrational modes.

**Justification**: Often one "breathing mode" dominates.

## Next Steps

Now that you understand the concepts:

- **[First Calculation](first-calculation.md)** - Detailed walkthrough
- **[User Guide: Potentials](../user-guide/potentials.md)** - Learn about potential types
- **[Theory: Multiphonon](../theory/multiphonon-theory.md)** - Deep dive into theory
- **[Tutorials](../tutorials/index.md)** - Hands-on examples

## Further Reading

- **Alkauskas et al.** (2014) - [Phys. Rev. B 90, 075202](https://doi.org/10.1103/PhysRevB.90.075202)
- **Huang & Rhys** (1950) - [Theory of Light Absorption and Non-Radiative Transitions in F-Centres, Proc. R. Soc. Lond. A **204**, 406](https://royalsocietypublishing.org/rspa/article/204/1078/406/8369/Theory-of-light-absorption-and-non-radiative)
- **Stoneham (1975)** - "Theory of Defects in Solids"
