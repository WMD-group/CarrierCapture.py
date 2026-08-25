"""
Physical constants used in carrier capture calculations.

All constants are in SI-compatible units suitable for semiconductor physics:
- Energy: eV (electronvolts)
- Length: Å (angstroms) or cm
- Mass: amu (atomic mass units)
- Temperature: K (kelvin)
"""

# Fundamental constants
AMU = 931.4940954e6  # eV/c² - atomic mass unit in energy units
HBAR_C = 0.19732697e-6  # eV·m - reduced Planck constant times speed of light
HBAR = 6.582119514e-16  # eV·s - reduced Planck constant
K_B = 8.6173303e-5  # eV/K - Boltzmann constant
M_E_C2 = 0.5109989461e6  # eV - electron rest energy m_e·c²
FINE_STRUCTURE = 7.2973525693e-3  # dimensionless fine-structure constant
RYDBERG_EV = 0.5 * M_E_C2 * FINE_STRUCTURE**2  # eV - Rydberg energy (13.6057 eV)

# Conversion factors
EV_TO_HARTREE = 1.0 / 27.21138602  # Hartree to eV conversion
HARTREE_TO_EV = 27.21138602

# Computational parameters
OCC_CUTOFF = 1e-5  # Convergence criterion for partition function
"""Maximum occupation of highest eigenvalue for convergence.

If occ(ε_max) >= OCC_CUTOFF, the partition function is not converged
and more eigenvalues (nev) should be computed.
"""
