"""
Sommerfeld (Coulomb enhancement) factor for capture by charged defects.

The capture coefficients computed by ConfigCoordinate assume a neutral
defect. For a charged defect, the Coulomb interaction between the free
carrier and the defect enhances (attractive) or suppresses (repulsive)
the carrier density at the defect site. The corrected coefficient is:

    C_charged(T) = s(T) * C(T)

following Pässler, phys. stat. sol. (b) 78, 625 (1976) and
Alkauskas et al., Phys. Rev. B 90, 075202 (2014), Sec. II.F.
"""

from typing import Union

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import quad

from .._constants import K_B, RYDBERG_EV


def sommerfeld_parameter(
    temperature: Union[float, NDArray[np.float64]],
    Z: int,
    m_eff: float,
    eps0: float,
    method: str = "Integrate",
) -> Union[float, NDArray[np.float64]]:
    """
    Calculate the Sommerfeld factor s(T) for capture by a charged defect.

    Parameters
    ----------
    temperature : float or NDArray[np.float64]
        Temperature(s) (K)
    Z : int
        Defect charge divided by carrier charge: Z < 0 for attractive
        centers, Z > 0 for repulsive, Z = 0 returns 1
    m_eff : float
        Carrier effective mass (units of the electron mass)
    eps0 : float
        Relative static dielectric constant
    method : str, default="Integrate"
        "Integrate": Maxwell-Boltzmann thermal average of the exact
        Coulomb enhancement factor. "Analytic": low-temperature limits
        of Pässler (1976).

    Returns
    -------
    s : float or NDArray[np.float64]
        Sommerfeld factor (same shape as temperature)

    Notes
    -----
    With the scaled Rydberg E_R = m_eff·Ry/eps0² and θ = π²Z²E_R/(k_B·T),
    the analytic limits are:

        attractive:  s = 4·sqrt(θ/π)
        repulsive:   s = (8/sqrt(3))·θ^(2/3)·exp(-3·θ^(1/3))

    The integral form averages f(E) = η/(1 - exp(-η)) (attractive) or
    η/(exp(η) - 1) (repulsive), with η(E) = 2π|Z|·sqrt(E_R/E), over a
    Maxwell-Boltzmann distribution.

    Examples
    --------
    >>> sommerfeld_parameter(300, Z=-1, m_eff=0.2, eps0=10, method="Analytic")
    7.27...
    """
    T = np.asarray(temperature, dtype=float)
    scalar_input = T.ndim == 0
    T = np.atleast_1d(T)

    if Z == 0:
        s = np.ones_like(T)
        return s[0] if scalar_input else s

    E_R = m_eff * RYDBERG_EV / eps0**2  # scaled Rydberg (eV)

    if method.lower().startswith("a"):
        theta = np.pi**2 * Z**2 * E_R / (K_B * T)
        if Z < 0:
            s = 4.0 * np.sqrt(theta / np.pi)
        else:
            s = (8.0 / np.sqrt(3.0)) * theta ** (2.0 / 3.0) * np.exp(-3.0 * theta ** (1.0 / 3.0))
    elif method.lower().startswith("i"):
        # s(T) = <f(E)> over MB distribution, in x = E/(k_B T):
        #   s = (2/sqrt(pi)) * ∫ f(x·kT) √x e^(-x) dx
        def enhancement(x: float, kT: float) -> float:
            eta = 2.0 * np.pi * abs(Z) * np.sqrt(E_R / (x * kT))
            if Z < 0:
                return eta / -np.expm1(-eta)
            return eta / np.expm1(eta)

        s = np.array(
            [
                (2.0 / np.sqrt(np.pi))
                * quad(lambda x: enhancement(x, K_B * t) * np.sqrt(x) * np.exp(-x), 0, np.inf)[0]
                for t in T
            ]
        )
    else:
        raise ValueError(f"Unknown method '{method}': use 'Integrate' or 'Analytic'")

    return s[0] if scalar_input else s
