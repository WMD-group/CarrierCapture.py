"""
Tests for the Sommerfeld (Coulomb enhancement) factor.
"""

import pytest
import numpy as np

from carriercapture import sommerfeld_parameter
from carriercapture._constants import K_B, RYDBERG_EV


class TestSommerfeldParameter:
    """Test sommerfeld_parameter behavior and limits."""

    def test_neutral_is_unity(self):
        """Z=0 gives exactly 1 for both methods, scalar and array."""
        for method in ("Analytic", "Integrate"):
            assert sommerfeld_parameter(300.0, Z=0, m_eff=0.2, eps0=10, method=method) == 1.0
            T = np.linspace(100, 500, 5)
            s = sommerfeld_parameter(T, Z=0, m_eff=0.2, eps0=10, method=method)
            np.testing.assert_array_equal(s, np.ones(5))

    def test_attractive_enhances_and_decreases_with_T(self):
        """Attractive center: s > 1 and monotonically decreasing with T."""
        T = np.linspace(50, 800, 20)
        for method in ("Analytic", "Integrate"):
            s = sommerfeld_parameter(T, Z=-1, m_eff=0.2, eps0=10, method=method)
            assert np.all(s > 1)
            assert np.all(np.diff(s) < 0)

    def test_repulsive_suppresses_and_increases_with_T(self):
        """Repulsive center: 0 < s < 1 and increasing with T."""
        T = np.linspace(50, 800, 20)
        for method in ("Analytic", "Integrate"):
            s = sommerfeld_parameter(T, Z=+1, m_eff=0.2, eps0=10, method=method)
            assert np.all(s > 0)
            assert np.all(s < 1)
            assert np.all(np.diff(s) > 0)

    def test_analytic_attractive_regression(self):
        """Hand-derived analytic value: s = 4*sqrt(theta/pi)."""
        T, Z, m_eff, eps0 = 300.0, -1, 0.2, 10.0
        E_R = m_eff * RYDBERG_EV / eps0**2
        theta = np.pi**2 * Z**2 * E_R / (K_B * T)
        expected = 4.0 * np.sqrt(theta / np.pi)
        s = sommerfeld_parameter(T, Z=Z, m_eff=m_eff, eps0=eps0, method="Analytic")
        assert s == pytest.approx(expected, rel=1e-12)
        assert s == pytest.approx(7.27, rel=1e-2)

    def test_analytic_matches_integrate_at_low_T(self):
        """The Pässler analytic form is the low-T limit of the integral."""
        s_a = sommerfeld_parameter(50.0, Z=-1, m_eff=0.2, eps0=10, method="Analytic")
        s_i = sommerfeld_parameter(50.0, Z=-1, m_eff=0.2, eps0=10, method="Integrate")
        assert s_a == pytest.approx(s_i, rel=1e-3)
        # The exact enhancement eta/(1 - exp(-eta)) >= eta, so the thermal
        # average always sits at or above the analytic (low-T) limit, with
        # the gap growing at high T.
        s_a_hot = sommerfeld_parameter(800.0, Z=-1, m_eff=0.2, eps0=10, method="Analytic")
        s_i_hot = sommerfeld_parameter(800.0, Z=-1, m_eff=0.2, eps0=10, method="Integrate")
        assert s_i_hot > s_a_hot
        assert (s_i_hot - s_a_hot) / s_a_hot > (s_i - s_a) / s_a

    def test_shapes_and_types(self):
        """Scalar in -> float out; array in -> same-shape array out."""
        s = sommerfeld_parameter(300.0, Z=-1, m_eff=0.2, eps0=10)
        assert np.isscalar(s) or s.ndim == 0
        T = np.linspace(100, 500, 7)
        s_arr = sommerfeld_parameter(T, Z=-1, m_eff=0.2, eps0=10)
        assert s_arr.shape == T.shape

    def test_unknown_method_raises(self):
        with pytest.raises(ValueError, match="method"):
            sommerfeld_parameter(300.0, Z=-1, m_eff=0.2, eps0=10, method="bogus")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
