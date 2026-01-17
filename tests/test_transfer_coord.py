"""
Tests for TransferCoordinate class.

Validates Marcus theory calculations for charge transfer rates and mobility.
"""

import pytest
import numpy as np
from carriercapture.core import Potential, TransferCoordinate


class TestTransferCoordinateCreation:
    """Test TransferCoordinate initialization."""

    def test_initialization(self):
        """Test basic initialization."""
        pot_1 = Potential.from_harmonic(hw=0.02, Q0=0.0, E0=0.0)
        pot_2 = Potential.from_harmonic(hw=0.02, Q0=5.0, E0=0.1)

        tc = TransferCoordinate(pot_1, pot_2, name="test_transfer")

        assert tc.name == "test_transfer"
        assert tc.pot_1 is pot_1
        assert tc.pot_2 is pot_2
        assert tc.Q_cross is None
        assert tc.E_cross is None
        assert tc.coupling is None
        assert tc.reorganization_energy is None
        assert tc.activation_energy is None
        assert tc.transfer_rate is None
        assert tc.temperature is None


class TestCouplingCalculation:
    """Test electronic coupling calculation."""

    def test_calculate_coupling_basic(self):
        """Test basic coupling calculation."""
        # Two harmonic potentials that cross
        pot_1 = Potential.from_harmonic(hw=0.03, Q0=0.0, E0=1.0)
        pot_2 = Potential.from_harmonic(hw=0.03, Q0=8.0, E0=0.0)

        tc = TransferCoordinate(pot_1, pot_2)

        # Calculate coupling (auto-detect crossing)
        coupling = tc.get_coupling()

        assert coupling is not None
        assert coupling > 0
        assert tc.Q_cross is not None
        assert tc.E_cross is not None

        # Crossing should be around midpoint
        assert tc.Q_cross == pytest.approx(4.0, abs=1.0)

    def test_calculate_coupling_specified_crossing(self):
        """Test coupling with specified crossing point."""
        pot_1 = Potential.from_harmonic(hw=0.03, Q0=0.0, E0=1.0)
        pot_2 = Potential.from_harmonic(hw=0.03, Q0=10.0, E0=0.0)

        tc = TransferCoordinate(pot_1, pot_2)

        # Specify crossing point manually
        Q_cross_manual = 5.0
        coupling = tc.get_coupling(Q_cross=Q_cross_manual)

        assert coupling > 0
        assert tc.Q_cross == Q_cross_manual

        # Energy should be evaluated at specified point
        E1 = pot_1(Q_cross_manual)
        E2 = pot_2(Q_cross_manual)
        assert tc.E_cross == pytest.approx(0.5 * (E1 + E2), rel=1e-6)

    def test_coupling_without_fit_raises(self):
        """Test that coupling calculation requires fitted potentials."""
        pot_1 = Potential(Q_data=np.linspace(0, 10, 50), E_data=np.linspace(0, 5, 50))
        pot_2 = Potential(Q_data=np.linspace(0, 10, 50), E_data=np.linspace(1, 4, 50))

        tc = TransferCoordinate(pot_1, pot_2)

        with pytest.raises(ValueError, match="must be fitted"):
            tc.get_coupling()


class TestReorganizationEnergy:
    """Test reorganization energy calculation."""

    def test_reorganization_energy_symmetric(self):
        """Test reorganization energy for symmetric case."""
        # Two identical parabolas, displaced
        hw = 0.02  # eV
        dQ = 8.0  # amu^0.5·Å

        pot_1 = Potential.from_harmonic(hw=hw, Q0=0.0, E0=0.0)
        pot_2 = Potential.from_harmonic(hw=hw, Q0=dQ, E0=0.0)

        tc = TransferCoordinate(pot_1, pot_2)

        lambda_reorg = tc.get_reorganization_energy()

        # For harmonic potentials: λ = 0.5 * k * dQ^2 = 0.5 * m * ω^2 * dQ^2
        # But our harmonic uses: E = a * (Q - Q0)^2 where a = (amu/2) * (hw/(hbar_c*1e10))^2
        # So λ = a * dQ^2

        # Let's just check it's positive and reasonable
        assert lambda_reorg > 0
        assert lambda_reorg < 1.0  # Should be less than 1 eV for typical parameters

    def test_reorganization_energy_asymmetric(self):
        """Test reorganization energy for asymmetric case."""
        # Different curvatures
        pot_1 = Potential.from_harmonic(hw=0.03, Q0=0.0, E0=0.0)
        pot_2 = Potential.from_harmonic(hw=0.02, Q0=10.0, E0=0.5)

        tc = TransferCoordinate(pot_1, pot_2)

        lambda_reorg = tc.get_reorganization_energy()

        assert lambda_reorg > 0
        # Should be larger due to larger displacement
        assert lambda_reorg > 0.1

    def test_reorganization_energy_consistency(self):
        """Test that reorganization energy is symmetric."""
        pot_1 = Potential.from_harmonic(hw=0.025, Q0=0.0, E0=0.0)
        pot_2 = Potential.from_harmonic(hw=0.025, Q0=6.0, E0=0.0)

        tc = TransferCoordinate(pot_1, pot_2)
        lambda_reorg = tc.get_reorganization_energy()

        # Reverse the potentials
        tc_rev = TransferCoordinate(pot_2, pot_1)
        lambda_reorg_rev = tc_rev.get_reorganization_energy()

        # Should be the same (symmetric)
        assert lambda_reorg == pytest.approx(lambda_reorg_rev, rel=0.01)

    def test_reorganization_energy_without_fit_raises(self):
        """Test that reorganization energy requires fitted potentials."""
        pot_1 = Potential(Q_data=np.linspace(0, 10, 50), E_data=np.linspace(0, 5, 50))
        pot_2 = Potential(Q_data=np.linspace(0, 10, 50), E_data=np.linspace(1, 4, 50))

        tc = TransferCoordinate(pot_1, pot_2)

        with pytest.raises(ValueError, match="must be fitted"):
            tc.get_reorganization_energy()


class TestActivationEnergy:
    """Test activation energy calculation."""

    def test_activation_energy_symmetric(self):
        """Test activation energy for symmetric case (ΔG = 0)."""
        pot_1 = Potential.from_harmonic(hw=0.02, Q0=0.0, E0=0.0)
        pot_2 = Potential.from_harmonic(hw=0.02, Q0=8.0, E0=0.0)

        tc = TransferCoordinate(pot_1, pot_2)
        lambda_reorg = tc.get_reorganization_energy()

        # Symmetric case: ΔG‡ = λ / 4
        barrier = tc.get_activation_energy(delta_G=0.0)

        expected = lambda_reorg / 4.0
        assert barrier == pytest.approx(expected, rel=1e-6)

    def test_activation_energy_downhill(self):
        """Test activation energy for exergonic (downhill) transfer."""
        pot_1 = Potential.from_harmonic(hw=0.02, Q0=0.0, E0=0.5)
        pot_2 = Potential.from_harmonic(hw=0.02, Q0=8.0, E0=0.0)

        tc = TransferCoordinate(pot_1, pot_2)
        lambda_reorg = tc.get_reorganization_energy()

        # Downhill transfer (ΔG < 0)
        delta_G = -0.2  # 0.2 eV downhill
        barrier = tc.get_activation_energy(delta_G=delta_G)

        # Should be less than symmetric case
        barrier_symmetric = lambda_reorg / 4.0
        assert barrier < barrier_symmetric

        # Marcus formula: ΔG‡ = (λ + ΔG)² / (4λ)
        expected = (lambda_reorg + delta_G) ** 2 / (4 * lambda_reorg)
        assert barrier == pytest.approx(expected, rel=1e-6)

    def test_activation_energy_uphill(self):
        """Test activation energy for endergonic (uphill) transfer."""
        pot_1 = Potential.from_harmonic(hw=0.02, Q0=0.0, E0=0.0)
        pot_2 = Potential.from_harmonic(hw=0.02, Q0=8.0, E0=0.3)

        tc = TransferCoordinate(pot_1, pot_2)
        lambda_reorg = tc.get_reorganization_energy()

        # Uphill transfer (ΔG > 0)
        delta_G = 0.2  # 0.2 eV uphill
        barrier = tc.get_activation_energy(delta_G=delta_G)

        # Should be greater than symmetric case
        barrier_symmetric = lambda_reorg / 4.0
        assert barrier > barrier_symmetric

    def test_activation_energy_without_reorganization_raises(self):
        """Test that activation energy requires reorganization energy."""
        pot_1 = Potential.from_harmonic(hw=0.02, Q0=0.0, E0=0.0)
        pot_2 = Potential.from_harmonic(hw=0.02, Q0=5.0, E0=0.0)

        tc = TransferCoordinate(pot_1, pot_2)

        with pytest.raises(ValueError, match="reorganization energy"):
            tc.get_activation_energy()


class TestTransferRate:
    """Test Marcus transfer rate calculation."""

    def test_transfer_rate_basic(self):
        """Test basic transfer rate calculation."""
        pot_1 = Potential.from_harmonic(hw=0.02, Q0=0.0, E0=0.0)
        pot_2 = Potential.from_harmonic(hw=0.02, Q0=8.0, E0=0.0)

        tc = TransferCoordinate(pot_1, pot_2)

        # Calculate prerequisites
        tc.get_coupling()
        tc.get_reorganization_energy()

        # Calculate transfer rate
        temperature = np.linspace(100, 500, 50)
        rate = tc.get_transfer_rate(temperature=temperature)

        assert rate is not None
        assert len(rate) == 50
        assert np.all(rate > 0)
        assert tc.temperature is not None
        assert tc.transfer_rate is not None

    def test_transfer_rate_temperature_dependence(self):
        """Test that transfer rate varies with temperature."""
        pot_1 = Potential.from_harmonic(hw=0.02, Q0=0.0, E0=0.0)
        pot_2 = Potential.from_harmonic(hw=0.02, Q0=10.0, E0=0.0)

        tc = TransferCoordinate(pot_1, pot_2)
        tc.get_coupling()
        tc.get_reorganization_energy()

        temperature = np.array([100.0, 200.0, 300.0, 400.0, 500.0])
        rate = tc.get_transfer_rate(temperature=temperature)

        # Rate should vary with temperature
        assert len(np.unique(rate)) > 1
        assert np.std(rate) > 0

        # Generally, rate increases with temperature for activated process
        # (though Marcus theory has a maximum at T = λ / (4kB))
        # At low T, should see increase
        assert rate[1] > rate[0]

    def test_transfer_rate_downhill(self):
        """Test transfer rate for downhill (exergonic) transfer."""
        pot_1 = Potential.from_harmonic(hw=0.02, Q0=0.0, E0=0.5)
        pot_2 = Potential.from_harmonic(hw=0.02, Q0=8.0, E0=0.0)

        tc = TransferCoordinate(pot_1, pot_2)
        tc.get_coupling()
        tc.get_reorganization_energy()

        temperature = np.array([300.0])

        # Symmetric case (ΔG = 0)
        rate_symmetric = tc.get_transfer_rate(temperature=temperature, delta_G=0.0)

        # Downhill case (ΔG < 0)
        rate_downhill = tc.get_transfer_rate(temperature=temperature, delta_G=-0.2)

        # Downhill should be faster
        assert rate_downhill[0] > rate_symmetric[0]

    def test_transfer_rate_without_prerequisites_raises(self):
        """Test that transfer rate requires coupling and reorganization energy."""
        pot_1 = Potential.from_harmonic(hw=0.02, Q0=0.0, E0=0.0)
        pot_2 = Potential.from_harmonic(hw=0.02, Q0=5.0, E0=0.0)

        tc = TransferCoordinate(pot_1, pot_2)

        temperature = np.array([300.0])

        # Without coupling
        with pytest.raises(ValueError, match="coupling"):
            tc.get_transfer_rate(temperature=temperature)

        # With coupling but without reorganization energy
        tc.get_coupling()
        with pytest.raises(ValueError, match="reorganization energy"):
            tc.get_transfer_rate(temperature=temperature)


class TestMobility:
    """Test Einstein mobility calculation."""

    def test_mobility_basic(self):
        """Test basic mobility calculation."""
        pot_1 = Potential.from_harmonic(hw=0.02, Q0=0.0, E0=0.0)
        pot_2 = Potential.from_harmonic(hw=0.02, Q0=8.0, E0=0.0)

        tc = TransferCoordinate(pot_1, pot_2)
        tc.get_coupling()
        tc.get_reorganization_energy()

        temperature = np.linspace(100, 500, 20)
        distance = 5.0  # Å (typical hopping distance)

        mobility = tc.calculate_mobility(
            temperature=temperature, distance=distance, delta_G=0.0
        )

        assert mobility is not None
        assert len(mobility) == 20
        assert np.all(mobility > 0)

    def test_mobility_temperature_dependence(self):
        """Test that mobility varies with temperature."""
        pot_1 = Potential.from_harmonic(hw=0.02, Q0=0.0, E0=0.0)
        pot_2 = Potential.from_harmonic(hw=0.02, Q0=10.0, E0=0.0)

        tc = TransferCoordinate(pot_1, pot_2)
        tc.get_coupling()
        tc.get_reorganization_energy()

        temperature = np.array([100.0, 200.0, 300.0, 400.0, 500.0])
        mobility = tc.calculate_mobility(
            temperature=temperature, distance=5.0, delta_G=0.0
        )

        # Mobility should vary with temperature
        assert len(np.unique(mobility)) > 1
        assert np.std(mobility) > 0

    def test_mobility_distance_dependence(self):
        """Test that mobility depends on hopping distance."""
        pot_1 = Potential.from_harmonic(hw=0.02, Q0=0.0, E0=0.0)
        pot_2 = Potential.from_harmonic(hw=0.02, Q0=8.0, E0=0.0)

        tc = TransferCoordinate(pot_1, pot_2)
        tc.get_coupling()
        tc.get_reorganization_energy()

        temperature = np.array([300.0])

        # Short hopping distance
        mobility_short = tc.calculate_mobility(
            temperature=temperature, distance=3.0, delta_G=0.0
        )

        # Long hopping distance
        mobility_long = tc.calculate_mobility(
            temperature=temperature, distance=8.0, delta_G=0.0
        )

        # Longer distance should give higher mobility (μ ∝ d²)
        assert mobility_long[0] > mobility_short[0]

        # Should scale as d²
        ratio = mobility_long[0] / mobility_short[0]
        expected_ratio = (8.0 / 3.0) ** 2
        assert ratio == pytest.approx(expected_ratio, rel=0.01)


class TestSerialization:
    """Test serialization and deserialization."""

    def test_to_dict_from_dict(self):
        """Test round-trip serialization."""
        pot_1 = Potential.from_harmonic(hw=0.02, Q0=0.0, E0=0.0)
        pot_2 = Potential.from_harmonic(hw=0.02, Q0=8.0, E0=0.1)

        tc = TransferCoordinate(pot_1, pot_2, name="test_transfer")

        # Calculate some properties
        tc.get_coupling()
        tc.get_reorganization_energy()
        tc.get_activation_energy(delta_G=0.0)
        tc.get_transfer_rate(temperature=np.linspace(100, 500, 20))

        # Serialize
        data = tc.to_dict()

        assert data["name"] == "test_transfer"
        assert "pot_1" in data
        assert "pot_2" in data
        assert "Q_cross" in data
        assert "coupling" in data
        assert "reorganization_energy" in data
        assert "transfer_rate" in data

        # Deserialize
        tc2 = TransferCoordinate.from_dict(data)

        assert tc2.name == tc.name
        assert tc2.Q_cross == pytest.approx(tc.Q_cross, rel=1e-6)
        assert tc2.coupling == pytest.approx(tc.coupling, rel=1e-6)
        assert tc2.reorganization_energy == pytest.approx(tc.reorganization_energy, rel=1e-6)
        np.testing.assert_allclose(tc2.transfer_rate, tc.transfer_rate, rtol=1e-6)


class TestFullWorkflow:
    """Test complete Marcus theory workflow."""

    def test_simple_charge_transfer(self):
        """Test full workflow for simple charge transfer."""
        # Two harmonic states with slight energy offset
        hw = 0.02  # eV
        dQ = 10.0  # amu^0.5·Å
        dE = 0.1  # eV energy offset

        pot_1 = Potential.from_harmonic(hw=hw, Q0=0.0, E0=dE)
        pot_2 = Potential.from_harmonic(hw=hw, Q0=dQ, E0=0.0)

        # Create transfer coordinate
        tc = TransferCoordinate(pot_1, pot_2, name="hole_transfer")

        # Calculate properties in order
        coupling = tc.get_coupling()
        lambda_reorg = tc.get_reorganization_energy()
        barrier = tc.get_activation_energy(delta_G=dE)

        # Calculate transfer rate
        temperature = np.linspace(100, 500, 50)
        rate = tc.get_transfer_rate(temperature=temperature, delta_G=dE)

        # Calculate mobility
        mobility = tc.calculate_mobility(
            temperature=temperature, distance=5.0, delta_G=dE
        )

        # Sanity checks
        assert coupling > 0
        assert lambda_reorg > 0
        assert barrier > 0
        assert np.all(rate > 0)
        assert np.all(mobility > 0)
        assert np.all(np.isfinite(rate))
        assert np.all(np.isfinite(mobility))

        # Physical reasonableness
        assert coupling < 0.5  # Typical coupling is << 1 eV
        assert lambda_reorg < 2.0  # Typical reorganization < few eV
        assert barrier < lambda_reorg  # Barrier should be < λ
        assert np.max(rate) < 1e15  # Transfer rate < phonon frequency

    def test_symmetric_vs_asymmetric_transfer(self):
        """Compare symmetric and asymmetric transfer rates."""
        # Symmetric case
        pot_1_sym = Potential.from_harmonic(hw=0.02, Q0=0.0, E0=0.0)
        pot_2_sym = Potential.from_harmonic(hw=0.02, Q0=8.0, E0=0.0)

        tc_sym = TransferCoordinate(pot_1_sym, pot_2_sym)
        tc_sym.get_coupling()
        tc_sym.get_reorganization_energy()

        temperature = np.array([300.0])
        rate_sym = tc_sym.get_transfer_rate(temperature=temperature, delta_G=0.0)

        # Asymmetric case (downhill)
        pot_1_asym = Potential.from_harmonic(hw=0.02, Q0=0.0, E0=0.3)
        pot_2_asym = Potential.from_harmonic(hw=0.02, Q0=8.0, E0=0.0)

        tc_asym = TransferCoordinate(pot_1_asym, pot_2_asym)
        tc_asym.get_coupling()
        tc_asym.get_reorganization_energy()

        rate_asym = tc_asym.get_transfer_rate(temperature=temperature, delta_G=0.3)

        # Downhill should be faster
        assert rate_asym[0] > rate_sym[0]


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
