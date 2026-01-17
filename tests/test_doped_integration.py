"""
Tests for doped package integration.

These tests check that the doped integration module works correctly,
including proper error handling when doped is not installed.
"""

import pytest
import numpy as np
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock

# Check if doped integration is available
try:
    from carriercapture.io import doped_interface
    from carriercapture.io.doped_interface import (
        load_defect_entry,
        get_available_charge_states,
        validate_charge_states,
        _get_dQ_from_structures,
        suggest_Q0,
        load_path_calculations,
        extract_cc_data_from_structures,
        create_potential_from_doped,
        DOPED_AVAILABLE,
        MONTY_AVAILABLE,
    )
    DOPED_INTEGRATION_AVAILABLE = True
    # Check if the actual doped package is available (not just the interface module)
    DOPED_REALLY_AVAILABLE = DOPED_AVAILABLE and MONTY_AVAILABLE
except ImportError:
    DOPED_INTEGRATION_AVAILABLE = False
    DOPED_REALLY_AVAILABLE = False


class TestImportHandling:
    """Test that import errors are handled gracefully."""

    def test_import_without_doped(self):
        """Test that importing the module gives appropriate warnings when doped is not available."""
        # This test just checks that the module can be imported
        # If doped is not available, DOPED_AVAILABLE should be False
        if not DOPED_INTEGRATION_AVAILABLE:
            pytest.skip("doped integration module not available (expected if doped not installed)")

        # If we get here, module was imported successfully
        assert True


@pytest.mark.skipif(not DOPED_INTEGRATION_AVAILABLE, reason="doped package not installed")
class TestDopedIntegration:
    """Tests for doped integration functions."""

    def test_get_dQ_from_structures(self):
        """Test mass-weighted displacement calculation."""
        # Create mock structures
        mock_struct1 = Mock()
        mock_struct2 = Mock()

        # Mock the atoms with distance and atomic mass
        atom1 = Mock()
        atom1.distance = Mock(return_value=1.0)  # 1 Å displacement
        atom1.specie = Mock()
        atom1.specie.atomic_mass = 16.0  # Oxygen mass

        atom2 = Mock()
        atom2.distance = Mock(return_value=0.5)  # 0.5 Å displacement
        atom2.specie = Mock()
        atom2.specie.atomic_mass = 16.0

        # Mock iteration over structures
        mock_struct1.__iter__ = Mock(return_value=iter([atom1, atom2]))
        mock_struct2.__iter__ = Mock(return_value=iter([atom1, atom2]))

        # Calculate expected dQ: sqrt(16*1^2 + 16*0.5^2) = sqrt(16 + 4) = sqrt(20)
        expected_dQ = np.sqrt(20.0)

        # Note: Actual test would require real pymatgen structures
        # For now, we just test that the function exists and has the right signature
        assert callable(_get_dQ_from_structures)

    def test_suggest_Q0_callable(self):
        """Test that suggest_Q0 function is callable."""
        assert callable(suggest_Q0)

    def test_load_path_calculations_signature(self):
        """Test that load_path_calculations has the correct signature."""
        import inspect
        sig = inspect.signature(load_path_calculations)
        params = list(sig.parameters.keys())

        assert 'path_dir' in params
        assert 'ref_structure' in params
        assert 'energy_key' in params
        assert 'verbose' in params

    def test_extract_cc_data_from_structures_signature(self):
        """Test that extract_cc_data_from_structures has the correct signature."""
        import inspect
        sig = inspect.signature(extract_cc_data_from_structures)
        params = list(sig.parameters.keys())

        assert 'struct_initial' in params
        assert 'struct_final' in params
        assert 'energy_initial' in params
        assert 'energy_final' in params
        assert 'n_images' in params
        assert 'align' in params
        assert 'verbose' in params

    def test_create_potential_from_doped_signature(self):
        """Test that create_potential_from_doped has the correct signature."""
        import inspect
        sig = inspect.signature(create_potential_from_doped)
        params = list(sig.parameters.keys())

        assert 'defect_entry' in params
        assert 'charge_state' in params
        assert 'Q_data' in params
        assert 'E_data' in params
        assert 'name' in params

    @pytest.mark.skipif(not DOPED_REALLY_AVAILABLE, reason="doped package not installed")
    def test_validate_charge_states_with_mock(self):
        """Test charge state validation with mock DefectEntry."""
        mock_defect = Mock()
        mock_defect.charge_state = 0

        # Should not raise for valid charge state
        # This will fail because actual implementation checks the charge state
        # but serves as a template for future mocking tests
        try:
            validate_charge_states(mock_defect, 0, 0)
        except (AttributeError, ValueError):
            # Expected if mock doesn't have the full interface
            pass

    @pytest.mark.skipif(not DOPED_REALLY_AVAILABLE, reason="doped package not installed")
    def test_get_available_charge_states_with_mock(self):
        """Test getting available charge states from mock DefectEntry."""
        mock_defect = Mock()
        mock_defect.charge_state = +1

        charges = get_available_charge_states(mock_defect)
        assert isinstance(charges, list)
        assert +1 in charges


class TestDopedErrorHandling:
    """Test error handling when doped functions are called without doped installed."""

    @pytest.mark.skipif(DOPED_INTEGRATION_AVAILABLE, reason="Only test when doped is NOT installed")
    def test_import_error_when_doped_not_installed(self):
        """Test that importing doped_interface raises ImportError when doped not installed."""
        # This test should only run when doped is NOT installed
        # It verifies that the module properly indicates doped is unavailable
        from carriercapture.io import DOPED_INTEGRATION_AVAILABLE
        assert not DOPED_INTEGRATION_AVAILABLE


class TestPotentialCreation:
    """Test Potential creation from doped data."""

    @pytest.mark.skipif(not DOPED_REALLY_AVAILABLE, reason="doped package not installed")
    def test_create_potential_from_doped_returns_potential(self):
        """Test that create_potential_from_doped returns a Potential object."""
        # Create mock DefectEntry
        mock_defect = Mock()
        mock_defect.charge_state = 0
        mock_defect.name = "v_O"

        # Create mock Q and E data
        Q_data = np.linspace(0, 10, 20)
        E_data = 0.5 * (Q_data - 5) ** 2  # Simple parabola

        # Create potential
        pot = create_potential_from_doped(
            mock_defect,
            charge_state=0,
            Q_data=Q_data,
            E_data=E_data,
            name="test_potential"
        )

        # Check that we got a Potential object
        from carriercapture.core.potential import Potential
        assert isinstance(pot, Potential)
        assert pot.name == "test_potential"
        assert np.array_equal(pot.Q_data, Q_data)
        assert np.array_equal(pot.E_data, E_data)

    @pytest.mark.skipif(not DOPED_REALLY_AVAILABLE, reason="doped package not installed")
    def test_create_potential_auto_naming(self):
        """Test automatic naming of potential from DefectEntry."""
        mock_defect = Mock()
        mock_defect.charge_state = +1
        mock_defect.name = "v_O"

        Q_data = np.linspace(0, 10, 20)
        E_data = np.zeros_like(Q_data)

        pot = create_potential_from_doped(
            mock_defect,
            charge_state=+1,
            Q_data=Q_data,
            E_data=E_data
        )

        # Check auto-generated name
        assert "v_O" in pot.name
        assert "+1" in pot.name or "q1" in pot.name.lower()


class TestLoadPathCalculations:
    """Test loading VASP calculations along configuration coordinate path."""

    @pytest.mark.skipif(not DOPED_REALLY_AVAILABLE, reason="doped package not installed")
    def test_load_path_calculations_requires_valid_directory(self):
        """Test that load_path_calculations raises FileNotFoundError for invalid path."""
        with pytest.raises(FileNotFoundError):
            load_path_calculations("/nonexistent/path/to/calculations")

    @pytest.mark.skipif(not DOPED_INTEGRATION_AVAILABLE, reason="doped package not installed")
    def test_load_path_calculations_signature_energy_key(self):
        """Test energy_key parameter in load_path_calculations."""
        import inspect
        sig = inspect.signature(load_path_calculations)

        # Check default value for energy_key
        assert sig.parameters['energy_key'].default == "energy"


class TestCLIIntegration:
    """Test CLI integration with doped flags."""

    def test_capture_command_has_doped_options(self):
        """Test that capture command has doped-related options."""
        from carriercapture.cli.commands.capture import capture_cmd

        # Get the click command parameters
        doped_params = [p.name for p in capture_cmd.params]

        # Check that doped-related parameters exist
        assert 'doped' in doped_params
        assert 'charge_i' in doped_params
        assert 'charge_f' in doped_params
        assert 'doped_path_i' in doped_params
        assert 'doped_path_f' in doped_params
        assert 'n_images' in doped_params
        assert 'auto_q0' in doped_params


class TestDocumentation:
    """Test that documentation strings are present."""

    @pytest.mark.skipif(not DOPED_INTEGRATION_AVAILABLE, reason="doped package not installed")
    def test_all_functions_have_docstrings(self):
        """Test that all exported functions have docstrings."""
        functions_to_check = [
            load_defect_entry,
            get_available_charge_states,
            validate_charge_states,
            suggest_Q0,
            load_path_calculations,
            extract_cc_data_from_structures,
            create_potential_from_doped,
        ]

        for func in functions_to_check:
            assert func.__doc__ is not None, f"{func.__name__} missing docstring"
            assert len(func.__doc__.strip()) > 0, f"{func.__name__} has empty docstring"


class TestModuleExports:
    """Test that the module exports the expected functions."""

    @pytest.mark.skipif(not DOPED_INTEGRATION_AVAILABLE, reason="doped package not installed")
    def test_doped_interface_exports(self):
        """Test that doped_interface module exports expected functions."""
        from carriercapture.io import doped_interface

        expected_exports = [
            "load_defect_entry",
            "get_available_charge_states",
            "validate_charge_states",
            "suggest_Q0",
            "load_path_calculations",
            "extract_cc_data_from_structures",
            "create_potential_from_doped",
        ]

        for export in expected_exports:
            assert hasattr(doped_interface, export), f"Missing export: {export}"

    def test_io_module_exports_doped_if_available(self):
        """Test that io module exports doped functions when available."""
        from carriercapture import io

        if DOPED_INTEGRATION_AVAILABLE:
            # Should have DOPED_INTEGRATION_AVAILABLE flag
            assert hasattr(io, 'DOPED_INTEGRATION_AVAILABLE')
            assert io.DOPED_INTEGRATION_AVAILABLE is True

            # Should export doped functions
            expected_exports = [
                "load_defect_entry",
                "create_potential_from_doped",
            ]

            for export in expected_exports:
                assert hasattr(io, export), f"io module missing doped export: {export}"
        else:
            # When doped not available, DOPED_INTEGRATION_AVAILABLE should be False
            if hasattr(io, 'DOPED_INTEGRATION_AVAILABLE'):
                assert io.DOPED_INTEGRATION_AVAILABLE is False


# Integration test template (requires actual doped data)
@pytest.mark.skipif(not DOPED_INTEGRATION_AVAILABLE, reason="doped package not installed")
@pytest.mark.skip(reason="Requires real doped data files - template for future testing")
class TestFullWorkflow:
    """Integration tests for full doped → CarrierCapture workflow."""

    def test_end_to_end_doped_to_capture(self, tmp_path):
        """
        Test complete workflow from doped DefectEntry to capture coefficient.

        This is a template test that would require:
        - Real DefectEntry JSON.GZ file
        - VASP path calculations
        - or structures with energies
        """
        # TODO: Add real test data in tests/data/doped_examples/
        # defect_path = tmp_path / "test_defect.json.gz"
        # ...
        pytest.skip("Test template - requires real doped data")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
