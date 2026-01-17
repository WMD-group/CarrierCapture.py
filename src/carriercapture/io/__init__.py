"""I/O utilities for reading and writing data."""

from .readers import (
    read_potential_data,
    read_json,
    read_yaml,
    read_csv,
    read_npz,
    load_potential_from_file,
)
from .writers import (
    write_json,
    write_yaml,
    write_csv,
    write_npz,
    write_potential_data,
    write_capture_results,
    save_potential,
)

# Optional doped integration (requires doped package)
try:
    from .doped_interface import (
        load_defect_entry,
        get_available_charge_states,
        validate_charge_states,
        suggest_Q0,
        load_path_calculations,
        extract_cc_data_from_structures,
        create_potential_from_doped,
    )
    DOPED_INTEGRATION_AVAILABLE = True
except ImportError:
    DOPED_INTEGRATION_AVAILABLE = False

__all__ = [
    # Readers
    "read_potential_data",
    "read_json",
    "read_yaml",
    "read_csv",
    "read_npz",
    "load_potential_from_file",
    # Writers
    "write_json",
    "write_yaml",
    "write_csv",
    "write_npz",
    "write_potential_data",
    "write_capture_results",
    "save_potential",
]

# Add doped functions to __all__ if available
if DOPED_INTEGRATION_AVAILABLE:
    __all__.extend([
        # doped integration
        "load_defect_entry",
        "get_available_charge_states",
        "validate_charge_states",
        "suggest_Q0",
        "load_path_calculations",
        "extract_cc_data_from_structures",
        "create_potential_from_doped",
    ])
