# Core API

Core classes for potential energy surfaces, capture coefficient calculations, and quantum mechanics.

## Potential

The `Potential` class represents 1D potential energy surfaces with fitting and quantum solving capabilities.

::: carriercapture.core.potential.Potential
    options:
      show_root_heading: true
      show_source: false
      members_order: source
      heading_level: 3

## ConfigCoordinate

The `ConfigCoordinate` class calculates carrier capture coefficients using configuration coordinate diagrams.

::: carriercapture.core.config_coord.ConfigCoordinate
    options:
      show_root_heading: true
      show_source: false
      members_order: source
      heading_level: 3

## TransferCoordinate

!!! warning "Experimental"
    The `TransferCoordinate` class implements Marcus theory for charge transfer. This is an experimental feature (Phase 3) and may have incomplete functionality.

::: carriercapture.core.transfer_coord.TransferCoordinate
    options:
      show_root_heading: true
      show_source: false
      members_order: source
      heading_level: 3

## Schrödinger Solver

Low-level functions for solving the 1D Schrödinger equation.

::: carriercapture.core.schrodinger
    options:
      show_root_heading: true
      show_source: false
      heading_level: 3

## Constants

Physical constants used throughout CarrierCapture.

::: carriercapture._constants
    options:
      show_root_heading: true
      show_source: false
      heading_level: 3
