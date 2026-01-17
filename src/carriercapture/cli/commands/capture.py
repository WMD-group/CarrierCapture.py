"""
Capture command: Calculate carrier capture coefficients.

Computes multiphonon carrier capture rates from configuration coordinate diagram.
"""

import click
from pathlib import Path
import numpy as np

from carriercapture.core import Potential, ConfigCoordinate
from carriercapture.io.readers import read_yaml, load_potential_from_file
from carriercapture.io.writers import write_capture_results


@click.command(name="capture")
@click.argument(
    "config_file",
    type=click.Path(exists=True, path_type=Path),
    required=False,
)
@click.option(
    "--pot-i",
    type=click.Path(exists=True, path_type=Path),
    help="Initial state potential file"
)
@click.option(
    "--pot-f",
    type=click.Path(exists=True, path_type=Path),
    help="Final state potential file"
)
@click.option(
    "-W", "--coupling",
    type=float,
    help="Electron-phonon coupling (eV)"
)
@click.option(
    "-g", "--degeneracy",
    type=int,
    default=1,
    help="Degeneracy factor"
)
@click.option(
    "-V", "--volume",
    type=float,
    help="Supercell volume (cm³)"
)
@click.option(
    "--temp-range",
    nargs=3,
    type=float,
    help="Temperature range: T_min T_max n_points (K)"
)
@click.option(
    "--Q0",
    type=float,
    help="Shift for coordinate operator (amu^0.5·Å)"
)
@click.option(
    "--cutoff",
    type=float,
    default=0.25,
    help="Energy cutoff for overlaps (eV)"
)
@click.option(
    "--sigma",
    type=float,
    default=0.025,
    help="Gaussian delta width (eV)"
)
@click.option(
    "-o", "--output",
    type=click.Path(path_type=Path),
    help="Output file (.json, .yaml, .csv, .npz)"
)
@click.option(
    "--plot",
    is_flag=True,
    help="Generate Arrhenius plot"
)
@click.option(
    "--plot-output",
    type=click.Path(path_type=Path),
    help="Save plot to file"
)
# doped integration options
@click.option(
    "--doped",
    type=click.Path(exists=True, path_type=Path),
    help="Path to doped DefectEntry JSON.GZ file (enables doped integration mode)"
)
@click.option(
    "--charge-i",
    type=int,
    help="Initial charge state (for doped mode)"
)
@click.option(
    "--charge-f",
    type=int,
    help="Final charge state (for doped mode)"
)
@click.option(
    "--doped-path-i",
    type=click.Path(exists=True, path_type=Path),
    help="Path directory for initial state VASP calculations (for doped mode)"
)
@click.option(
    "--doped-path-f",
    type=click.Path(exists=True, path_type=Path),
    help="Path directory for final state VASP calculations (for doped mode)"
)
@click.option(
    "--n-images",
    type=int,
    default=10,
    help="Number of interpolation images (for doped structure generation)"
)
@click.option(
    "--auto-Q0",
    is_flag=True,
    help="Automatically suggest Q0 from structure displacement (doped mode)"
)
@click.pass_context
def capture_cmd(ctx, config_file, pot_i, pot_f, coupling, degeneracy, volume,
                temp_range, q0, cutoff, sigma, output, plot, plot_output,
                doped, charge_i, charge_f, doped_path_i, doped_path_f, n_images, auto_q0):
    """
    Calculate carrier capture coefficients.

    Can use either:
    1. YAML config file with all parameters
    2. Command-line options for quick calculations
    3. doped integration mode (with --doped flag)

    \b
    Examples:
      # From config file
      $ carriercapture capture config.yaml -o results.json

      # From command line
      $ carriercapture capture --pot-i excited.json --pot-f ground.json \\
          -W 0.205 -V 1e-21 --temp-range 100 500 50 --Q0 10.0

      # With plotting
      $ carriercapture capture config.yaml --plot --plot-output arrhenius.png

      # doped integration mode (from VASP path calculations)
      $ carriercapture capture --doped defect.json.gz \\
          --charge-i 0 --charge-f +1 \\
          --doped-path-i path_q0/ --doped-path-f path_q1/ \\
          -W 0.205 -V 1e-21 --temp-range 100 500 50 \\
          --auto-Q0 -o results.json
    """
    verbose = ctx.obj.get('verbose', 0)

    # Load from config file if provided
    if config_file:
        if verbose > 0:
            click.echo(f"Loading configuration from: {config_file}")

        try:
            config = read_yaml(config_file)
        except Exception as e:
            click.echo(f"Error reading config file: {e}", err=True)
            ctx.exit(1)

        # Extract parameters from config
        pot_i_config = config.get('potential_initial', {})
        pot_f_config = config.get('potential_final', {})
        capture_config = config.get('capture', {})

        # Override with command-line options if provided
        if not pot_i:
            pot_i = pot_i_config.get('file')
        if not pot_f:
            pot_f = pot_f_config.get('file')
        if not coupling:
            coupling = capture_config.get('W')
        if degeneracy == 1:  # Default value
            degeneracy = capture_config.get('degeneracy', 1)
        if not volume:
            volume = capture_config.get('volume')
        if not temp_range:
            temp_config = capture_config.get('temperature', {})
            if isinstance(temp_config, dict):
                t_min = temp_config.get('min', 100)
                t_max = temp_config.get('max', 500)
                n_points = temp_config.get('n_points', 50)
                temp_range = (t_min, t_max, n_points)
        if not q0:
            q0 = capture_config.get('Q0')
        if cutoff == 0.25:  # Default value
            cutoff = capture_config.get('cutoff', 0.25)
        if sigma == 0.025:  # Default value
            sigma = capture_config.get('sigma', 0.025)

    # Handle doped integration mode
    if doped:
        if verbose > 0:
            click.echo(f"\n=== doped Integration Mode ===")
            click.echo(f"Loading DefectEntry from: {doped}")

        # Check if doped integration is available
        try:
            from carriercapture.io.doped_interface import (
                load_defect_entry,
                load_path_calculations,
                extract_cc_data_from_structures,
                create_potential_from_doped,
                suggest_Q0 as doped_suggest_Q0,
            )
        except ImportError:
            click.echo(
                "Error: doped integration requires doped package. "
                "Install with: pip install carriercapture[doped]",
                err=True
            )
            ctx.exit(1)

        # Validate charge states provided
        if charge_i is None or charge_f is None:
            click.echo(
                "Error: --charge-i and --charge-f required in doped mode",
                err=True
            )
            ctx.exit(1)

        # Load DefectEntry
        try:
            defect_entry = load_defect_entry(doped)
            if verbose > 0:
                click.echo(f"✓ Loaded DefectEntry: {getattr(defect_entry, 'name', 'unnamed')}")
                click.echo(f"  Charge states: initial={charge_i:+d}, final={charge_f:+d}")
        except Exception as e:
            click.echo(f"Error loading DefectEntry: {e}", err=True)
            ctx.exit(1)

        # Load Q-E data from doped
        # Option 1: Load from VASP path calculations
        if doped_path_i and doped_path_f:
            if verbose > 0:
                click.echo(f"\nLoading VASP path calculations...")
                click.echo(f"  Initial state: {doped_path_i}")
                click.echo(f"  Final state: {doped_path_f}")

            try:
                Q_i, E_i = load_path_calculations(doped_path_i, verbose=verbose > 1)
                Q_f, E_f = load_path_calculations(doped_path_f, verbose=verbose > 1)

                if verbose > 0:
                    click.echo(f"✓ Loaded {len(Q_i)} initial state images")
                    click.echo(f"✓ Loaded {len(Q_f)} final state images")

            except Exception as e:
                click.echo(f"Error loading path calculations: {e}", err=True)
                ctx.exit(1)

        # Option 2: Extract from structures (simple harmonic approximation)
        else:
            if verbose > 0:
                click.echo(f"\nExtracting Q-E data from DefectEntry structures...")
                click.echo(f"  Using harmonic approximation with {n_images} interpolation points")

            # Note: This requires DefectEntry to have structure information
            # For now, print a helpful message
            click.echo(
                "Warning: Full structure-based CC data extraction not yet implemented. "
                "Please provide VASP path calculations with --doped-path-i and --doped-path-f",
                err=True
            )
            click.echo(
                "Alternatively, use the standard workflow with --pot-i and --pot-f",
                err=True
            )
            ctx.exit(1)

        # Create Potential objects from doped data
        if verbose > 0:
            click.echo(f"\nCreating Potential objects...")

        try:
            potential_i = create_potential_from_doped(
                defect_entry, charge_i, Q_data=Q_i, E_data=E_i
            )
            potential_f = create_potential_from_doped(
                defect_entry, charge_f, Q_data=Q_f, E_data=E_f
            )

            # Fit and solve potentials
            if verbose > 0:
                click.echo(f"  Fitting potentials with spline...")

            potential_i.fit(fit_type="spline", order=4, smoothness=0.001)
            potential_f.fit(fit_type="spline", order=4, smoothness=0.001)

            if verbose > 0:
                click.echo(f"  Solving Schrödinger equation...")

            # Use higher nev for initial state if not specified
            nev_i = 180 if charge_i == 0 else 60
            nev_f = 60

            potential_i.solve(nev=nev_i)
            potential_f.solve(nev=nev_f)

            if verbose > 0:
                click.echo(f"✓ Initial state: {len(potential_i.eigenvalues)} states")
                click.echo(f"✓ Final state: {len(potential_f.eigenvalues)} states")

        except Exception as e:
            click.echo(f"Error creating potentials from doped data: {e}", err=True)
            ctx.exit(1)

        # Auto-suggest Q0 if requested
        if auto_q0 and q0 is None:
            try:
                # Get structures from defect_entry
                # This is a placeholder - actual implementation depends on DefectEntry structure
                if verbose > 0:
                    click.echo(f"\nAuto-calculating Q0 from structure displacement...")

                # For now, use midpoint of Q range
                q0 = 0.5 * (Q_i[-1] + Q_f[-1])

                if verbose > 0:
                    click.echo(f"✓ Suggested Q0: {q0:.2f} amu^0.5·Å")
            except Exception as e:
                if verbose > 0:
                    click.echo(f"Warning: Could not auto-calculate Q0: {e}")

        # Set pot_i and pot_f to signal that we have potentials loaded
        pot_i = "doped"
        pot_f = "doped"

    # Validate required parameters
    if not pot_i or not pot_f:
        click.echo("Error: Must provide --pot-i and --pot-f (or config file, or --doped)", err=True)
        ctx.exit(1)
    if coupling is None:
        click.echo("Error: Must provide -W/--coupling (or config file)", err=True)
        ctx.exit(1)
    if volume is None:
        click.echo("Error: Must provide -V/--volume (or config file)", err=True)
        ctx.exit(1)
    if q0 is None:
        click.echo("Error: Must provide --Q0 (or config file)", err=True)
        ctx.exit(1)
    if not temp_range:
        temp_range = (100, 500, 50)  # Default

    # Load potentials (skip if already loaded from doped)
    if not doped:
        if verbose > 0:
            click.echo(f"\nLoading initial potential: {pot_i}")

        try:
            data_i = load_potential_from_file(pot_i)
            potential_i = Potential.from_dict(data_i)
        except Exception as e:
            click.echo(f"Error loading initial potential: {e}", err=True)
            ctx.exit(1)

        if verbose > 0:
            click.echo(f"Loading final potential: {pot_f}")

        try:
            data_f = load_potential_from_file(pot_f)
            potential_f = Potential.from_dict(data_f)
        except Exception as e:
            click.echo(f"Error loading final potential: {e}", err=True)
            ctx.exit(1)

        # Check that potentials are solved
        if potential_i.eigenvalues is None:
            click.echo("Error: Initial potential must be solved (use 'solve' command)", err=True)
            ctx.exit(1)
        if potential_f.eigenvalues is None:
            click.echo("Error: Final potential must be solved (use 'solve' command)", err=True)
            ctx.exit(1)

    if verbose > 0:
        click.echo(f"\nInitial potential: {len(potential_i.eigenvalues)} states")
        click.echo(f"Final potential: {len(potential_f.eigenvalues)} states")
        click.echo(f"\nCapture parameters:")
        click.echo(f"  W (coupling): {coupling} eV")
        click.echo(f"  g (degeneracy): {degeneracy}")
        click.echo(f"  V (volume): {volume:.2e} cm³")
        click.echo(f"  Q0: {q0} amu^0.5·Å")
        click.echo(f"  Energy cutoff: {cutoff} eV")
        click.echo(f"  Delta width: {sigma} eV")

    # Create configuration coordinate
    cc = ConfigCoordinate(
        pot_i=potential_i,
        pot_f=potential_f,
        W=coupling,
        degeneracy=degeneracy,
    )

    # Calculate overlaps
    if verbose > 0:
        click.echo(f"\nCalculating wavefunction overlaps...")

    try:
        cc.calculate_overlap(Q0=q0, cutoff=cutoff, sigma=sigma)
    except Exception as e:
        click.echo(f"Error calculating overlaps: {e}", err=True)
        ctx.exit(1)

    if verbose > 0:
        n_nonzero = np.count_nonzero(cc.overlap_matrix)
        n_total = cc.overlap_matrix.size
        click.echo(f"✓ Computed {n_nonzero}/{n_total} non-zero overlaps "
                  f"({100*n_nonzero/n_total:.1f}%)")

    # Calculate capture coefficient
    if verbose > 0:
        click.echo(f"\nCalculating capture coefficient...")
        click.echo(f"  Temperature: {temp_range[0]:.0f} - {temp_range[1]:.0f} K "
                  f"({int(temp_range[2])} points)")

    temperature = np.linspace(temp_range[0], temp_range[1], int(temp_range[2]))

    try:
        cc.calculate_capture_coefficient(volume=volume, temperature=temperature)
    except Exception as e:
        click.echo(f"Error calculating capture coefficient: {e}", err=True)
        ctx.exit(1)

    if verbose > 0:
        click.echo("✓ Calculation completed successfully")
        click.echo(f"\nResults:")
        click.echo(f"  C(T={temperature[0]:.0f}K) = {cc.capture_coefficient[0]:.3e} cm³/s")
        click.echo(f"  C(T={temperature[len(temperature)//2]:.0f}K) = "
                  f"{cc.capture_coefficient[len(temperature)//2]:.3e} cm³/s")
        click.echo(f"  C(T={temperature[-1]:.0f}K) = {cc.capture_coefficient[-1]:.3e} cm³/s")

    # Save output
    if output:
        if verbose > 0:
            click.echo(f"\nSaving results to: {output}")

        try:
            # Detect format from extension
            ext = output.suffix.lower()
            format_map = {
                '.json': 'json',
                '.yaml': 'yaml',
                '.yml': 'yaml',
                '.csv': 'csv',
                '.npz': 'npz',
            }
            file_format = format_map.get(ext, 'json')

            write_capture_results(cc, output, file_format=file_format)
            if verbose > 0:
                click.echo("✓ Saved successfully")
        except Exception as e:
            click.echo(f"Error saving results: {e}", err=True)
            ctx.exit(1)
    else:
        # Print summary if no output file
        click.echo(f"\nCapture coefficient calculated for {len(temperature)} temperatures")
        click.echo(f"Use -o/--output to save results")

    # Plot if requested
    if plot:
        try:
            import matplotlib.pyplot as plt

            fig, ax = plt.subplots(figsize=(10, 6))

            # Arrhenius plot: log(C) vs 1000/T
            x = 1000.0 / temperature
            y = np.log10(cc.capture_coefficient)

            ax.plot(x, y, 'o-', linewidth=2, markersize=5)

            ax.set_xlabel("1000/T (K$^{-1}$)", fontsize=12)
            ax.set_ylabel("log$_{10}$(C) [cm$^3$/s]", fontsize=12)
            ax.set_title("Capture Coefficient (Arrhenius Plot)", fontsize=14)
            ax.grid(True, alpha=0.3)

            # Add temperature labels on top axis
            ax2 = ax.twiny()
            temps_label = [100, 200, 300, 400, 500]
            ax2.set_xlim(ax.get_xlim())
            ax2.set_xticks([1000/t for t in temps_label if temp_range[0] <= t <= temp_range[1]])
            ax2.set_xticklabels([f"{t}K" for t in temps_label if temp_range[0] <= t <= temp_range[1]])

            plt.tight_layout()

            if plot_output:
                plt.savefig(plot_output, dpi=300, bbox_inches='tight')
                if verbose > 0:
                    click.echo(f"✓ Plot saved to: {plot_output}")
            else:
                plt.show()

        except ImportError:
            click.echo("Warning: matplotlib not available for plotting", err=True)
        except Exception as e:
            click.echo(f"Error during plotting: {e}", err=True)
