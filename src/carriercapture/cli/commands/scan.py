"""
Scan command: High-throughput parameter screening.

Performs systematic parameter sweeps over ΔQ, ΔE, and ℏω to screen
materials for carrier capture properties. Supports parallel execution.
"""

import click
from pathlib import Path
import numpy as np


@click.command(name="scan")
@click.option(
    "--dQ-min",
    type=float,
    required=True,
    help="Minimum ΔQ value (amu^0.5·Å)"
)
@click.option(
    "--dQ-max",
    type=float,
    required=True,
    help="Maximum ΔQ value (amu^0.5·Å)"
)
@click.option(
    "--dQ-points",
    type=int,
    required=True,
    help="Number of ΔQ points"
)
@click.option(
    "--dE-min",
    type=float,
    required=True,
    help="Minimum ΔE value (eV)"
)
@click.option(
    "--dE-max",
    type=float,
    required=True,
    help="Maximum ΔE value (eV)"
)
@click.option(
    "--dE-points",
    type=int,
    required=True,
    help="Number of ΔE points"
)
@click.option(
    "--hbar-omega-i",
    type=float,
    default=0.008,
    help="ℏω for initial state (eV). Default: 0.008 (8 meV)"
)
@click.option(
    "--hbar-omega-f",
    type=float,
    default=0.008,
    help="ℏω for final state (eV). Default: 0.008 (8 meV)"
)
@click.option(
    "-T", "--temperature",
    type=float,
    default=300.0,
    help="Temperature (K). Default: 300"
)
@click.option(
    "-V", "--volume",
    type=float,
    default=1e-21,
    help="Supercell volume (cm³). Default: 1e-21"
)
@click.option(
    "-g", "--degeneracy",
    type=int,
    default=1,
    help="Degeneracy factor. Default: 1"
)
@click.option(
    "--sigma",
    type=float,
    default=0.01,
    help="Gaussian delta width (eV). Default: 0.01"
)
@click.option(
    "--cutoff",
    type=float,
    default=0.25,
    help="Energy cutoff for overlaps (eV). Default: 0.25"
)
@click.option(
    "--nev-i",
    type=int,
    default=180,
    help="Number of eigenvalues for initial state. Default: 180"
)
@click.option(
    "--nev-f",
    type=int,
    default=60,
    help="Number of eigenvalues for final state. Default: 60"
)
@click.option(
    "-j", "--n-jobs",
    type=int,
    default=1,
    help="Number of parallel jobs. Use -1 for all cores. Default: 1"
)
@click.option(
    "-o", "--output",
    type=click.Path(path_type=Path),
    required=True,
    help="Output file (.npz or .h5)"
)
@click.option(
    "--no-progress",
    is_flag=True,
    help="Disable progress bar"
)
@click.pass_context
def scan_cmd(ctx, dq_min, dq_max, dq_points, de_min, de_max, de_points,
             hbar_omega_i, hbar_omega_f, temperature, volume, degeneracy,
             sigma, cutoff, nev_i, nev_f, n_jobs, output, no_progress):
    """
    Run high-throughput parameter scan.

    Systematically scans parameter space (ΔQ, ΔE) to compute capture
    coefficients for many material configurations. Uses harmonic
    potentials for fast screening.

    \\b
    Examples:
      # Basic scan over ΔQ and ΔE
      $ carriercapture scan --dQ-min 0 --dQ-max 25 --dQ-points 25 \\
                           --dQ-min 0 --dE-max 2.5 --dE-points 10 \\
                           -o scan_results.npz

      # Parallel scan with 4 cores
      $ carriercapture scan --dQ-min 0 --dQ-max 25 --dQ-points 25 \\
                           --dE-min 0 --dE-max 2.5 --dE-points 10 \\
                           -j 4 -o results.npz

      # Use all available cores
      $ carriercapture scan --dQ-min 0 --dQ-max 25 --dQ-points 50 \\
                           --dE-min 0 --dE-max 2.5 --dE-points 20 \\
                           -j -1 -o results.npz

      # Custom phonon frequencies
      $ carriercapture scan --dQ-min 0 --dQ-max 25 --dQ-points 25 \\
                           --dE-min 0 --dE-max 2.5 --dE-points 10 \\
                           --hbar-omega-i 0.010 --hbar-omega-f 0.010 \\
                           -o results.npz
    """
    verbose = ctx.obj.get('verbose', 0)

    if verbose > 0:
        click.echo("Setting up parameter scan...")
        click.echo(f"  ΔQ: {dq_points} points from {dq_min} to {dq_max} amu^0.5·Å")
        click.echo(f"  ΔE: {de_points} points from {de_min} to {de_max} eV")
        click.echo(f"  ℏω_i = {hbar_omega_i:.4f} eV, ℏω_f = {hbar_omega_f:.4f} eV")
        click.echo(f"  Temperature: {temperature} K")
        click.echo(f"  Volume: {volume:.2e} cm³")
        click.echo(f"  Parallel jobs: {n_jobs}")

    try:
        from carriercapture.analysis.parameter_scan import (
            ParameterScanner,
            ScanParameters
        )

        # Create scan parameters
        params = ScanParameters(
            dQ_range=(dq_min, dq_max, dq_points),
            dE_range=(de_min, de_max, de_points),
            hbar_omega_i=hbar_omega_i,
            hbar_omega_f=hbar_omega_f,
            temperature=temperature,
            volume=volume,
            degeneracy=degeneracy,
            sigma=sigma,
            cutoff=cutoff,
            nev_initial=nev_i,
            nev_final=nev_f,
        )

        # Create scanner
        scanner = ParameterScanner(params, verbose=verbose > 0)

        # Run scan
        if verbose > 0:
            click.echo("\nStarting scan...")

        results = scanner.run_harmonic_scan(
            n_jobs=n_jobs,
            show_progress=not no_progress
        )

        # Save results
        if verbose > 0:
            click.echo(f"\nSaving results to: {output}")

        # Detect format from extension
        if output.suffix.lower() in ['.h5', '.hdf5']:
            file_format = 'hdf5'
        else:
            file_format = 'npz'

        results.save(output, format=file_format)

        if verbose > 0:
            click.echo("✓ Scan complete!")
            click.echo(f"\nResults summary:")
            click.echo(f"  Grid size: {dq_points} × {de_points} = {dq_points * de_points} points")

            # Calculate statistics
            valid_mask = ~np.isnan(results.capture_coefficients)
            if np.any(valid_mask):
                valid_coeffs = results.capture_coefficients[valid_mask]
                click.echo(f"  Successful calculations: {np.sum(valid_mask)}/{results.capture_coefficients.size}")
                click.echo(f"  Capture coefficient range:")
                click.echo(f"    Min: {valid_coeffs.min():.3e} cm³/s")
                click.echo(f"    Max: {valid_coeffs.max():.3e} cm³/s")
                click.echo(f"    Mean: {valid_coeffs.mean():.3e} cm³/s")
            else:
                click.echo(f"  No successful calculations!")

    except ImportError as e:
        click.echo(
            f"Error: Required dependencies not available.\n"
            f"Install with: pip install joblib\n"
            f"Details: {e}",
            err=True
        )
        ctx.exit(1)
    except Exception as e:
        click.echo(f"Error during parameter scan: {e}", err=True)
        if verbose > 1:
            import traceback
            traceback.print_exc()
        ctx.exit(1)


@click.command(name="scan-plot")
@click.argument(
    "scan_file",
    type=click.Path(exists=True, path_type=Path),
)
@click.option(
    "--type",
    "plot_type",
    type=click.Choice(["heatmap", "contour", "both"]),
    default="heatmap",
    help="Type of plot to generate"
)
@click.option(
    "--log-scale",
    is_flag=True,
    help="Use log scale for capture coefficients"
)
@click.option(
    "-o", "--output",
    type=click.Path(path_type=Path),
    help="Output file (.html, .png)"
)
@click.option(
    "--show",
    is_flag=True,
    help="Display plot in browser"
)
@click.pass_context
def scan_plot_cmd(ctx, scan_file, plot_type, log_scale, output, show):
    """
    Plot results from parameter scan.

    Generates heatmaps or contour plots of capture coefficients
    as a function of ΔQ and ΔE.

    \\b
    Examples:
      # Plot heatmap
      $ carriercapture scan-plot results.npz --show

      # Plot with log scale
      $ carriercapture scan-plot results.npz --log-scale -o heatmap.html

      # Generate contour plot
      $ carriercapture scan-plot results.npz --type contour --show
    """
    verbose = ctx.obj.get('verbose', 0)

    if verbose > 0:
        click.echo(f"Loading scan results from: {scan_file}")

    try:
        from carriercapture.analysis.parameter_scan import ScanResult
        import plotly.graph_objects as go

        # Load results
        if scan_file.suffix.lower() in ['.h5', '.hdf5']:
            results = ScanResult.load(scan_file, format='hdf5')
        else:
            results = ScanResult.load(scan_file, format='npz')

        if verbose > 0:
            click.echo(f"Loaded grid: {results.dQ_grid.shape[0]} × {results.dE_grid.shape[0]}")

        # Prepare data
        Z = results.capture_coefficients
        if log_scale:
            Z = np.log10(Z + 1e-30)  # Add small epsilon to avoid log(0)
            colorbar_title = "log₁₀(C) [cm³/s]"
        else:
            colorbar_title = "C [cm³/s]"

        # Create figure
        fig = go.Figure()

        if plot_type in ["heatmap", "both"]:
            fig.add_trace(go.Heatmap(
                x=results.dE_grid,
                y=results.dQ_grid,
                z=Z,
                colorscale='Viridis',
                colorbar=dict(title=colorbar_title),
            ))

        if plot_type in ["contour", "both"]:
            fig.add_trace(go.Contour(
                x=results.dE_grid,
                y=results.dQ_grid,
                z=Z,
                colorscale='Viridis',
                colorbar=dict(title=colorbar_title),
            ))

        fig.update_layout(
            title="Parameter Scan: Capture Coefficient",
            xaxis_title="ΔE (eV)",
            yaxis_title="ΔQ (amu<sup>0.5</sup>·Å)",
            template="plotly_white",
            width=800,
            height=700,
        )

        # Save or show
        if output:
            fig.write_html(str(output))
            if verbose > 0:
                click.echo(f"✓ Saved plot to: {output}")

        if show:
            fig.show()

        if not output and not show:
            click.echo("Note: Use --output to save or --show to display plot")

    except Exception as e:
        click.echo(f"Error plotting scan results: {e}", err=True)
        if verbose > 1:
            import traceback
            traceback.print_exc()
        ctx.exit(1)
