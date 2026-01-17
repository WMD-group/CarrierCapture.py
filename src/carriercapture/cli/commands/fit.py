"""
Fit command: Fit potential energy surface to data.

Supports multiple fitting methods:
- spline: Cubic spline interpolation (FITPACK)
- harmonic: Harmonic oscillator
- morse: Morse potential
- polynomial: Polynomial of specified degree
- morse_poly: Morse + polynomial correction
"""

import click
from pathlib import Path
import numpy as np

from carriercapture.core import Potential
from carriercapture.io.readers import read_potential_data, load_potential_from_file
from carriercapture.io.writers import save_potential


@click.command(name="fit")
@click.argument(
    "data_file",
    type=click.Path(exists=True, path_type=Path),
)
@click.option(
    "-f", "--fit-type",
    type=click.Choice(["spline", "harmonic", "morse", "polynomial", "morse_poly"]),
    default="spline",
    help="Fitting method"
)
@click.option(
    "-o", "--output",
    type=click.Path(path_type=Path),
    help="Output file (auto-detect format from extension: .json, .yaml, .npz, .dat)"
)
@click.option(
    "--order",
    type=int,
    default=4,
    help="Spline order (for spline fitting)"
)
@click.option(
    "--smoothness",
    type=float,
    default=0.001,
    help="Spline smoothness parameter (for spline fitting)"
)
@click.option(
    "--degree",
    type=int,
    default=4,
    help="Polynomial degree (for polynomial fitting)"
)
@click.option(
    "--hw",
    type=float,
    help="Phonon frequency in eV (for harmonic fitting)"
)
@click.option(
    "--Q0",
    type=float,
    default=0.0,
    help="Equilibrium position (amu^0.5·Å)"
)
@click.option(
    "--E0",
    type=float,
    default=0.0,
    help="Energy offset (eV)"
)
@click.option(
    "--plot",
    is_flag=True,
    help="Generate plot of fitted potential"
)
@click.option(
    "--plot-output",
    type=click.Path(path_type=Path),
    help="Save plot to file (requires --plot)"
)
@click.pass_context
def fit_cmd(ctx, data_file, fit_type, output, order, smoothness, degree, hw, q0, e0, plot, plot_output):
    """
    Fit potential energy surface to data.

    \b
    DATA_FILE: Path to potential data file (Q, E format)

    \b
    Examples:
      # Spline fit with custom parameters
      $ carriercapture fit excited.dat -f spline --order 4 --smoothness 0.001

      # Harmonic fit
      $ carriercapture fit data.dat -f harmonic --hw 0.03 --Q0 5.0

      # Morse potential fit
      $ carriercapture fit data.dat -f morse -o morse_fit.json

      # Polynomial fit
      $ carriercapture fit data.dat -f polynomial --degree 6

      # Fit and plot
      $ carriercapture fit data.dat -f spline --plot --plot-output fit.png
    """
    verbose = ctx.obj.get('verbose', 0)

    # Read data
    if verbose > 0:
        click.echo(f"Reading data from: {data_file}")

    Q_data, E_data = read_potential_data(data_file)

    if verbose > 0:
        click.echo(f"Loaded {len(Q_data)} data points")
        click.echo(f"Q range: [{Q_data.min():.3f}, {Q_data.max():.3f}] amu^0.5·Å")
        click.echo(f"E range: [{E_data.min():.3f}, {E_data.max():.3f}] eV")

    # Create potential
    pot = Potential(Q_data=Q_data, E_data=E_data, Q0=q0, E0=e0)

    # Prepare fit parameters
    fit_kwargs = {}
    if fit_type == "spline":
        fit_kwargs = {"order": order, "smoothness": smoothness}
    elif fit_type == "harmonic":
        if hw is None:
            click.echo("Error: --hw (phonon frequency) required for harmonic fitting", err=True)
            ctx.exit(1)
        fit_kwargs = {"hw": hw}
    elif fit_type == "polynomial":
        fit_kwargs = {"degree": degree}
    elif fit_type == "morse_poly":
        fit_kwargs = {"poly_degree": degree}

    # Fit
    if verbose > 0:
        click.echo(f"\nFitting with {fit_type} method...")
        if fit_kwargs:
            params_str = ", ".join(f"{k}={v}" for k, v in fit_kwargs.items())
            click.echo(f"Parameters: {params_str}")

    try:
        pot.fit(fit_type=fit_type, **fit_kwargs)
    except Exception as e:
        click.echo(f"Error during fitting: {e}", err=True)
        ctx.exit(1)

    if verbose > 0:
        click.echo("✓ Fitting completed successfully")
        if pot.fit_params:
            click.echo("\nFit parameters:")
            for key, val in pot.fit_params.items():
                if isinstance(val, (list, np.ndarray)) and len(val) > 5:
                    click.echo(f"  {key}: [{len(val)} values]")
                else:
                    click.echo(f"  {key}: {val}")

    # Save output
    if output:
        if verbose > 0:
            click.echo(f"\nSaving to: {output}")

        try:
            save_potential(pot, output)
            if verbose > 0:
                click.echo("✓ Saved successfully")
        except Exception as e:
            click.echo(f"Error saving output: {e}", err=True)
            ctx.exit(1)
    else:
        if verbose == 0:
            # If no output file and not verbose, print basic info
            click.echo(f"Fitted {len(Q_data)} points with {fit_type} method")

    # Plot if requested
    if plot:
        try:
            import matplotlib.pyplot as plt

            plt.figure(figsize=(10, 6))

            # Original data
            plt.scatter(Q_data, E_data, alpha=0.5, label="Data", s=30)

            # Fitted curve
            Q_fit = np.linspace(Q_data.min(), Q_data.max(), 500)
            E_fit = pot(Q_fit)
            plt.plot(Q_fit, E_fit, 'r-', linewidth=2, label=f"{fit_type.capitalize()} fit")

            plt.xlabel("Q (amu$^{0.5}$·Å)", fontsize=12)
            plt.ylabel("E (eV)", fontsize=12)
            plt.title(f"Potential Fit: {fit_type}", fontsize=14)
            plt.legend()
            plt.grid(True, alpha=0.3)
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
