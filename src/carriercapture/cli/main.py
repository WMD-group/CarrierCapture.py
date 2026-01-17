"""
Main CLI entry point for CarrierCapture.

Provides a Click-based command-line interface with subcommands for
fitting, solving, and analyzing carrier capture rates.
"""

import click
from pathlib import Path


@click.group()
@click.version_option()
@click.option(
    "-v", "--verbose",
    count=True,
    help="Increase verbosity (can be repeated: -v, -vv, -vvv)"
)
@click.pass_context
def cli(ctx, verbose):
    """
    CarrierCapture: Carrier capture and non-radiative recombination calculations.

    A modern Python package for computing carrier capture rates in semiconductors
    using multiphonon theory.

    \b
    Common workflows:
      1. Fit potential:     carriercapture fit data.dat -f spline -o potential.json
      2. Solve Schrödinger: carriercapture solve potential.json -n 60
      3. Calculate capture: carriercapture capture config.yaml -o results.json
      4. Parameter scan:    carriercapture scan --dQ-min 0 --dQ-max 25 ...
      5. Visualize:         carriercapture viz --port 8050
      6. Full workflow:     carriercapture capture config.yaml --plot

    \b
    Examples:
      # Fit potential energy surface
      $ carriercapture fit excited.dat -f spline --order 4 --smoothness 0.001

      # Solve for phonon states
      $ carriercapture solve potential.json -n 180 -o solved.json

      # Calculate capture coefficient
      $ carriercapture capture config.yaml -V 1e-21 --temp-range 100 500 50

      # Launch interactive dashboard
      $ carriercapture viz --port 8050

      # Generate static plot
      $ carriercapture plot potential.json --show-wf --show

    For more information on each command, use:
      carriercapture COMMAND --help
    """
    # Store verbosity in context for subcommands
    ctx.ensure_object(dict)
    ctx.obj['verbose'] = verbose


# Import and register subcommands
from .commands import fit, solve, capture, viz, scan

cli.add_command(fit.fit_cmd)
cli.add_command(solve.solve_cmd)
cli.add_command(capture.capture_cmd)
cli.add_command(viz.viz_cmd)
cli.add_command(viz.plot_cmd)
cli.add_command(scan.scan_cmd)
cli.add_command(scan.scan_plot_cmd)


if __name__ == "__main__":
    cli()
