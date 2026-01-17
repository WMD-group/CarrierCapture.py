"""
Visualization command: Launch interactive dashboard or generate static plots.

Provides options to:
- Launch Dash web dashboard for interactive visualization
- Generate static plots from data files
- Export publication-quality figures
"""

import click
from pathlib import Path


@click.command(name="viz")
@click.option(
    "--port",
    type=int,
    default=8050,
    help="Port for Dash server (default: 8050)"
)
@click.option(
    "--host",
    type=str,
    default="127.0.0.1",
    help="Host address (default: 127.0.0.1)"
)
@click.option(
    "--debug",
    is_flag=True,
    help="Run in debug mode with hot reloading"
)
@click.option(
    "--no-browser",
    is_flag=True,
    help="Don't automatically open browser"
)
@click.option(
    "--data",
    type=click.Path(exists=True, path_type=Path),
    help="Load data file on startup"
)
@click.pass_context
def viz_cmd(ctx, port, host, debug, no_browser, data):
    """
    Launch interactive visualization dashboard.

    Opens a web-based dashboard for exploring potential energy surfaces,
    fitting potentials, solving Schrödinger equation, and visualizing
    carrier capture rates.

    \\b
    Examples:
      # Launch dashboard on default port
      $ carriercapture viz

      # Launch on custom port
      $ carriercapture viz --port 8080

      # Launch with data file preloaded
      $ carriercapture viz --data potential.json

      # Run in debug mode
      $ carriercapture viz --debug
    """
    verbose = ctx.obj.get('verbose', 0)

    if verbose > 0:
        click.echo("Starting CarrierCapture visualization dashboard...")
        if data:
            click.echo(f"Preloading data from: {data}")

    try:
        from carriercapture.visualization.interactive import run_server

        # Open browser if not disabled
        if not no_browser:
            import webbrowser
            import threading

            def open_browser():
                import time
                time.sleep(1.5)  # Wait for server to start
                webbrowser.open(f"http://{host}:{port}")

            threading.Thread(target=open_browser, daemon=True).start()

        # Start server
        run_server(port=port, debug=debug, host=host)

    except ImportError as e:
        click.echo(
            f"Error: Dash dependencies not available. Install with: pip install dash\n"
            f"Details: {e}",
            err=True
        )
        ctx.exit(1)
    except Exception as e:
        click.echo(f"Error starting visualization dashboard: {e}", err=True)
        if debug:
            raise
        ctx.exit(1)


@click.command(name="plot")
@click.argument(
    "potential_file",
    type=click.Path(exists=True, path_type=Path),
)
@click.option(
    "--type",
    "plot_type",
    type=click.Choice(["potential", "spectrum", "both"]),
    default="potential",
    help="Type of plot to generate"
)
@click.option(
    "--show-wf",
    is_flag=True,
    help="Show wavefunctions on potential plot"
)
@click.option(
    "--max-wf",
    type=int,
    default=20,
    help="Maximum number of wavefunctions to plot"
)
@click.option(
    "-o", "--output",
    type=click.Path(path_type=Path),
    help="Output file (.html, .png, .pdf, .svg)"
)
@click.option(
    "--width",
    type=int,
    default=900,
    help="Figure width in pixels"
)
@click.option(
    "--height",
    type=int,
    default=600,
    help="Figure height in pixels"
)
@click.option(
    "--show",
    is_flag=True,
    help="Display plot in browser"
)
@click.pass_context
def plot_cmd(ctx, potential_file, plot_type, show_wf, max_wf, output, width, height, show):
    """
    Generate static plots from potential data.

    \\b
    Examples:
      # Plot potential energy surface
      $ carriercapture plot potential.json --show

      # Plot with wavefunctions
      $ carriercapture plot potential.json --show-wf -o figure.html

      # Plot eigenvalue spectrum
      $ carriercapture plot potential.json --type spectrum -o spectrum.png

      # Generate both plots
      $ carriercapture plot potential.json --type both --show-wf
    """
    verbose = ctx.obj.get('verbose', 0)

    if verbose > 0:
        click.echo(f"Loading potential from: {potential_file}")

    try:
        from carriercapture.io.readers import load_potential_from_file
        from carriercapture.core.potential import Potential
        from carriercapture.visualization.static import (
            plot_potential,
            plot_eigenvalue_spectrum,
        )

        # Load potential
        data = load_potential_from_file(potential_file)
        pot = Potential.from_dict(data)

        if verbose > 0:
            click.echo(f"Loaded potential: {pot.name or 'unnamed'}")
            if pot.eigenvalues is not None:
                click.echo(f"  {len(pot.eigenvalues)} eigenvalues")

        # Generate plots
        if plot_type in ["potential", "both"]:
            if verbose > 0:
                click.echo("Generating potential energy surface plot...")

            fig = plot_potential(
                pot,
                show_wavefunctions=show_wf,
                max_wf_to_plot=max_wf if show_wf else None,
            )
            fig.update_layout(width=width, height=height)

            if output:
                output_file = output if plot_type == "potential" else output.with_stem(output.stem + "_potential")
                fig.write_html(str(output_file))
                if verbose > 0:
                    click.echo(f"✓ Saved to: {output_file}")

            if show:
                fig.show()

        if plot_type in ["spectrum", "both"]:
            if pot.eigenvalues is None:
                click.echo("Error: Potential must be solved to plot spectrum", err=True)
                ctx.exit(1)

            if verbose > 0:
                click.echo("Generating eigenvalue spectrum plot...")

            fig = plot_eigenvalue_spectrum(pot)
            fig.update_layout(width=width, height=height)

            if output:
                output_file = output if plot_type == "spectrum" else output.with_stem(output.stem + "_spectrum")
                fig.write_html(str(output_file))
                if verbose > 0:
                    click.echo(f"✓ Saved to: {output_file}")

            if show:
                fig.show()

        if not output and not show:
            click.echo("Note: Use --output to save or --show to display plots")

    except Exception as e:
        click.echo(f"Error generating plot: {e}", err=True)
        if ctx.obj.get('verbose', 0) > 1:
            raise
        ctx.exit(1)
