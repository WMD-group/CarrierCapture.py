"""
Solve command: Solve Schrödinger equation for phonon states.

Takes a fitted potential and computes eigenvalues and eigenvectors.
"""

import click
from pathlib import Path
import numpy as np

from carriercapture.core import Potential
from carriercapture.io.readers import load_potential_from_file
from carriercapture.io.writers import save_potential


@click.command(name="solve")
@click.argument(
    "potential_file",
    type=click.Path(exists=True, path_type=Path),
)
@click.option(
    "-n", "--nev",
    type=int,
    default=60,
    help="Number of eigenvalues to compute"
)
@click.option(
    "-o", "--output",
    type=click.Path(path_type=Path),
    help="Output file (default: add '_solved' suffix)"
)
@click.option(
    "--Q-range",
    nargs=2,
    type=float,
    help="Q grid range: Q_min Q_max (amu^0.5·Å)"
)
@click.option(
    "--npoints",
    type=int,
    default=3001,
    help="Number of grid points for solving"
)
@click.option(
    "--plot",
    is_flag=True,
    help="Plot eigenvalue spectrum"
)
@click.option(
    "--plot-output",
    type=click.Path(path_type=Path),
    help="Save plot to file"
)
@click.option(
    "--plot-wavefunctions",
    type=int,
    help="Number of wavefunctions to plot"
)
@click.pass_context
def solve_cmd(ctx, potential_file, nev, output, q_range, npoints, plot, plot_output, plot_wavefunctions):
    """
    Solve Schrödinger equation for phonon states.

    \b
    POTENTIAL_FILE: Path to fitted potential (JSON, YAML, NPZ)

    \b
    Examples:
      # Solve with default parameters (60 states)
      $ carriercapture solve potential.json

      # Solve for more states with custom grid
      $ carriercapture solve potential.json -n 180 --npoints 5000

      # Solve and plot eigenvalues
      $ carriercapture solve potential.json --plot

      # Solve and plot wavefunctions
      $ carriercapture solve potential.json --plot-wavefunctions 10
    """
    verbose = ctx.obj.get('verbose', 0)

    # Load potential
    if verbose > 0:
        click.echo(f"Loading potential from: {potential_file}")

    try:
        data = load_potential_from_file(potential_file)
        pot = Potential.from_dict(data)
    except Exception as e:
        click.echo(f"Error loading potential: {e}", err=True)
        ctx.exit(1)

    # Check if potential is fitted
    if pot.fit_func is None:
        click.echo("Error: Potential must be fitted before solving", err=True)
        ctx.exit(1)

    if verbose > 0:
        click.echo(f"Loaded potential: {pot.name or 'unnamed'}")
        if pot.fit_type:
            click.echo(f"Fit type: {pot.fit_type}")

    # Set Q range if specified
    if q_range:
        pot.Q = np.linspace(q_range[0], q_range[1], npoints)
        pot.E = pot.fit_func(pot.Q)
        if verbose > 0:
            click.echo(f"Using Q range: [{q_range[0]:.3f}, {q_range[1]:.3f}]")

    # Solve
    if verbose > 0:
        click.echo(f"\nSolving Schrödinger equation for {nev} states...")
        click.echo(f"Grid: {len(pot.Q)} points from {pot.Q.min():.3f} to {pot.Q.max():.3f}")

    try:
        pot.solve(nev=nev)
    except Exception as e:
        click.echo(f"Error during solving: {e}", err=True)
        ctx.exit(1)

    if verbose > 0:
        click.echo("✓ Solving completed successfully")
        click.echo(f"\nEigenvalue spectrum:")
        click.echo(f"  E₀ = {pot.eigenvalues[0]:.6f} eV (ground state)")
        click.echo(f"  E₁ = {pot.eigenvalues[1]:.6f} eV")
        if len(pot.eigenvalues) > 2:
            click.echo(f"  E₂ = {pot.eigenvalues[2]:.6f} eV")
        click.echo(f"  ...")
        click.echo(f"  E_{nev-1} = {pot.eigenvalues[-1]:.6f} eV (highest)")
        spacing = np.diff(pot.eigenvalues)
        click.echo(f"\nAverage spacing: {spacing.mean():.6f} eV")
        click.echo(f"Min/Max spacing: {spacing.min():.6f} / {spacing.max():.6f} eV")

    # Save output
    if output is None:
        # Auto-generate output filename
        output = potential_file.with_stem(potential_file.stem + "_solved")

    if verbose > 0:
        click.echo(f"\nSaving to: {output}")

    try:
        save_potential(pot, output)
        if verbose > 0:
            click.echo("✓ Saved successfully")
    except Exception as e:
        click.echo(f"Error saving output: {e}", err=True)
        ctx.exit(1)

    # Plot if requested
    if plot or plot_wavefunctions:
        try:
            import matplotlib.pyplot as plt

            if plot_wavefunctions:
                # Plot potential + wavefunctions
                fig, ax = plt.subplots(figsize=(12, 8))

                # Plot potential
                ax.plot(pot.Q, pot.E, 'k-', linewidth=2, label="Potential", alpha=0.5)

                # Plot selected wavefunctions (scaled and shifted)
                n_wf = min(plot_wavefunctions, len(pot.eigenvalues))
                colors = plt.cm.viridis(np.linspace(0, 1, n_wf))

                for i in range(n_wf):
                    # Scale wavefunction for visibility
                    psi = pot.eigenvectors[i, :]
                    scale = 0.2 * (pot.E.max() - pot.E.min())
                    psi_scaled = psi * scale + pot.eigenvalues[i]

                    ax.plot(pot.Q, psi_scaled, color=colors[i], linewidth=1.5,
                            label=f"ψ_{i} (E={pot.eigenvalues[i]:.3f} eV)")

                    # Draw eigenvalue line
                    ax.axhline(pot.eigenvalues[i], color=colors[i], linestyle='--',
                              alpha=0.3, linewidth=0.8)

                ax.set_xlabel("Q (amu$^{0.5}$·Å)", fontsize=12)
                ax.set_ylabel("E (eV)", fontsize=12)
                ax.set_title(f"Potential with {n_wf} Wavefunctions", fontsize=14)
                ax.legend(loc='best', fontsize=8, ncol=2)
                ax.grid(True, alpha=0.3)

            else:
                # Just plot eigenvalue spectrum
                fig, ax = plt.subplots(figsize=(10, 6))

                # Plot potential
                ax.plot(pot.Q, pot.E, 'k-', linewidth=2, label="Potential")

                # Plot eigenvalues as horizontal lines
                for i, E_i in enumerate(pot.eigenvalues):
                    color = 'blue' if i < 10 else 'gray'
                    alpha = 0.8 if i < 10 else 0.3
                    ax.axhline(E_i, color=color, linestyle='--', alpha=alpha, linewidth=1)

                # Label first few levels
                for i in range(min(5, len(pot.eigenvalues))):
                    ax.text(pot.Q.max() * 0.95, pot.eigenvalues[i],
                           f' n={i}', fontsize=9, va='center')

                ax.set_xlabel("Q (amu$^{0.5}$·Å)", fontsize=12)
                ax.set_ylabel("E (eV)", fontsize=12)
                ax.set_title(f"Eigenvalue Spectrum ({nev} states)", fontsize=14)
                ax.legend()
                ax.grid(True, alpha=0.3)

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
