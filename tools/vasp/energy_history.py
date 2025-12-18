#!/usr/bin/env python3
"""
Extract energy at each SCF step from VASP OSZICAR file
and plot energy difference from final energy (in meV).
"""

import argparse
import re
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def parse_oszicar_energy(filepath):
    """
    Parse OSZICAR file and extract energy for each SCF step.

    Parameters:
    -----------
    filepath : str
        Path to OSZICAR file

    Returns:
    --------
    energies : list
        List of energies (in eV) for each SCF step
    steps : list
        List of SCF step numbers
    """
    with open(filepath, 'r') as f:
        lines = f.readlines()

    energies = []
    steps = []

    for i, line in enumerate(lines):
        # Detect SCF step line (3 uppercase letters followed by ':')
        # Format: "DAV:   1    -0.217621620828E+03   -0.21762E+03   -0.95322E-02 54425   0.254E+00"
        # Energy is the 3rd field (index 2) after splitting
        # Examples: DAV:, DMP:, RMM:, etc.
        step_match = re.match(r'^([A-Z]{3}):\s*(\d+)', line)
        if step_match:
            try:
                parts = line.split()
                # Energy is the 3rd field (index 2)
                if len(parts) >= 3:
                    step_num = int(step_match.group(2))
                    energy_str = parts[2]
                    energy = float(energy_str)
                    steps.append(step_num)
                    energies.append(energy)
            except (ValueError, IndexError):
                continue

    return energies, steps


def plot_energy_history(energies, steps, output_file=None, show_plot=True, ylim=None):
    """
    Plot energy history as difference from final energy (in meV).

    Parameters:
    -----------
    energies : list
        List of energies (in eV)
    steps : list
        List of SCF step numbers
    output_file : str, optional
        Output filename (if not specified, only display)
    show_plot : bool
        Whether to display the plot
    ylim : float, optional
        Symmetric y-axis limit (in meV). Sets range to [-ylim, ylim]
    """
    if not energies:
        print("No energy data found.")
        return

    # Convert to numpy array
    energies = np.array(energies)

    # Use final energy as reference (set to 0)
    final_energy = energies[-1]

    # Calculate energy difference from final energy (convert eV to meV)
    # 1 eV = 1000 meV
    energy_diff_mev = (energies - final_energy) * 1000.0

    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(steps, energy_diff_mev, marker='o', markersize=4, linewidth=1.5, color='blue')
    plt.axhline(y=0, color='r', linestyle='--', linewidth=1, alpha=0.5, label='Final Energy')
    plt.xlabel('SCF Step', fontsize=12)
    plt.ylabel('Energy Difference from Final Energy (meV)', fontsize=12)
    plt.title('Energy History (SCF Convergence)', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=10)

    # Set symmetric y-axis limits if specified
    if ylim is not None:
        plt.ylim(-ylim, ylim)

    plt.tight_layout()

    # Print statistics
    print(f"Number of SCF steps: {len(energies)}")
    print(f"Final energy: {final_energy:.6f} eV")
    print(f"Energy range: {np.min(energies):.6f} - {np.max(energies):.6f} eV")
    print(f"Energy difference range: {np.min(energy_diff_mev):.3f} - {np.max(energy_diff_mev):.3f} meV")

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {output_file}.")

    if show_plot:
        plt.show()
    else:
        plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Plot energy history from VASP OSZICAR file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python energy_history.py OSZICAR
  python energy_history.py OSZICAR --output energy_history.png
  python energy_history.py OSZICAR --no-show
  python energy_history.py OSZICAR --ylim 10
        """
    )

    parser.add_argument(
        'oszicar_file', type=str,
        help='Path to OSZICAR file'
    )
    parser.add_argument(
        '--output', '-o', type=str, default=None,
        help='Output filename (if not specified, only display)'
    )
    parser.add_argument(
        '--no-show', action='store_true',
        help='Do not display plot (useful when --output is specified)'
    )
    parser.add_argument(
        '--ylim', type=float, default=None,
        help='Symmetric y-axis limit (in meV). Sets range to [-ylim, ylim]'
    )

    args = parser.parse_args()

    # Check if file exists
    oszicar_path = Path(args.oszicar_file)
    if not oszicar_path.exists():
        print(f"Error: File not found: {args.oszicar_file}")
        return 1

    # Parse OSZICAR file
    print(f"Loading OSZICAR file: {args.oszicar_file}")
    energies, steps = parse_oszicar_energy(args.oszicar_file)

    if not energies:
        print("Error: Failed to extract energy data.")
        return 1

    # Plot
    show_plot = not args.no_show
    plot_energy_history(energies, steps, args.output, show_plot, args.ylim)

    return 0


if __name__ == '__main__':
    exit(main())
