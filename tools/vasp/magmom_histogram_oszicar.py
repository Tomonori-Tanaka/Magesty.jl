#!/usr/bin/env python3
"""
Extract magnetic moment magnitudes for the last SCF step from VASP OSZICAR file
and plot them as a histogram.
"""

import argparse
import re
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def parse_atom_spec(spec_str):
    """
    Parse atom specification string into a set of atom indices.

    Supports:
    - Range specification: "1-32" -> {1, 2, ..., 32}
    - Individual atoms: "1,2,5" -> {1, 2, 5}
    - Mixed: "1-5,10,15-20" -> {1,2,3,4,5,10,15,16,17,18,19,20}

    Parameters:
    -----------
    spec_str : str
        Atom specification string (e.g., "1-32" or "1,2,5" or "1-5,10,15-20")

    Returns:
    --------
    atoms : set
        Set of atom indices
    """
    atoms = set()
    if not spec_str:
        return atoms

    # Split by comma
    parts = spec_str.split(',')
    for part in parts:
        part = part.strip()
        if '-' in part:
            # Range specification
            try:
                start, end = part.split('-')
                start = int(start.strip())
                end = int(end.strip())
                atoms.update(range(start, end + 1))
            except ValueError:
                raise ValueError(f"Invalid range specification: {part}")
        else:
            # Individual atom
            try:
                atoms.add(int(part.strip()))
            except ValueError:
                raise ValueError(f"Invalid atom specification: {part}")

    return atoms


def parse_oszicar(filepath, magmom_type='MW_int', atom_filter=None):
    """
    Parse OSZICAR file and extract magnetic moments for each atom at each SCF step.

    Parameters:
    -----------
    filepath : str
        Path to OSZICAR file
    magmom_type : str
        Which magnetic moment type to extract: 'MW_int' or 'M_int'
    atom_filter : set, optional
        Set of atom indices to include (if None, include all atoms)

    Returns:
    --------
    data : dict
        Dictionary in the form {step: {ion_index: magnitude}}
    """
    with open(filepath, 'r') as f:
        lines = f.readlines()

    data = {}
    current_step = None
    in_magmom_section = False
    pending_data = {}  # Store data before DAV: line is encountered

    for i, line in enumerate(lines):
        # Detect start of SCF step (3 uppercase letters followed by ':')
        # Examples: DAV:, DMP:, RMM:, etc.
        match = re.match(r'^([A-Z]{3}):\s+(\d+)', line)
        if match:
            # When we encounter a DAV: line, assign pending data to this step
            # The data before DAV: N belongs to step N
            step_number = int(match.group(2))
            if pending_data:
                data[step_number] = pending_data.copy()
                pending_data.clear()
            else:
                # Initialize empty dict if no pending data
                data[step_number] = {}

            # Set the current step number for next iteration
            current_step = step_number
            in_magmom_section = False
            continue

        # Detect start of magnetic moment section
        if 'ion' in line and 'MW_int' in line and 'M_int' in line:
            in_magmom_section = True
            continue

        # Read magnetic moment data
        if in_magmom_section:
            # Extract ion number and numerical values
            # Format: " 1  0.342194521 -0.714762764 -0.200309003    0.492262481 -1.030728193 -0.287560733"
            parts = line.split()
            if len(parts) >= 7:
                try:
                    ion_index = int(parts[0])

                    # Filter atoms if specified
                    if atom_filter is not None and ion_index not in atom_filter:
                        continue

                    if magmom_type == 'MW_int':
                        # MW_int x, y, z components (indices 1, 2, 3)
                        mx = float(parts[1])
                        my = float(parts[2])
                        mz = float(parts[3])
                    elif magmom_type == 'M_int':
                        # M_int x, y, z components (indices 4, 5, 6)
                        mx = float(parts[4])
                        my = float(parts[5])
                        mz = float(parts[6])
                    else:
                        raise ValueError(f"Unknown magmom_type: {magmom_type}")

                    # Calculate magnitude
                    magnitude = np.sqrt(mx**2 + my**2 + mz**2)
                    # Store in pending_data (will be assigned to current_step when DAV: is encountered)
                    pending_data[ion_index] = magnitude

                except (ValueError, IndexError):
                    # Skip if numerical conversion fails
                    continue

        # Detect end of section or next SCF step
        # Check for 3 uppercase letters followed by ':' pattern
        if in_magmom_section and (re.match(r'^[A-Z]{3}:', line) or line.strip() == ''):
            in_magmom_section = False

    # Assign any remaining pending data to the last step
    if pending_data and current_step is not None:
        data[current_step] = pending_data.copy()

    return data


def plot_magmom_histogram(data, output_file=None, show_plot=True, bins='auto', density=False):
    """
    Plot histogram of magnetic moment magnitudes for the last SCF step.

    Parameters:
    -----------
    data : dict
        Data obtained from parse_oszicar()
    output_file : str, optional
        Output filename (if not specified, only display)
    show_plot : bool
        Whether to display the plot
    bins : int or str, optional
        Number of bins for histogram (default: 'auto')
    density : bool, optional
        If True, normalize the histogram (default: False)
    """
    if not data:
        print("No data found.")
        return

    # Get the last step with data
    steps = sorted(data.keys())
    if not steps:
        print("No SCF steps found.")
        return

    # Find the last step that has data
    last_step = None
    last_step_data = None
    for step in reversed(steps):
        if data[step]:
            last_step = step
            last_step_data = data[step]
            break

    if last_step is None or not last_step_data:
        print("No data found in any SCF step.")
        return

    # Extract magnitudes
    magnitudes = list(last_step_data.values())

    if not magnitudes:
        print("No magnetic moment data found.")
        return

    # Create histogram
    plt.figure(figsize=(10, 6))
    plt.hist(magnitudes, bins=bins, density=density, alpha=0.7, edgecolor='black')

    plt.xlabel('Magnetic Moment Magnitude (μB)', fontsize=12)
    ylabel = 'Density' if density else 'Frequency'
    plt.ylabel(ylabel, fontsize=12)
    plt.title(f'Magnetic Moment Histogram (SCF Step {last_step})', fontsize=14)
    plt.grid(True, alpha=0.3, axis='y')

    # Add statistics
    mean_val = np.mean(magnitudes)
    std_val = np.std(magnitudes)
    median_val = np.median(magnitudes)
    stats_text = f'Mean: {mean_val:.3f} μB\nStd: {std_val:.3f} μB\nMedian: {median_val:.3f} μB\nN: {len(magnitudes)}'
    plt.text(0.98, 0.98, stats_text, transform=plt.gca().transAxes,
             fontsize=10, verticalalignment='top', horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()

    # Save file first if specified
    if output_file:
        try:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Plot saved to {output_file}.")
            # Close figure after saving to ensure file is written
            if not show_plot:
                plt.close()
                return
        except Exception as e:
            print(f"Error saving plot to {output_file}: {e}")
            plt.close()
            return

    # Show plot if requested
    if show_plot:
        plt.show()
    else:
        # Close figure if not showing to free memory
        plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Plot histogram of magnetic moment magnitudes from the last SCF step in VASP OSZICAR file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python magmom_histogram_oszicar.py OSZICAR --type MW_int
  python magmom_histogram_oszicar.py OSZICAR --type M_int --output histogram.png
  python magmom_histogram_oszicar.py OSZICAR --atoms 1-32 --bins 50
  python magmom_histogram_oszicar.py OSZICAR --atoms 1,2,5 --density
        """
    )

    parser.add_argument(
        'oszicar_file', type=str,
        help='Path to OSZICAR file'
    )
    parser.add_argument(
        '--type', '-t', type=str, choices=['MW_int', 'M_int'],
        default='MW_int',
        help='Magnetic moment type to use (MW_int or M_int) [default: MW_int]'
    )
    parser.add_argument(
        '--atoms', '-a', type=str, default=None,
        help='Atom specification (e.g., "1-32", "1,2,5", or "1-5,10,15-20"). '
             'If not specified, all atoms are included.'
    )
    parser.add_argument(
        '--output', '-o', type=str, default=None,
        help='Output filename (if not specified, only display)'
    )
    parser.add_argument(
        '--bins', '-b', type=str, default=None,
        help='Number of bins for histogram (default: auto). Use "auto" for automatic binning.'
    )
    parser.add_argument(
        '--density', action='store_true',
        help='Normalize the histogram (show density instead of frequency)'
    )
    parser.add_argument(
        '--no-show', action='store_true',
        help='Do not display plot (useful when --output is specified)'
    )

    args = parser.parse_args()

    # Handle bins argument
    if args.bins is None or args.bins.lower() == 'auto':
        bins = 'auto'
    else:
        try:
            bins = int(args.bins)
        except ValueError:
            print(f"Error: Invalid bins value: {args.bins}. Must be an integer or 'auto'.")
            return 1

    # Check if file exists
    oszicar_path = Path(args.oszicar_file)
    if not oszicar_path.exists():
        print(f"Error: File not found: {args.oszicar_file}")
        return 1

    # Parse atom specification
    atom_filter = None
    if args.atoms:
        try:
            atom_filter = parse_atom_spec(args.atoms)
            print(f"Atom filter: {sorted(atom_filter)}")
        except ValueError as e:
            print(f"Error: Invalid atom specification: {e}")
            return 1

    # Parse OSZICAR file
    print(f"Loading OSZICAR file: {args.oszicar_file}")
    print(f"Magnetic moment type: {args.type}")

    data = parse_oszicar(args.oszicar_file, args.type, atom_filter)

    if not data:
        print("Error: Failed to extract data.")
        return 1

    steps = sorted(data.keys())
    print(f"Number of SCF steps extracted: {len(steps)}")
    if steps:
        # Find the last step that has data
        last_step = None
        for step in reversed(steps):
            if data[step]:
                last_step = step
                break

        if last_step is not None:
            print(f"Using last SCF step with data: {last_step}")
            n_atoms = len(data[last_step])
            print(f"Number of atoms in step {last_step}: {n_atoms}")
            if atom_filter:
                print(f"Filtered atoms: {sorted(atom_filter)}")
        else:
            print("Warning: No data found in any SCF step.")

    # Plot histogram
    show_plot = not args.no_show
    plot_magmom_histogram(data, args.output, show_plot, bins, args.density)

    return 0


if __name__ == '__main__':
    exit(main())
