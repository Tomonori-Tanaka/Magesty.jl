#!/usr/bin/env python3
"""
Extract magnetic moment magnitudes for each atom at each SCF step from VASP OSZICAR file
and plot them as a graph.
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


def plot_magmom_history(data, output_file=None, show_plot=True):
    """
    Plot magnetic moment history.

    Parameters:
    -----------
    data : dict
        Data obtained from parse_oszicar()
    output_file : str, optional
        Output filename (if not specified, only display)
    show_plot : bool
        Whether to display the plot
    """
    if not data:
        print("No data found.")
        return

    # Sort step numbers
    steps = sorted(data.keys())

    # Get all ion indices
    all_ions = set()
    for step_data in data.values():
        all_ions.update(step_data.keys())
    all_ions = sorted(all_ions)

    # Convert data to array
    n_steps = len(steps)
    n_ions = len(all_ions)

    magmom_array = np.zeros((n_steps, n_ions))

    for i, step in enumerate(steps):
        for j, ion in enumerate(all_ions):
            if ion in data[step]:
                magmom_array[i, j] = data[step][ion]
            else:
                magmom_array[i, j] = np.nan

    # Plot
    plt.figure(figsize=(12, 8))

    # Plot for each ion
    for j, ion in enumerate(all_ions):
        plt.plot(steps, magmom_array[:, j], marker='o', markersize=3,
                 label=f'Ion {ion}', alpha=0.7, linewidth=1)

    plt.xlabel('SCF Step', fontsize=12)
    plt.ylabel('Magnetic Moment Magnitude (μB)', fontsize=12)
    plt.title('Magnetic Moment History', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8, ncol=2)
    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {output_file}.")

    if show_plot:
        plt.show()
    else:
        plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Plot magnetic moment history from VASP OSZICAR file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python magmom_history.py OSZICAR --type MW_int
  python magmom_history.py OSZICAR --type M_int --output magmom_history.png
  python magmom_history.py OSZICAR --atoms 1-32
  python magmom_history.py OSZICAR --atoms 1,2,5
  python magmom_history.py OSZICAR --atoms 1-5,10,15-20
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
        '--no-show', action='store_true',
        help='Do not display plot (useful when --output is specified)'
    )

    args = parser.parse_args()

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

    print(f"Number of SCF steps extracted: {len(data)}")
    if data:
        n_ions = len(data[list(data.keys())[0]])
        print(f"Number of atoms: {n_ions}")
        if atom_filter:
            print(f"Filtered atoms: {sorted(atom_filter)}")

    # Plot
    show_plot = not args.no_show
    plot_magmom_history(data, args.output, show_plot)

    return 0


if __name__ == '__main__':
    exit(main())
