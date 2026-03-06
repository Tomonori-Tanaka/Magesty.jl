#!/usr/bin/env python3
"""
Torque comparison scatter plot (Python port of FitCheck_torque.jl).

Input format (per line, whitespace-separated; comments starting with '#' are ignored):

    atom_index  element  DFT_tx  DFT_ty  DFT_tz  SCE_tx  SCE_ty  SCE_tz

Energies are given in eV and converted to meV for plotting, except for
`--type dir` where unit vectors (direction components, unitless) are used.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

# Allow importing shared helpers from tools/utils when run as a script
TOOLS_DIR = os.path.dirname(__file__)  # .../tools
if TOOLS_DIR not in sys.path:
    sys.path.append(TOOLS_DIR)

from utils.utils import parse_atom_indices

MARKER_SIZE = 5
MARKER_ALPHA = 0.8
AXIS_PADDING = 10.0  # meV

# Bright default colors for scatter (one per file), shared conceptually with FitCheck_energy.py
DEFAULT_COLORS = [
    "#17BECF",  # cyan
    "#EA580C",  # orange
    "#2563EB",  # blue
    "#DC2626",  # red
    "#16A34A",  # green
    "#7C3AED",  # purple
]


def calculate_statistics(y1: np.ndarray, y2: np.ndarray) -> dict[str, float]:
    """Compute summary statistics between observed and predicted values."""
    diff = y1 - y2
    rmse = np.sqrt(np.mean(diff ** 2))
    max_error = float(np.max(np.abs(diff)))
    mean_error = float(np.mean(np.abs(diff)))
    std_error = float(np.std(diff))
    ss_res = float(np.sum(diff ** 2))
    ss_tot = float(np.sum((y1 - np.mean(y1)) ** 2))
    r_squared = 1.0 - ss_res / ss_tot if ss_tot != 0 else float("nan")
    return {
        "RMSE": rmse,
        "Max Error": max_error,
        "Mean Error": mean_error,
        "Std Error": std_error,
        "R²": r_squared,
    }


def parse_torque_file(path: str) -> list[dict]:
    """
    Parse a torque file and return list of dicts with:
        atom_index, element, dft_torque (np.array(3)), sce_torque (np.array(3))
    """
    data: list[dict] = []
    with open(path) as f:
        for li, line in enumerate(f, 1):
            line_str = line.strip()
            if not line_str or line_str.startswith("#"):
                continue
            cols = line_str.split()
            if len(cols) < 8:
                raise ValueError(
                    f"Line {li} has {len(cols)} elements, expected >=8: {line!r}"
                )
            atom_index = int(cols[0])
            element = cols[1]
            dft_vals = [float(cols[2]), float(cols[3]), float(cols[4])]
            sce_vals = [float(cols[5]), float(cols[6]), float(cols[7])]
            data.append(
                {
                    "atom_index": atom_index,
                    "element": element,
                    "dft_torque": np.array(dft_vals, dtype=float),
                    "sce_torque": np.array(sce_vals, dtype=float),
                }
            )
    return data


def filter_data(
    data: list[dict],
    atom_indices: list[int] | None,
    elements: list[str] | None,
) -> list[dict]:
    """
    Filter torque data by atom indices or elements.
    `atom_indices` and `elements` are mutually exclusive.
    """
    if atom_indices is not None and elements is not None:
        raise ValueError("atom_indices and elements options are mutually exclusive")

    if atom_indices is not None:
        atom_set = set(atom_indices)
        return [x for x in data if x["atom_index"] in atom_set]
    if elements is not None:
        elem_set = set(elements)
        return [x for x in data if x["element"] in elem_set]
    return data


def create_plot(plot_type: str) -> tuple[plt.Figure, plt.Axes, str]:
    """Create a matplotlib figure/axes configured for torque comparison."""
    if plot_type == "all":
        title = "Torque Component Comparison"
        xlabel = "DFT Torque (meV)"
        ylabel = "SCE Torque (meV)"
        unit = "meV"
    elif plot_type == "norm":
        title = "Torque Magnitude Comparison"
        xlabel = "DFT Torque Magnitude (meV)"
        ylabel = "SCE Torque Magnitude (meV)"
        unit = "meV"
    elif plot_type == "dir":
        title = "Torque Direction Comparison"
        xlabel = "DFT Torque Direction Component"
        ylabel = "SCE Torque Direction Component"
        unit = "unitless"
    else:
        raise ValueError(f"Invalid plot_type: {plot_type}")

    fig, ax = plt.subplots()
    ax.set_aspect("equal")
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_box_aspect(1)

    # Reference lines behind the scatter
    big = 100_000.0
    ax.plot([-big, big], [-big, big], color="black", lw=1, label="y = x", zorder=0)
    ax.axvline(0, color="gray", ls="--", lw=1, zorder=0)
    ax.axhline(0, color="gray", ls="--", lw=1, zorder=0)

    return fig, ax, unit


def plot_torque(
    files: list[str],
    plot_type: str,
    *,
    output: str | None = None,
    lim: float | None = None,
    lim_min: float | None = None,
    lim_max: float | None = None,
    atom_indices: list[int] | None = None,
    elements: list[str] | None = None,
    no_legend: bool = False,
    marker_size: float = MARKER_SIZE,
    marker_alpha: float = MARKER_ALPHA,
    tick_interval: float | None = None,
) -> None:
    """
    Load torque data from files (convert eV -> meV for 'all' and 'norm'),
    shift each file by its own observed-value center, then draw scatter and y=x line.

    If `lim` is provided, X/Y axes are fixed to [-lim, lim] (meV).
    If `lim_min` and/or `lim_max` are provided, X/Y axes are fixed to [lim_min, lim_max] (meV).
    When only one of `lim_min` or `lim_max` is given, the other is set to the same absolute value
    with the opposite sign.
    """
    observed_lists: list[np.ndarray] = []
    predicted_lists: list[np.ndarray] = []

    for path in files:
        data = parse_torque_file(path)
        filtered = filter_data(data, atom_indices, elements)
        if not filtered:
            print(f"Warning: No data found in {path} after filtering")
            continue

        observed: list[float] = []
        predicted: list[float] = []

        if plot_type == "all":
            # Plot all x, y, z components (eV -> meV)
            for item in filtered:
                obs = item["dft_torque"] * 1000.0
                pre = item["sce_torque"] * 1000.0
                observed.extend(obs.tolist())
                predicted.extend(pre.tolist())
        elif plot_type == "norm":
            # Plot magnitude (eV -> meV)
            for item in filtered:
                obs_norm = float(np.linalg.norm(item["dft_torque"])) * 1000.0
                pre_norm = float(np.linalg.norm(item["sce_torque"])) * 1000.0
                observed.append(obs_norm)
                predicted.append(pre_norm)
        elif plot_type == "dir":
            # Plot direction components (unit vectors, unitless)
            for item in filtered:
                obs_vec = item["dft_torque"]
                pre_vec = item["sce_torque"]
                obs_norm = float(np.linalg.norm(obs_vec))
                pre_norm = float(np.linalg.norm(pre_vec))
                if obs_norm > 1e-10 and pre_norm > 1e-10:
                    obs_unit = obs_vec / obs_norm
                    pre_unit = pre_vec / pre_norm
                    observed.extend(obs_unit.tolist())
                    predicted.extend(pre_unit.tolist())
        else:
            raise ValueError(f"Invalid plot_type: {plot_type}")

        observed_lists.append(np.array(observed, dtype=float))
        predicted_lists.append(np.array(predicted, dtype=float))

    if not observed_lists:
        raise RuntimeError("No data found in any file after filtering")

    fig, ax, unit = create_plot(plot_type)

    # Optionally hide legend
    if no_legend:
        ax.legend([], [], frameon=False)

    # Print statistics per dataset
    for i, path in enumerate(files):
        if i >= len(observed_lists) or observed_lists[i].size == 0:
            continue
        stats = calculate_statistics(observed_lists[i], predicted_lists[i])
        print("-" * 30)
        print(f"File Index: {i}")
        print(f"File Name: {path}")
        print(f"RMSE: {stats['RMSE']:.4f} {unit}")
        print(f"Max Error: {stats['Max Error']:.4f} {unit}")
        print(f"Mean Error: {stats['Mean Error']:.4f} {unit}")
        print(f"Std Error: {stats['Std Error']:.4f} {unit}")
        print(f"R²: {stats['R²']:.4f}")
        print("-" * 30)

    # Determine axis limits (all modes use same logic; for 'dir' this is in unitless space)
    if lim is not None:
        xmin, xmax = -lim, lim
        ymin, ymax = -lim, lim
    elif lim_min is not None or lim_max is not None:
        local_min = lim_min
        local_max = lim_max
        if local_min is None and local_max is not None:
            local_min = -abs(local_max)
        elif local_min is not None and local_max is None:
            local_max = abs(local_min)
        xmin, xmax = local_min, local_max  # type: ignore[arg-type]
        ymin, ymax = local_min, local_max  # type: ignore[arg-type]
    else:
        all_observed = np.concatenate(observed_lists)
        all_predicted = np.concatenate(predicted_lists)
        minv = min(float(all_observed.min()), float(all_predicted.min()))
        maxv = max(float(all_observed.max()), float(all_predicted.max()))
        xmin = minv - AXIS_PADDING
        xmax = maxv + AXIS_PADDING
        ymin = minv - AXIS_PADDING
        ymax = maxv + AXIS_PADDING

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.grid(False)
    ax.tick_params(direction="in")
    if tick_interval is not None:
        ax.xaxis.set_major_locator(MultipleLocator(tick_interval))
        ax.yaxis.set_major_locator(MultipleLocator(tick_interval))

    # Scatter for each file dataset
    for i, (obs, pre) in enumerate(zip(observed_lists, predicted_lists)):
        if obs.size == 0:
            continue
        stats = calculate_statistics(obs, pre)
        series_label = (
            f"{i}: {Path(files[i]).name} (RMSE: {stats['RMSE']:.2f} {unit})"
        )
        ax.scatter(
            obs,
            pre,
            s=marker_size,
            alpha=marker_alpha,
            zorder=5,
            label=series_label,
            color=DEFAULT_COLORS[i % len(DEFAULT_COLORS)],
            edgecolors="black",
            linewidths=0.5,
        )

    if not no_legend:
        ax.legend(loc="upper left")

    plt.tight_layout()

    if output:
        fig.savefig(output, dpi=300, bbox_inches="tight")
        print(f"\nPlot saved as '{output}'")
    else:
        print("\nPlot not saved.")

    plt.show()


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Draw a scatter plot between observed and predicted magnetic torque (given in eV), "
            "after shifting each file by its own observed-value center. Input is whitespace-separated text. "
            "Format: atom_index element DFT_torque_x DFT_torque_y DFT_torque_z "
            "SCE_torque_x SCE_torque_y SCE_torque_z. "
            "If --lim is provided, axes are fixed to [-lim, lim] (meV). "
            "Alternatively, use --lim-min/--lim-max to fix axes to [lim-min, lim-max] (meV). "
            "Use --atom-indices or --elements to filter data (mutually exclusive)."
        )
    )
    parser.add_argument("files", nargs="+", help="Input data files")
    parser.add_argument(
        "-t",
        "--type",
        dest="plot_type",
        choices=["all", "norm", "dir"],
        default="all",
        help="Type of data to plot (all: components, norm: magnitude, dir: direction)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output filename (png, svg, pdf; format inferred from extension)",
    )
    parser.add_argument(
        "-l",
        "--lim",
        type=float,
        default=None,
        help="X/Y axis limits (meV). Fix to [-lim, lim]",
    )
    parser.add_argument(
        "--lim-min",
        type=float,
        default=None,
        help="Minimum X/Y axis limit (meV). Used together with --lim-max; "
        "if one side is omitted, it is set to the opposite sign with the same absolute value.",
    )
    parser.add_argument(
        "--lim-max",
        type=float,
        default=None,
        help="Maximum X/Y axis limit (meV). Used together with --lim-min; "
        "if one side is omitted, it is set to the opposite sign with the same absolute value.",
    )
    parser.add_argument(
        "-a",
        "--atom-indices",
        type=str,
        default=None,
        help="Atom indices to plot (comma-separated, e.g. '1,2,3' or ranges like '1-5,7,10-12')",
    )
    parser.add_argument(
        "-e",
        "--elements",
        type=str,
        default=None,
        help="Elements to plot (comma-separated, e.g. 'Fe,Co')",
    )
    parser.add_argument(
        "-m",
        "--marker-size",
        type=float,
        default=MARKER_SIZE,
        help=f"Marker size for scatter points (default: {MARKER_SIZE})",
    )
    parser.add_argument(
        "-A",
        "--marker-alpha",
        type=float,
        default=MARKER_ALPHA,
        help=f"Marker transparency (0-1) for scatter points (default: {MARKER_ALPHA})",
    )
    parser.add_argument(
        "-L",
        "--no-legend",
        action="store_true",
        help="Disable legend in the plot",
    )
    parser.add_argument(
        "-T",
        "--tick-interval",
        type=float,
        default=None,
        metavar="meV",
        help="Major tick interval (meV). Same for X and Y. Default: auto.",
    )

    args = parser.parse_args()

    atom_indices = None
    elements = None
    if args.atom_indices is not None:
        atom_indices = parse_atom_indices([args.atom_indices])
    if args.elements is not None:
        elements = [s.strip() for s in args.elements.split(",") if s.strip()]

    plot_torque(
        args.files,
        args.plot_type,
        output=args.output,
        lim=args.lim,
        lim_min=args.lim_min,
        lim_max=args.lim_max,
        atom_indices=atom_indices,
        elements=elements,
        no_legend=args.no_legend,
        marker_size=args.marker_size,
        marker_alpha=args.marker_alpha,
        tick_interval=args.tick_interval,
    )


if __name__ == "__main__":
    main()

