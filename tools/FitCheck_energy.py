#!/usr/bin/env python3
"""
Energy comparison scatter plot (Python port of FitCheck_energy.jl).
Uses matplotlib. Input: whitespace-separated text, 2 columns (observed, predicted)
or >=3 columns (use 2nd and 3rd). Energies in eV, displayed in meV.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

MARKER_SIZE = 5
MARKER_ALPHA = 0.8
AXIS_PADDING = 10.0  # meV

# Bright default colors for scatter (one per file): blue, red, then others
DEFAULT_COLORS = [
    "#17BECF",  # cyan
    "#EA580C",  # orange
    "#2563EB",  # blue (clear blue)
    "#DC2626",  # red (clear red)
    "#16A34A",  # green
    "#7C3AED",  # purple
]


def parse_index_list(spec: str) -> list[int]:
    """Parse a string like '1,3,5-10' into a list of integer indices."""
    indices = []
    for token in spec.split(","):
        t = token.strip()
        if not t:
            continue
        if "-" in t:
            parts = t.split("-", 1)
            if len(parts) != 2:
                raise ValueError(f"Invalid range token: {t}")
            try:
                start_idx = int(parts[0].strip())
                end_idx = int(parts[1].strip())
            except ValueError as e:
                raise ValueError(f"Invalid integer in range: {t}") from e
            if start_idx > end_idx:
                raise ValueError(f"Range start greater than end: {t}")
            indices.extend(range(start_idx, end_idx + 1))
        else:
            try:
                indices.append(int(t))
            except ValueError as e:
                raise ValueError(f"Invalid integer token: {t}") from e
    return indices


def calculate_statistics(y1: np.ndarray, y2: np.ndarray) -> dict[str, float]:
    """Compute summary statistics between observed and predicted values."""
    diff = y1 - y2
    rmse = np.sqrt(np.mean(diff**2))
    max_error = np.max(np.abs(diff))
    mean_error = np.mean(np.abs(diff))
    std_error = np.std(diff)
    ss_res = np.sum(diff**2)
    ss_tot = np.sum((y1 - np.mean(y1)) ** 2)
    r_squared = 1.0 - ss_res / ss_tot if ss_tot != 0 else float("nan")
    return {
        "RMSE": rmse,
        "Max Error": max_error,
        "Mean Error": mean_error,
        "Std Error": std_error,
        "R²": r_squared,
    }


def parse_file(path: str) -> tuple[np.ndarray, np.ndarray]:
    """
    Parse a text file; return (observed, predicted) arrays.
    Skip lines starting with # and empty lines.
    Either 2 columns (observed, predicted) or >=3 (use 2nd and 3rd).
    """
    observed, predicted = [], []
    with open(path) as f:
        for li, line in enumerate(f, 1):
            line_str = line.strip()
            if not line_str or line_str.startswith("#"):
                continue
            cols = line_str.split()
            if len(cols) == 2:
                o, p = float(cols[0]), float(cols[1])
            elif len(cols) >= 3:
                o, p = float(cols[1]), float(cols[2])
            else:
                raise ValueError(
                    f"Line {li} has {len(cols)} elements, expected 2 or >=3: {line!r}"
                )
            observed.append(o)
            predicted.append(p)
    return np.array(observed), np.array(predicted)


def plot_energy(
    files: list[str],
    *,
    output: str | None = None,
    lim: float | None = None,
    lim_min: float | None = None,
    lim_max: float | None = None,
    zero_min: bool = False,
    zero_at: float | None = None,
    colored_indices: list[int] | None = None,
    no_legend: bool = False,
    marker_size: float = MARKER_SIZE,
    marker_alpha: float = MARKER_ALPHA,
    tick_interval: float | None = None,
    per: int | None = None,
) -> None:
    unit = "meV"
    if per is not None and per != 0:
        unit = f"meV/{per}"

    observed_lists = []
    predicted_lists = []
    for path in files:
        obs, pre = parse_file(path)
        obs = obs * 1000.0  # eV -> meV
        pre = pre * 1000.0
        if per is not None and per != 0:
            obs = obs / per
            pre = pre / per
        observed_lists.append(obs)
        predicted_lists.append(pre)

    if zero_min and zero_at is not None:
        raise ValueError("Options --zero-min and --zero-at are mutually exclusive.")

    if zero_at is not None:
        for i in range(len(observed_lists)):
            observed_lists[i] = observed_lists[i] - zero_at
            predicted_lists[i] = predicted_lists[i] - zero_at
    elif zero_min:
        global_min = np.concatenate(observed_lists).min()
        for i in range(len(observed_lists)):
            observed_lists[i] = observed_lists[i] - global_min
            predicted_lists[i] = predicted_lists[i] - global_min
    else:
        for i in range(len(observed_lists)):
            center = (observed_lists[i].min() + observed_lists[i].max()) / 2
            observed_lists[i] = observed_lists[i] - center
            predicted_lists[i] = predicted_lists[i] - center

    for i, path in enumerate(files):
        stats = calculate_statistics(observed_lists[i], predicted_lists[i])
        print("-" * 30)
        print(f"File Index: {i}")
        print(f"File Name: {path}")
        print(f"RMSE: {stats['RMSE']:.4f} meV")
        print(f"R²: {stats['R²']:.4f}")
        print(f"Max Error: {stats['Max Error']:.4f} meV")
        print("-" * 30)

    # fig, ax = plt.subplots(figsize=(10, 10))
    fig, ax = plt.subplots()
    ax.set_aspect("equal")
    ax.set_title("Energy Comparison")
    ax.set_xlabel(f"DFT Energy ({unit})")
    ax.set_ylabel(f"SCE Energy ({unit})")
    ax.set_box_aspect(1)

    # Reference lines behind the scatter (zorder=0)
    big = 100_000.0
    ax.plot([-big, big], [-big, big], color="black", lw=1, zorder=0)
    ax.axvline(0, color="gray", ls="--", lw=1, zorder=0)
    ax.axhline(0, color="gray", ls="--", lw=1, zorder=0)

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
        xmin, xmax = local_min, local_max
        ymin, ymax = local_min, local_max
    else:
        all_obs = np.concatenate(observed_lists)
        all_pre = np.concatenate(predicted_lists)
        minv = min(all_obs.min(), all_pre.min())
        maxv = max(all_obs.max(), all_pre.max())
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

    for i, (obs, pre) in enumerate(zip(observed_lists, predicted_lists)):
        stats = calculate_statistics(obs, pre)
        series_label = f"{i}: {Path(files[i]).name} (RMSE: {stats['RMSE']:.2f} {unit})"

        if colored_indices is None:
            ax.scatter(
                obs,
                pre,
                s=marker_size,
                alpha=marker_alpha,
                label=series_label,
                zorder=5,
                color=DEFAULT_COLORS[i % len(DEFAULT_COLORS)],
                edgecolors="black",
                linewidths=0.5,
            )
        else:
            n = len(obs)
            mask_colored = np.zeros(n, dtype=bool)
            for idx in colored_indices:
                if 1 <= idx <= n:
                    mask_colored[idx - 1] = True
            mask_uncolored = ~mask_colored

            obs_u, pre_u = obs[mask_uncolored], pre[mask_uncolored]
            obs_c, pre_c = obs[mask_colored], pre[mask_colored]

            if len(obs_u) > 0:
                ax.scatter(
                    obs_u,
                    pre_u,
                    s=marker_size,
                    alpha=marker_alpha,
                    label=series_label,
                    zorder=5,
                    color=DEFAULT_COLORS[i % len(DEFAULT_COLORS)],
                    edgecolors="black",
                    linewidths=0.5,
                )
            if len(obs_c) > 0:
                ax.scatter(
                    obs_c,
                    pre_c,
                    s=marker_size,
                    alpha=marker_alpha,
                    color="orange",
                    label=f"{series_label} (colored)",
                    zorder=5,
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
            "Draw a meV scatter plot between observed and predicted energies (given in eV), "
            "after shifting each file by its own observed-value center (or optionally so that "
            "the global minimum or a specified energy becomes 0). Input: whitespace-separated text. "
            "Accepts 2 columns (observed predicted) or 3+ (index observed predicted). "
            "--lim fixes axes to [-lim, lim] in eV; --lim-min/--lim-max also in eV."
        ),
    )
    parser.add_argument("files", nargs="+", help="Input data files")
    parser.add_argument(
        "-o", "--output", default=None, help="Output filename (png, svg, pdf)"
    )
    parser.add_argument(
        "-l",
        "--lim",
        type=float,
        default=None,
        help="X/Y axis limits (eV). Fix to [-lim, lim]",
    )
    parser.add_argument(
        "--lim-min",
        type=float,
        default=None,
        help="Minimum X/Y axis limit (eV). With --lim-max.",
    )
    parser.add_argument(
        "--lim-max",
        type=float,
        default=None,
        help="Maximum X/Y axis limit (eV). With --lim-min.",
    )
    parser.add_argument(
        "-z",
        "--zero-min",
        action="store_true",
        help="Shift so global minimum observed energy is 0",
    )
    parser.add_argument(
        "--zero-at",
        type=float,
        default=None,
        help="Shift so this observed energy (eV) is 0 (mutually exclusive with --zero-min)",
    )
    parser.add_argument(
        "-c",
        "--colored",
        type=str,
        default=None,
        help="Indices to highlight (e.g. '5,7-10')",
    )
    parser.add_argument(
        "-m",
        "--marker-size",
        type=float,
        default=MARKER_SIZE,
        help="Marker size for scatter",
    )
    parser.add_argument(
        "-A",
        "--marker-alpha",
        type=float,
        default=MARKER_ALPHA,
        help="Marker transparency (0-1)",
    )
    parser.add_argument(
        "-L",
        "--no-legend",
        action="store_true",
        help="Disable legend",
    )
    parser.add_argument(
        "-T",
        "--tick-interval",
        type=float,
        default=None,
        metavar="meV",
        help="Major tick interval (meV). Same for X and Y. Default: auto.",
    )
    parser.add_argument(
        "--per",
        type=int,
        default=None,
        help="Divide all energies by this integer (e.g. number of atoms per formula unit) and show results per that unit.",
    )
    args = parser.parse_args()

    lim_mev = args.lim * 1000.0 if args.lim is not None else None
    lim_min_mev = args.lim_min * 1000.0 if args.lim_min is not None else None
    lim_max_mev = args.lim_max * 1000.0 if args.lim_max is not None else None
    zero_at_mev = args.zero_at * 1000.0 if args.zero_at is not None else None

    colored_indices = None
    if args.colored is not None:
        colored_indices = parse_index_list(args.colored)

    plot_energy(
        args.files,
        output=args.output,
        lim=lim_mev,
        lim_min=lim_min_mev,
        lim_max=lim_max_mev,
        zero_min=args.zero_min,
        zero_at=zero_at_mev,
        colored_indices=colored_indices,
        no_legend=args.no_legend,
        marker_size=args.marker_size,
        marker_alpha=args.marker_alpha,
        tick_interval=args.tick_interval,
        per=args.per,
    )


if __name__ == "__main__":
    main()
