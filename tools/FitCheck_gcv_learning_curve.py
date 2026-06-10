#!/usr/bin/env python3
"""
Data-sufficiency GCV learning-curve plot for Magesty.jl.

Plots the generalized cross-validation (GCV) score against the training-set
size, with error bars from the spread over random subset draws. A curve that
flattens with size indicates that enough training data is present.

Input: the whitespace-separated text written by `write_gcv_learning_curve`, with
columns `size  GCV_mean  GCV_std` (lines starting with # are comments).
The GCV score is in the weighted-objective unit, not eV^2.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


def parse_file(path: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Parse columns size, gcv_mean, gcv_std, skipping # comments / blanks."""
    size, mean, std = [], [], []
    with open(path) as f:
        for li, line in enumerate(f, 1):
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            cols = s.split()
            if len(cols) < 2:
                raise ValueError(
                    f"Line {li} has {len(cols)} columns, expected >= 2: {line!r}"
                )
            size.append(int(float(cols[0])))
            mean.append(float(cols[1]))
            std.append(float(cols[2]) if len(cols) >= 3 else 0.0)
    return np.array(size), np.array(mean), np.array(std)


def plot_gcv_learning_curve(
    files: list[str],
    *,
    output: str | None = None,
    logy: bool = False,
) -> None:
    fig, ax = plt.subplots()
    ax.set_title("Data-sufficiency GCV sweep")
    ax.set_xlabel("training-set size (number of configurations)")
    ax.set_ylabel("GCV score (weighted-objective unit)")
    if logy:
        ax.set_yscale("log")
    ax.tick_params(direction="in")

    for path in files:
        size, mean, std = parse_file(path)
        name = Path(path).name
        ax.errorbar(
            size,
            mean,
            yerr=std,
            marker="o",
            ms=5,
            capsize=3,
            label=name,
            zorder=5,
        )
        print("-" * 30)
        print(f"File: {name}")
        for s, m, sd in zip(size, mean, std):
            print(f"  size={s:5d}  GCV={m:.6g} +/- {sd:.3g}")

    ax.legend(loc="best")
    plt.tight_layout()
    if output:
        fig.savefig(output, dpi=300, bbox_inches="tight")
        print(f"\nPlot saved as '{output}'")
    plt.show()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("files", nargs="+", help="gcv_learning_curve.txt file(s)")
    parser.add_argument("-o", "--output", default=None, help="Output image path")
    parser.add_argument(
        "--logy", action="store_true", help="Use a logarithmic y axis"
    )
    args = parser.parse_args()
    plot_gcv_learning_curve(args.files, output=args.output, logy=args.logy)


if __name__ == "__main__":
    main()
