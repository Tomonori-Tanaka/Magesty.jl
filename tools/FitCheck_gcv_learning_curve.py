#!/usr/bin/env python3
"""
Data-sufficiency GCV learning-curve plot for Magesty.jl.

Plots the generalized cross-validation (GCV) score against the training-set
size, with error bars from the spread over random subset draws. A curve that
flattens with size indicates that enough training data is present. With --r2 the
GCV-based predictive R^2 is plotted instead (a fixed 0..1 scale, easier to read).

Input: the whitespace-separated text written by `write_gcv_learning_curve`, with
columns `size  GCV_mean  GCV_std  GCV_R2_mean  GCV_R2_std` (lines starting with #
are comments). The GCV score is in the weighted-objective unit, not eV^2; GCV_R2
= 1 - GCV/MSY reads on a fixed scale (1 perfect, 0 matches the null model).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


def parse_file(
    path: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Parse size, gcv_mean/std, gcv_r2_mean/std, skipping # comments / blanks."""
    size, mean, std, r2_mean, r2_std = [], [], [], [], []
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
            r2_mean.append(float(cols[3]) if len(cols) >= 4 else np.nan)
            r2_std.append(float(cols[4]) if len(cols) >= 5 else 0.0)
    return (
        np.array(size),
        np.array(mean),
        np.array(std),
        np.array(r2_mean),
        np.array(r2_std),
    )


def plot_gcv_learning_curve(
    files: list[str],
    *,
    output: str | None = None,
    logy: bool = False,
    use_r2: bool = False,
) -> None:
    fig, ax = plt.subplots()
    ax.set_xlabel("training-set size (number of configurations)")
    if use_r2:
        ax.set_title("Data-sufficiency curve (predictive $R^2$)")
        ax.set_ylabel("GCV predictive $R^2$  (1 perfect, 0 null)")
    else:
        ax.set_title("Data-sufficiency GCV sweep")
        ax.set_ylabel("GCV score (weighted-objective unit)")
    if logy:
        ax.set_yscale("log")
    ax.tick_params(direction="in")

    for path in files:
        size, mean, std, r2_mean, r2_std = parse_file(path)
        name = Path(path).name
        y, yerr = (r2_mean, r2_std) if use_r2 else (mean, std)
        ax.errorbar(
            size,
            y,
            yerr=yerr,
            marker="o",
            ms=5,
            capsize=3,
            label=name,
            zorder=5,
        )
        print("-" * 30)
        print(f"File: {name}")
        for s, m, sd, rm, rsd in zip(size, mean, std, r2_mean, r2_std):
            print(
                f"  size={s:5d}  GCV={m:.6g} +/- {sd:.3g}"
                f"   R^2={rm:.6g} +/- {rsd:.3g}"
            )

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
    parser.add_argument(
        "-r",
        "--r2",
        action="store_true",
        help="Plot predictive R^2 (fixed scale) instead of the raw GCV score",
    )
    args = parser.parse_args()
    plot_gcv_learning_curve(
        args.files, output=args.output, logy=args.logy, use_r2=args.r2
    )


if __name__ == "__main__":
    main()
