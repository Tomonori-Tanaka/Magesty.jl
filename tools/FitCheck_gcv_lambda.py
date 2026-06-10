#!/usr/bin/env python3
"""
Ridge GCV penalty-sweep plot for Magesty.jl.

Plots the generalized cross-validation (GCV) score against the ridge penalty
lambda (log-log), marks the GCV minimizer, and optionally overlays the
effective degrees of freedom on a secondary axis.

Input: the whitespace-separated text written by `write_gcv_lambda`, with
columns `lambda  GCV  effective_dof` (lines starting with # are comments).
The GCV score is in the weighted-objective unit, not eV^2.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


def parse_file(path: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Parse columns lambda, gcv, dof, skipping # comments and blank lines."""
    lam, gcv, dof = [], [], []
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
            lam.append(float(cols[0]))
            gcv.append(float(cols[1]))
            dof.append(float(cols[2]) if len(cols) >= 3 else np.nan)
    return np.array(lam), np.array(gcv), np.array(dof)


def plot_gcv_lambda(
    files: list[str],
    *,
    output: str | None = None,
    show_dof: bool = False,
) -> None:
    fig, ax = plt.subplots()
    ax.set_title("Ridge GCV penalty sweep")
    ax.set_xlabel(r"penalty $\lambda$")
    ax.set_ylabel("GCV score (weighted-objective unit)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.tick_params(direction="in")

    ax_dof = ax.twinx() if show_dof else None
    if ax_dof is not None:
        ax_dof.set_ylabel("effective dof  tr(H)")

    for i, path in enumerate(files):
        lam, gcv, dof = parse_file(path)
        name = Path(path).name
        finite = np.isfinite(gcv)
        if finite.any():
            j = int(np.nanargmin(np.where(finite, gcv, np.inf)))
            best_lam, best_gcv = lam[j], gcv[j]
            print("-" * 30)
            print(f"File: {name}")
            print(f"lambda_best = {best_lam:.6g}  (GCV = {best_gcv:.6g})")
        else:
            best_lam = best_gcv = None
            print(f"File: {name}: no finite GCV values")

        line = ax.plot(lam, gcv, marker="o", ms=4, label=name, zorder=5)[0]
        if best_lam is not None:
            ax.scatter(
                [best_lam],
                [best_gcv],
                s=90,
                facecolors="none",
                edgecolors=line.get_color(),
                linewidths=1.8,
                zorder=6,
            )
        if ax_dof is not None:
            ax_dof.plot(lam, dof, ls="--", lw=1, color=line.get_color(), alpha=0.6)

    ax.legend(loc="best")
    plt.tight_layout()
    if output:
        fig.savefig(output, dpi=300, bbox_inches="tight")
        print(f"\nPlot saved as '{output}'")
    plt.show()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("files", nargs="+", help="gcv_lambda.txt file(s)")
    parser.add_argument("-o", "--output", default=None, help="Output image path")
    parser.add_argument(
        "-d",
        "--dof",
        action="store_true",
        help="Overlay effective dof on a secondary axis",
    )
    args = parser.parse_args()
    plot_gcv_lambda(args.files, output=args.output, show_dof=args.dof)


if __name__ == "__main__":
    main()
