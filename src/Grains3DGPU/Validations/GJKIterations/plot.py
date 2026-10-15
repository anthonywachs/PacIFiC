#!/usr/bin/env python3
"""Plot mean GJK iterations vs distance (log bins) from gjk_iterations.csv.

Usage: python3 plot.py --csv data/gjk_iterations.csv --plot-dir data
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

plt.rcParams.update({"text.usetex": True, "font.family": "Helvetica"})

# ---------------------------------------------------------------------------
# Style map: all subplots share the same colour.  #2282A4 is used for both
# the line/markers and (with transparency) the +/-1 std shaded band.
_LINE_COLOR  = "tab:blue"
_SHADE_COLOR = "tab:blue"

PAIR_STYLES = {
    "box-box":                   {"color": _LINE_COLOR, "shade": _SHADE_COLOR, "marker": "o", "label": r"Box--Box"},
    "cylinder-cylinder":         {"color": _LINE_COLOR, "shade": _SHADE_COLOR, "marker": "o", "label": r"Cylinder--Cylinder"},
    "superquadric-superquadric": {"color": _LINE_COLOR, "shade": _SHADE_COLOR, "marker": "o", "label": r"Superquadric--Superquadric"},
    "cylinder-box":              {"color": _LINE_COLOR, "shade": _SHADE_COLOR, "marker": "o", "label": r"Cylinder--Box"},
}

PAIR_ORDER = ["box-box", "cylinder-cylinder", "superquadric-superquadric", "cylinder-box"]


# ---------------------------------------------------------------------------
def bin_data(df_pair: pd.DataFrame, bins_per_decade: int = 5):
    """Bin by log10(distance); return geometric bin centres and mean/std/count of NbIter."""
    log_dist = np.log10(df_pair["ActualDistance"].values)
    n_iter   = df_pair["NbIter"].values.astype(float)

    log_min = np.floor(log_dist.min())
    log_max = np.ceil(log_dist.max())

    edges = np.arange(log_min, log_max + 1.0 / bins_per_decade,
                      1.0 / bins_per_decade)

    bin_centers, means, stds, counts = [], [], [], []
    for i in range(len(edges) - 1):
        mask = (log_dist >= edges[i]) & (log_dist < edges[i + 1])
        if mask.sum() == 0:
            continue
        vals = n_iter[mask]
        bin_centers.append(10.0 ** (0.5 * (edges[i] + edges[i + 1])))
        means.append(vals.mean())
        stds.append(vals.std())
        counts.append(int(mask.sum()))

    return (np.array(bin_centers),
            np.array(means),
            np.array(stds),
            np.array(counts))


# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Plot GJK iteration count vs distance from gjk_iterations.csv"
    )
    parser.add_argument(
        "--csv",
        dest="csv_path",
        type=str,
        required=True,
        help="Path to gjk_iterations.csv",
    )
    parser.add_argument(
        "--plot-dir",
        dest="plot_dir",
        type=str,
        required=True,
        help="Directory where the output PDF figure will be written",
    )
    args = parser.parse_args()

    if not os.path.exists(args.csv_path):
        msg = f"Missing results file '{args.csv_path}'"
        print(msg, file=sys.stderr)
        raise FileNotFoundError(msg)

    df = pd.read_csv(args.csv_path)

    required = {"ShapePair", "ActualDistance", "NbIter"}
    if not required.issubset(df.columns):
        raise ValueError(f"CSV is missing columns: {required - set(df.columns)}")

    # Keep only separated (non-overlapping), converged pairs.
    # NbIter == 1000 indicates the GJK iteration cap was hit (non-convergence);
    # those rows are excluded so they do not inflate the statistics.
    df = df[(df["ActualDistance"] > 0.0) & (df["NbIter"] < 1000)].copy()

    # -----------------------------------------------------------------------
    fig, axes = plt.subplots(1, 4, figsize=(16, 4),
                             sharex=True, sharey=True)

    for ax_idx, pair_name in enumerate(PAIR_ORDER):
        ax    = axes[ax_idx]
        style = PAIR_STYLES[pair_name]

        df_pair = df[df["ShapePair"] == pair_name]

        if df_pair.empty:
            ax.text(0.5, 0.5, "No data", ha="center", va="center",
                    transform=ax.transAxes)
            continue

        centers, means, stds, _ = bin_data(df_pair, bins_per_decade=5)

        sigma = 1.0
        means_s = gaussian_filter1d(means, sigma=sigma)
        stds_s  = gaussian_filter1d(stds,  sigma=sigma)

        ax.plot(
            centers, means_s,
            linestyle="-",
            linewidth=1.5,
            color=style["color"],
            marker=style["marker"],
            markersize=5,
            zorder=3,
            label=style["label"],
        )
        ax.fill_between(
            centers,
            means_s - stds_s,
            means_s + stds_s,
            color=style["shade"],
            alpha=0.3,
            linewidth=0,
            zorder=2,
        )

        ax.set_xscale("log")
        ax.set_box_aspect(1)
        ax.grid(color="lightgrey", linestyle="--", linewidth=0.5)

        ax.set_xlabel(r"Distance, $d$ [m]")
        if ax_idx == 0:
            ax.set_ylabel(r"GJK iterations $[-]$")

        ax.legend(loc="upper right", fontsize=8)

    fig.tight_layout()

    os.makedirs(args.plot_dir, exist_ok=True)
    out_path = os.path.join(args.plot_dir, "GJKIterationsVsDistance.pdf")
    plt.savefig(out_path, format="pdf", bbox_inches="tight")
    print(f"Figure saved to '{out_path}'")
    plt.close(fig)


if __name__ == "__main__":
    main()
