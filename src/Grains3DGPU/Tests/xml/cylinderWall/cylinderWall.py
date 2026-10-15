#!/usr/bin/env python3
import os
import sys
import math
import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# Use non-interactive backend
matplotlib.use("Agg")
plt.rcParams.update({"text.usetex": True, "font.family": "Helvetica", "font.size": 16})


def readResults(thetaArray, results_root):
    simResults = np.zeros((np.size(thetaArray), 1))
    for t, theta in np.ndenumerate(thetaArray):
        theta_val = 90 - int(theta)
        # Expected filename may be under results_root or inside a subdirectory
        # depending on how the writer `Directory` was interpreted by the simulation.
        expected_name = f"cylinderWall{theta_val}_translational_velocity_z.dat"
        fname = os.path.join(results_root, expected_name)
        if not os.path.exists(fname):
            msg = f"Missing results file '{expected_name}' in directory '{results_root}'"
            print(msg, file=sys.stderr)
            raise FileNotFoundError(msg)
        # read file lines, split, take second column
        try:
            with open(fname, "r") as fh:
                lines = [ln.strip().split() for ln in fh.readlines() if ln.strip()]
        except Exception as e:
            print(f"Failed to read {fname}: {e}", file=sys.stderr)
            raise

        if not lines:
            print(f"Empty results file: {fname}", file=sys.stderr)
            raise ValueError(fname)

        try:
            velocity = np.array([float(row[2]) for row in lines])
        except Exception as e:
            print(f"Failed to parse velocities in {fname}: {e}", file=sys.stderr)
            raise

        # detect change points
        dif = velocity - np.roll(velocity, 1)
        idx = np.where(dif != 0)[0]
        if idx.size == 0:
            # fallback: take last value
            simResults[t] = -velocity[-1]
            continue

        # find the last index of the first contiguous block
        # remove the first index if it's 0 (because roll introduces a difference at index 0)
        if idx[0] == 0 and idx.size > 1:
            idx = idx[1:]

        first_max_idx = idx[0]
        for i in range(idx.size - 1):
            if idx[i] + 1 != idx[i + 1]:
                first_max_idx = idx[i]
                break
            else:
                first_max_idx = idx[i]

        simResults[t] = -velocity[first_max_idx + 1]

    return simResults


def analyticalSolution(theta):
    R = 4.e-3
    l = 6.e-3
    rho = 1163
    pre_vel = 1.0
    epsilon = 0.85
    mass = rho * math.pi * R ** 2 * l
    Iyy = 1.0 / 12.0 * mass * (3 * R ** 2 + l ** 2)
    r = math.sqrt(R ** 2 + 0.25 * l ** 2)
    alpha = np.arctan(l / R / 2) * np.ones(np.size(theta))
    alpha[0] = math.pi / 2
    alpha[-1] = 0
    post_omega = (mass * pre_vel * (1 + epsilon) * r * np.cos(alpha + theta)) / (
        Iyy + mass * r ** 2 * (np.cos(alpha + theta)) ** 2
    )
    post_vel = post_omega * r * np.cos(alpha + theta) - epsilon * pre_vel
    return post_vel


def main():
    parser = argparse.ArgumentParser(description="Generate cylinderWall plots from result files")
    parser.add_argument("--result-dir", dest="results_root", required=True,
                        help="Path to the per-test results directory (must contain results files)")
    parser.add_argument("--plot-dir", dest="plots_root", required=True,
                        help="Path where plots will be written")
    args = parser.parse_args()
    results_root = args.results_root
    plots_root = args.plots_root

    thetaArray = np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 90])
    try:
        simResults = readResults(thetaArray, results_root)
    except Exception as e:
        print(f"Error reading results: {e}", file=sys.stderr)
        return 2

    theta = np.linspace(0, math.pi / 2, num=90, endpoint=True)
    analResults = analyticalSolution(theta)

    fig = plt.figure(figsize=(5, 5))
    ax1 = fig.add_subplot(111)
    ax1.plot(theta / math.pi * 180, analResults, linewidth=1, linestyle="--", c="k", zorder=1, label="Analytical")
    ax1.plot(
        thetaArray,
        simResults.flatten(),
        linewidth=0,
        c="tab:blue",
        marker="o",
        markerfacecolor="tab:blue",
        zorder=2,
        label="Numerical",
    )
    ax1.set_xlabel(r"Impact angle, $\theta$ [deg]")
    ax1.set_ylabel(r"Dimensionless rebound speed, $V^+_z/V^-_z$ [-]")
    ax1.set_xlim(-5, 95)
    ax1.set_xticks(np.arange(0, 91, step=10))
    ax1.set_ylim(-0.9, 0.9)
    ax1.set_box_aspect(1)
    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles[1:] + [handles[0]], labels[1:] + [labels[0]])
    # Place grid lines beneath plot elements
    ax1.set_axisbelow(True)
    ax1.grid(color="lightgrey", linestyle="--", linewidth=0.5, zorder=0)

    out_dir = os.path.join(plots_root)
    os.makedirs(out_dir, exist_ok=True)
    out_eps = os.path.join(out_dir, "CylWallImpactTranslational.eps")
    out_png = os.path.join(out_dir, "CylWallImpactTranslational.png")

    fig.savefig(out_eps, format="eps", bbox_inches="tight")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
