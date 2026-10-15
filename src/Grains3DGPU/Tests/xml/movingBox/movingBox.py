#!/usr/bin/env python3
"""Trajectory and convergence plots for movingBox (constant-g z vs analytical)."""
import os
import sys
import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use("Agg")
plt.rcParams.update({"text.usetex": True, "font.family": "Helvetica"})

# ---------------------------------------------------------------------------
# Physical parameters (must match the YAML)
# ---------------------------------------------------------------------------
V0 = 5.0      # initial vertical velocity  [m/s]
G  = -9.81    # gravitational acceleration [m/s^2]
Z0 = 0.0      # initial z position         [m]
T_END = 1.0   # simulation end time        [s]

DT_VALUES  = [1e-1, 1e-2, 1e-3, 1e-4]
DT_LABELS  = ["dt1e-1", "dt1e-2", "dt1e-3", "dt1e-4"]
# Same colour; markers + offset slices avoid overlap on trajectory plots.
INTEGRATORS = {
    "FirstOrder":  (r"$1^{\mathrm{st}}$-order explicit",  "tab:blue", "o", slice(0,  None, 10)),
    "SecondOrder": (r"$2^{\mathrm{nd}}$-order leap-frog", "tab:blue", "X", slice(5,  None, 10)),
}
FINEST_DT_LABEL = "dt1e-4"


def analytical_z(t):
    return G / 2.0 * t ** 2 + V0 * t + Z0


def load_position(results_root: str, root_name: str) -> np.ndarray:
    """Load *_position_z.dat -> [time, z]."""
    fname = os.path.join(results_root, f"{root_name}_position_z.dat")
    if not os.path.exists(fname):
        msg = f"Missing results file '{os.path.basename(fname)}' in '{results_root}'"
        print(msg, file=sys.stderr)
        raise FileNotFoundError(msg)
    return np.loadtxt(fname)


def plot_trajectory(results_root: str, plots_root: str) -> None:
    """z(t) at finest dt for both integrators vs analytical."""
    fig, ax = plt.subplots(figsize=(5, 5))

    t_fine = np.linspace(0.0, T_END, 500)
    ax.plot(t_fine, analytical_z(t_fine),
                         linewidth=1, linestyle='--', c="k",
                         zorder=1, label="Analytical solution")

    for key, (label, color, marker, slc) in INTEGRATORS.items():
        root_name = f"movingBox_{FINEST_DT_LABEL}_{key}"
        data = load_position(results_root, root_name)
        t, z = data[:, 0], data[:, 1]
        ax.plot(t[slc], z[slc],
                linewidth=0, c=color, marker=marker, markersize=6,
                markerfacecolor=color, zorder=2, label=label)

    handles, labels = ax.get_legend_handles_labels()
    handles = handles[1:] + [handles[0]]
    labels  = labels[1:]  + [labels[0]]
    ax.legend(handles, labels)

    ax.set_xlabel(r"Time [s]")
    ax.set_ylabel(r"Height [m]")
    ax.set_xlim(0.0, T_END)
    ax.set_box_aspect(1)
    ax.grid(color='lightgrey', linestyle='--', linewidth=0.5)

    out = os.path.join(plots_root, "MovingBoxTrajectory.eps")
    os.makedirs(plots_root, exist_ok=True)
    plt.savefig(out, format='eps', bbox_inches='tight')
    plt.close(fig)
    print(f"Wrote: {out}")


def plot_convergence(results_root: str, plots_root: str) -> None:
    """Linf error along trajectory vs dt (log-log) with reference slopes."""
    fig, ax = plt.subplots(figsize=(5, 5))

    dt_labels = DT_LABELS
    dt_values = np.array(DT_VALUES)

    for key, (label, _, marker, _slc) in INTEGRATORS.items():
        errors = np.zeros(len(dt_labels))
        for i, dt_label in enumerate(dt_labels):
            root_name = f"movingBox_{dt_label}_{key}"
            data = load_position(results_root, root_name)
            t, z_num = data[:, 0], data[:, 1]
            z_anal = analytical_z(t)
            errors[i] = np.max(np.abs(z_num - z_anal))
        ax.loglog(dt_values, errors,
                  linewidth=0, c='tab:blue', marker=marker, markerfacecolor='tab:blue',
                  label=label)

    dt_ref = dt_values[0]
    err_ref_1 = np.abs(analytical_z(T_END) - load_position(results_root,
                       f"movingBox_dt1e-1_FirstOrder")[-1, 1])
    ax.loglog(dt_values, err_ref_1 * (dt_values / dt_ref) ** 1,
              '--', linewidth=1, c="k",
              label=r"Slope 1 ($\mathcal{O}(\Delta t)$)")

    err_ref_2 = np.abs(analytical_z(T_END) - load_position(results_root,
                       f"movingBox_dt1e-1_SecondOrder")[-1, 1])
    ax.loglog(dt_values, err_ref_2 * (dt_values / dt_ref) ** 2,
              '-.', linewidth=1, c="k",
              label=r"Slope 2 ($\mathcal{O}(\Delta t^2)$)")

    ax.set_xlabel(r"$\Delta t$, time step [s]")
    ax.set_ylabel(r"$L^\infty$ error [m]")
    ax.set_box_aspect(1)
    ax.legend(fontsize=8)
    plt.grid(True, which="both", color='lightgrey', ls="--", linewidth=0.5)

    out = os.path.join(plots_root, "MovingBoxConvergence.eps")
    os.makedirs(plots_root, exist_ok=True)
    plt.savefig(out, format='eps', bbox_inches='tight')
    plt.close(fig)
    print(f"Wrote: {out}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate trajectory and convergence plots for the movingBox test"
    )
    parser.add_argument('--result-dir', dest='results_root', type=str, required=True,
                        help='Directory containing movingBox_*.dat result files')
    parser.add_argument('--plot-dir', dest='plots_root', type=str, required=True,
                        help='Directory where plots will be written')
    args = parser.parse_args()

    plot_trajectory(args.results_root, args.plots_root)
    plot_convergence(args.results_root, args.plots_root)


if __name__ == '__main__':
    main()
