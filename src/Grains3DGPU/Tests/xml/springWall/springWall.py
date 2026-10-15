#!/usr/bin/env python3
"""Trajectory and convergence plots for the springWall test vs analytical z(t)."""
import os
import sys
import argparse
import numpy as np
from scipy.optimize import brentq
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use("Agg")
plt.rcParams.update({"text.usetex": True, "font.family": "Helvetica"})

# ---------------------------------------------------------------------------
# Physical parameters (must match the YAML)
# ---------------------------------------------------------------------------
R    = 0.1     # sphere radius              [m]
RHO  = 1000.0  # sphere density             [kg/m^3]
KN   = 5.0e4   # normal spring stiffness    [N/m]
V0   = -1.0    # initial velocity (toward wall) [m/s]
G    = 9.81    # gravitational acceleration [m/s^2]
T_END = 0.25   # simulation end time -- one bounce only, avoids phase-drift saturation [s]

_MASS  = RHO * (4.0 / 3.0) * np.pi * R ** 3   # ~ 4.1888 kg
OMEGA  = np.sqrt(KN / _MASS)                   # ~ 109.25 rad/s
_Z_EQ  = R - _MASS * G / KN                    # shifted equilibrium z during contact

def _z_contact(t):
    return _Z_EQ + (R - _Z_EQ) * np.cos(OMEGA * t) + (V0 / OMEGA) * np.sin(OMEGA * t)

def _zdot_contact(t):
    return -(R - _Z_EQ) * OMEGA * np.sin(OMEGA * t) + V0 * np.cos(OMEGA * t)

T_C    = brentq(lambda t: _z_contact(t) - R, 0.5 * np.pi / OMEGA, 1.9 * np.pi / OMEGA)
_V1    = _zdot_contact(T_C)                    # departure velocity ~ +|V0| (en=1)
T_FLT  = 2.0 * abs(_V1) / G                   # free-flight duration ~ 0.2039 s
T_PER  = T_C + T_FLT                           # full bounce period   ~ 0.234 s

DT_VALUES = [1e-3, 1e-4, 1e-5, 1e-6]
DT_LABELS = ["dt1e-3", "dt1e-4", "dt1e-5", "dt1e-6"]

# Same colour; markers distinguish integrators. Staggered slices (~50 markers at finest dt).
INTEGRATORS = {
    "FirstOrder":  (r"$1^{\mathrm{st}}$-order explicit",
                    "tab:blue", "o", slice(0,     None, 5000)),
    "SecondOrder": (r"$2^{\mathrm{nd}}$-order leap-frog",
                    "tab:blue", "X", slice(2500,  None, 5000)),
}
FINEST_DT_LABEL = "dt1e-6"


# ---------------------------------------------------------------------------
def analytical_z(t: np.ndarray) -> np.ndarray:
    """Analytical sphere-centre z for one bounce (contact + free flight)."""
    in_contact = t <= T_C
    z_contact = _z_contact(t)
    s = t - T_C
    z_free = R + _V1 * s - 0.5 * G * s ** 2
    return np.where(in_contact, z_contact, z_free)


# ---------------------------------------------------------------------------
def load_sphere_z(results_root: str, root_name: str) -> np.ndarray:
    """Load *_position_z.dat -> columns [time, z_obstacle, z_sphere]."""
    fname = os.path.join(results_root, f"{root_name}_position_z.dat")
    if not os.path.exists(fname):
        msg = (f"Missing results file '{os.path.basename(fname)}' "
               f"in '{results_root}'")
        print(msg, file=sys.stderr)
        raise FileNotFoundError(msg)
    return np.loadtxt(fname)


# ---------------------------------------------------------------------------
# Trajectory plot  (finest dt, both integrators)
# ---------------------------------------------------------------------------
def plot_trajectory(results_root: str, plots_root: str) -> None:
    """z(t) for the finest dt, both integrators vs analytical solution."""
    fig, ax = plt.subplots(figsize=(5, 5))

    t_fine = np.linspace(0.0, T_END, 1000)
    ax.plot(t_fine, analytical_z(t_fine),
                         linewidth=1, linestyle="--", c="k",
                         zorder=1, label="Analytical solution")

    for key, (label, color, marker, slc) in INTEGRATORS.items():
        root_name = f"springWall_{FINEST_DT_LABEL}_{key}"
        data = load_sphere_z(results_root, root_name)
        t = data[:, 0]
        z = data[:, 2]
        ax.plot(t[slc], z[slc],
                linewidth=0, c=color, marker=marker, markersize=6,
                markerfacecolor=color, zorder=2, label=label)

    handles, labels = ax.get_legend_handles_labels()
    handles = handles[1:] + [handles[0]]
    labels  = labels[1:]  + [labels[0]]
    ax.legend(handles, labels, fontsize=8)

    ax.set_xlabel(r"Time [s]")
    ax.set_ylabel(r"Mass centre $z$ [m]")
    ax.set_xlim(0.0, T_END)
    ax.set_box_aspect(1)
    ax.grid(color="lightgrey", linestyle="--", linewidth=0.5)

    os.makedirs(plots_root, exist_ok=True)
    out = os.path.join(plots_root, "SpringWallTrajectory.eps")
    plt.savefig(out, format="eps", bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote: {out}")


# ---------------------------------------------------------------------------
# Convergence plot  (Linf error vs dt)
# ---------------------------------------------------------------------------
def plot_convergence(results_root: str, plots_root: str) -> None:
    """Log-log Linf position error vs dt for both integrators."""
    fig, ax = plt.subplots(figsize=(5, 5))

    dt_arr = np.array(DT_VALUES)

    first_errors: dict = {}
    for key, (label, color, marker, _slc) in INTEGRATORS.items():
        errors = np.zeros(len(DT_LABELS))
        for i, dt_label in enumerate(DT_LABELS):
            root_name = f"springWall_{dt_label}_{key}"
            data = load_sphere_z(results_root, root_name)
            t     = data[:, 0]
            z_num = data[:, 2]
            # First bounce only (t <= T_PER) to avoid phase-drift across cycles.
            mask  = t <= T_PER
            z_ref = analytical_z(t[mask])
            errors[i] = np.max(np.abs(z_num[mask] - z_ref))
        first_errors[key] = errors[0]
        ax.loglog(dt_arr, errors,
                  linewidth=0, c=color, marker=marker, markersize=7,
                  markerfacecolor=color, zorder=2, label=label)

    dt_ref = dt_arr[0]
    err_fo = first_errors.get(
        "FirstOrder",
        first_errors[next(iter(first_errors))]
    )
    err_so = first_errors.get(
        "SecondOrder",
        first_errors[next(iter(first_errors))]
    )
    ax.loglog(dt_arr, err_fo * (dt_arr / dt_ref) ** 1,
              "--", linewidth=1, c="k", zorder=1,
              label=r"Slope 1 ($\mathcal{O}(\Delta t)$)")
    ax.loglog(dt_arr, err_so * (dt_arr / dt_ref) ** 2,
              "-.", linewidth=1, c="k", zorder=1,
              label=r"Slope 2 ($\mathcal{O}(\Delta t^2)$)")

    ax.set_xlabel(r"$\Delta t$, time step [s]")
    ax.set_ylabel(r"$L^\infty$ error [m]")
    ax.set_box_aspect(1)
    ax.legend(fontsize=8)
    ax.grid(True, which="both", color="lightgrey", ls="--", linewidth=0.5)

    os.makedirs(plots_root, exist_ok=True)
    out = os.path.join(plots_root, "SpringWallConvergence.eps")
    plt.savefig(out, format="eps", bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote: {out}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def main() -> None:
    parser = argparse.ArgumentParser(
        description="Trajectory and convergence plots for the springWall test"
    )
    parser.add_argument(
        "--result-dir",
        default="./results/springWall",
        help="Directory containing the *_position_z.dat files",
    )
    parser.add_argument(
        "--plot-dir",
        default="./plots/springWall",
        help="Directory for output .eps plot files",
    )
    args = parser.parse_args()

    print(f"Spring-wall parameters (with gravity):")
    print(f"  r      = {R} m")
    print(f"  rho    = {RHO} kg/m^3  =>  m = {_MASS:.5f} kg")
    print(f"  kn     = {KN:.1e} N/m  =>  omega = {OMEGA:.4f} rad/s")
    print(f"  g      = {G} m/s^2")
    print(f"  z_eq   = {_Z_EQ:.6f} m  (shifted equilibrium during contact)")
    print(f"  T_c    = {T_C*1e3:.3f} ms  (contact duration, numerical)")
    print(f"  V1     = {_V1:.6f} m/s  (departure speed)")
    print(f"  T_flt  = {T_FLT*1e3:.3f} ms  (free-flight duration)")
    print(f"  T_per  = {T_PER*1e3:.3f} ms  (bounce period)")
    print(f"  V0     = {V0} m/s")
    print()

    plot_trajectory(args.result_dir, args.plot_dir)
    plot_convergence(args.result_dir, args.plot_dir)


if __name__ == "__main__":
    main()
