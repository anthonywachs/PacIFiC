#!/usr/bin/env python3
"""Trajectory and convergence plots for collidingSpheresConv vs analytical two-sphere z(t)."""
import os
import sys
import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use("Agg")
plt.rcParams.update({"text.usetex": True, "font.family": "Helvetica", "font.size": 16})

# ---------------------------------------------------------------------------
# Physical parameters (must match the YAML)
# ---------------------------------------------------------------------------
R      = 0.1       # sphere radius               [m]
RHO    = 1000.0    # sphere density              [kg/m^3]
KN     = 5.0e4     # normal spring stiffness     [N/m]
V0     = -1.0      # initial velocity of sphere 1 [m/s]
V2I    = 0.0       # initial velocity of sphere 2 [m/s]
G      = -9.81     # gravitational acceleration  [m/s^2]
T_END  = 0.3       # simulation end time          [s]

Z1_INIT = 0.4      # initial z of sphere 1
Z2_INIT = 0.0      # initial z of sphere 2

M     = RHO * (4.0 / 3.0) * np.pi * R ** 3   # ~ 4.18879 kg
OMEGA = np.sqrt(2.0 * KN / M)                 # relative-motion angular freq ~ 154.53 rad/s
T_C   = np.pi / OMEGA                         # contact duration ~ 0.0203 s

T0    = (Z1_INIT - Z2_INIT - 2.0 * R) / (V2I - V0)

V1_AT_T0  = V0  + G * T0
V2_AT_T0  = V2I + G * T0
Z1_AT_T0  = Z1_INIT + V0  * T0 + 0.5 * G * T0 ** 2
Z2_AT_T0  = Z2_INIT + V2I * T0 + 0.5 * G * T0 ** 2

V_CM0 = (V1_AT_T0 + V2_AT_T0) / 2.0
Z_CM0 = (Z1_AT_T0 + Z2_AT_T0) / 2.0
V_REL = V1_AT_T0 - V2_AT_T0

Z_CM_END = Z_CM0 + V_CM0 * T_C + 0.5 * G * T_C ** 2
V_CM_END = V_CM0 + G * T_C
Z1_END   = Z_CM_END + R
Z2_END   = Z_CM_END - R
V1_POST  = V_CM_END + (-V_REL) / 2.0
V2_POST  = V_CM_END - (-V_REL) / 2.0

DT_VALUES  = [1e-3, 1e-4, 1e-5, 1e-6]
DT_LABELS  = ["dt1e-3", "dt1e-4", "dt1e-5", "dt1e-6"]
FINEST_DT_LABEL = "dt1e-6"

# Staggered slices so markers stay sparse at finest dt (~60 markers over 0.3 s).
INTEGRATORS = {
    "FirstOrder":  (r"$1^{\mathrm{st}}$-order explicit",
                    "tab:blue",   "o", slice(0,    None, 5000)),
    "SecondOrder": (r"$2^{\mathrm{nd}}$-order leap-frog",
                    "tab:blue", "X", slice(2500, None, 5000)),
}


# ---------------------------------------------------------------------------
def analytical_z(t: np.ndarray):
    """Piecewise z1, z2 with gravity: approach, harmonic contact, separation."""
    z1 = np.empty_like(t)
    z2 = np.empty_like(t)

    for i, ti in enumerate(t):
        if ti <= T0:
            # Phase 1: free approach until contact
            z1[i] = Z1_INIT + V0  * ti + 0.5 * G * ti ** 2
            z2[i] = Z2_INIT + V2I * ti + 0.5 * G * ti ** 2
        elif ti <= T0 + T_C:
            # Phase 2: harmonic relative motion + CM under g
            tau   = ti - T0
            x_rel = (V_REL / OMEGA) * np.sin(OMEGA * tau)
            z_cm  = Z_CM0 + V_CM0 * tau + 0.5 * G * tau ** 2
            z1[i] = z_cm + R + x_rel / 2.0
            z2[i] = z_cm - R - x_rel / 2.0
        else:
            # Phase 3: post-contact free flight
            s     = ti - T0 - T_C
            z1[i] = Z1_END + V1_POST * s + 0.5 * G * s ** 2
            z2[i] = Z2_END + V2_POST * s + 0.5 * G * s ** 2

    return z1, z2


# ---------------------------------------------------------------------------
def load_data(results_root: str, root_name: str) -> np.ndarray:
    """Load *_position_z.dat -> [time, z_sphere1, z_sphere2]."""
    fname = os.path.join(results_root, f"{root_name}_position_z.dat")
    if not os.path.exists(fname):
        msg = (f"Missing results file '{os.path.basename(fname)}' "
               f"in '{results_root}'")
        print(msg, file=sys.stderr)
        raise FileNotFoundError(msg)
    return np.loadtxt(fname)


# ---------------------------------------------------------------------------
def plot_trajectory(results_root: str, plots_root: str) -> None:
    """z(t) for both spheres (finest dt) vs analytical."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    t_fine = np.linspace(0.0, T_END, 2000)
    z1_anal, z2_anal = analytical_z(t_fine)

    sphere_info = [
        (1, "Sphere 1", z1_anal, r"Mass centre $z_1$ [m]"),
        (2, "Sphere 2", z2_anal, r"Mass centre $z_2$ [m]"),
    ]

    for ax, (col_idx, title, z_anal, ylabel) in zip(axes, sphere_info):
        ax.plot(t_fine, z_anal,
                linewidth=1, linestyle="--", c="k",
                zorder=1, label="Analytical solution")

        for key, (label, color, marker, slc) in INTEGRATORS.items():
            root_name = f"collidingSpheresConv_{FINEST_DT_LABEL}_{key}"
            data = load_data(results_root, root_name)
            t = data[:, 0]
            z = data[:, col_idx]
            ax.plot(t[slc], z[slc],
                    linewidth=0, c=color, marker=marker, markersize=5,
                    markerfacecolor=color, zorder=2, label=label)

        handles, labels = ax.get_legend_handles_labels()
        handles = handles[1:] + [handles[0]]
        labels  = labels[1:]  + [labels[0]]
        ax.legend(handles, labels, fontsize=8)

        ax.set_xlabel(r"Time [s]")
        ax.set_ylabel(ylabel)
        ax.set_xlim(0.0, T_END)
        ax.set_title(title, fontsize=10)
        ax.set_box_aspect(1)
        ax.set_axisbelow(True)
        ax.grid(color="lightgrey", linestyle="--", linewidth=0.5, zorder=0)

    os.makedirs(plots_root, exist_ok=True)
    fig.tight_layout()
    out = os.path.join(plots_root, "CollidingSpheresConvTrajectory.eps")
    plt.savefig(out, format="eps", bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote: {out}")


# ---------------------------------------------------------------------------
def plot_convergence(results_root: str, plots_root: str) -> None:
    """Log-log Linf position error vs dt for both integrators."""
    fig, ax = plt.subplots(figsize=(5, 5))
    dt_arr = np.array(DT_VALUES)
    first_errors: dict = {}

    for key, (label, color, marker, _) in INTEGRATORS.items():
        errors = np.zeros(len(DT_LABELS))
        for i, dt_label in enumerate(DT_LABELS):
            root_name = f"collidingSpheresConv_{dt_label}_{key}"
            data       = load_data(results_root, root_name)
            t, z1_num, z2_num = data[:, 0], data[:, 1], data[:, 2]
            z1_ref, z2_ref    = analytical_z(t)
            e1 = np.max(np.abs(z1_num - z1_ref))
            e2 = np.max(np.abs(z2_num - z2_ref))
            errors[i] = max(e1, e2)
        first_errors[key] = errors[0]
        ax.loglog(dt_arr, errors,
                  linewidth=0, c=color, marker=marker, markersize=7,
                  markerfacecolor=color, zorder=2, label=label)

    dt_ref = dt_arr[0]
    ax.loglog(dt_arr, first_errors["FirstOrder"]  * (dt_arr / dt_ref) ** 1,
              "--", linewidth=1, c="k", zorder=1,
              label=r"Slope 1 ($\mathcal{O}(\Delta t)$)")
    ax.loglog(dt_arr, first_errors["SecondOrder"] * (dt_arr / dt_ref) ** 2,
              "-.", linewidth=1, c="k", zorder=1,
              label=r"Slope 2 ($\mathcal{O}(\Delta t^2)$)")

    ax.set_xlabel(r"Time step, $\Delta t$ [s]")
    ax.set_ylabel(r"$L^\infty$ position error [m]")
    ax.set_box_aspect(1)
    ax.legend(fontsize=8)
    ax.set_axisbelow(True)
    ax.grid(True, which="both", color="lightgrey", ls="--", linewidth=0.5, zorder=0)

    os.makedirs(plots_root, exist_ok=True)
    out = os.path.join(plots_root, "CollidingSpheresConvConvergence.eps")
    plt.savefig(out, format="eps", bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote: {out}")


# ---------------------------------------------------------------------------
def main() -> None:
    parser = argparse.ArgumentParser(
        description="Trajectory and convergence plots for the collidingSpheresConv test"
    )
    parser.add_argument(
        "--result-dir",
        default="./results/collidingSpheresConv",
        help="Directory containing the *_position_z.dat files",
    )
    parser.add_argument(
        "--plot-dir",
        default="./plots/collidingSpheresConv",
        help="Directory for output .eps plot files",
    )
    args = parser.parse_args()

    print("Colliding Spheres Convergence -- physical parameters:")
    print(f"  R         = {R} m")
    print(f"  rho       = {RHO} kg/m^3  =>  m = {M:.5f} kg")
    print(f"  kn        = {KN:.1e} N/m  =>  omega = {OMEGA:.4f} rad/s")
    print(f"  g         = {G} m/s^2")
    print(f"  V0        = {V0} m/s  (initial velocity of sphere 1)")
    print(f"  t0        = {T0:.4f} s   (free-approach duration until contact)")
    print(f"  v1(t0)    = {V1_AT_T0:.4f} m/s  (sphere 1 velocity at contact)")
    print(f"  v2(t0)    = {V2_AT_T0:.4f} m/s  (sphere 2 velocity at contact)")
    print(f"  Tc        = {T_C * 1e3:.3f} ms  (contact duration)")
    print(f"  z1_end    = {Z1_END:.6f} m   z2_end    = {Z2_END:.6f} m")
    print(f"  v1_post   = {V1_POST:.4f} m/s  v2_post   = {V2_POST:.4f} m/s")
    print()

    plot_trajectory(args.result_dir, args.plot_dir)
    plot_convergence(args.result_dir, args.plot_dir)


if __name__ == "__main__":
    main()
