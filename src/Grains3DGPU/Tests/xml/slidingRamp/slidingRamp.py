#!/usr/bin/env python3
"""Post-processing: slidingRamp -- Hooke vs HookeMemory x(t) vs tilted-gravity analytics."""
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
R     = 0.1       # sphere radius              [m]
RHO   = 1000.0    # sphere density             [kg/m^3]
KN    = 5.0e4     # normal spring stiffness    [N/m]
KT    = 2.0e4     # tangential spring (HookeMemory only) [N/m]
EN    = 0.9       # normal restitution
MUC   = 0.8       # Coulomb friction coefficient
THETA = 20.0      # inclination angle          [ deg]
G     = 9.81      # gravitational acceleration [m/s^2]
T_END = 0.3       # simulation end time        [s]

M       = RHO * (4.0 / 3.0) * np.pi * R ** 3          # ~ 4.189 kg
THETA_R = np.radians(THETA)
G_TAN   = G * np.sin(THETA_R)                          # ~ 3.3552 m/s^2 (driving)
G_NOR   = G * np.cos(THETA_R)                          # ~ 9.2183 m/s^2 (normal)
F_N     = M * G_NOR                                    # ~ 38.60 N
OMEGA_T = np.sqrt(KT / M)                              # ~ 69.1 rad/s
X_EQ    = M * G_TAN / KT                               # ~ 7.03e-4 m (spring equil.)
T_OSC   = 2.0 * np.pi / OMEGA_T                       # ~ 0.0909 s (oscillation period)

_ft_peak = KT * 2.0 * X_EQ  # = 2*m*g*sin(theta)
_coulomb  = MUC * F_N
assert _ft_peak < _coulomb, (
    f"Peak spring force {_ft_peak:.2f} N exceeds Coulomb limit {_coulomb:.2f} N -- "
    "choose larger muc or smaller theta"
)


# ---------------------------------------------------------------------------
def analytical_hooke(t: np.ndarray):
    """Hooke (no memory, etat=0): pure free sliding under tangential gravity."""
    x = 0.5 * G_TAN * t ** 2
    z = np.full_like(t, R - M * G_NOR / KN)
    return x, z


def analytical_hookememory(t: np.ndarray):
    """HookeMemory (kt, etat=0): undamped harmonic oscillator around x_eq."""
    x = X_EQ * (1.0 - np.cos(OMEGA_T * t))
    z = np.full_like(t, R - M * G_NOR / KN)
    return x, z


# ---------------------------------------------------------------------------
def load_position(results_root: str, root_name: str, axis: str) -> np.ndarray:
    """Load *_position_{axis}.dat; columns: [time, z_obstacle, z_sphere]."""
    fname = os.path.join(results_root, f"{root_name}_position_{axis}.dat")
    if not os.path.exists(fname):
        msg = f"Missing '{os.path.basename(fname)}' in '{results_root}'"
        print(msg, file=sys.stderr)
        raise FileNotFoundError(msg)
    return np.loadtxt(fname)


# ---------------------------------------------------------------------------
def plot_comparison(results_root: str, plots_root: str) -> None:
    """Two figures: Hooke sliding and HookeMemory oscillation in x (cm)."""
    t_fine = np.linspace(0.0, T_END, 2000)
    x_hooke_anal, _ = analytical_hooke(t_fine)
    x_mem_anal, _ = analytical_hookememory(t_fine)

    marker_dt = 0.01

    os.makedirs(plots_root, exist_ok=True)

    fig1, ax1 = plt.subplots(figsize=(5, 5))
    ax1.plot(t_fine, x_hooke_anal * 100,
             linewidth=1, linestyle="--", c="k", zorder=1, label="Analytical")
    data = load_position(results_root, "slidingRamp_Hooke", "x")
    t, x_num = data[:, 0], data[:, 2]
    dt_data = t[1] - t[0]
    stride = max(1, int(round(marker_dt / dt_data)))
    ax1.plot(t[::stride], (x_num * 100)[::stride],
             linewidth=0, c="tab:blue", marker="o", markersize=5,
             markerfacecolor="tab:blue", zorder=2, label=r"Simulation ($k_t = 0$)")
    ax1.set_xlabel(r"Time, $t$ [s]")
    ax1.set_ylabel(r"Position, $x$ [cm]")
    ax1.set_xlim(0.0, T_END)
    ax1.set_box_aspect(1)
    handles, labels = ax1.get_legend_handles_labels()
    try:
        idx_anal = labels.index("Analytical")
        idx_sim = next(i for i, l in enumerate(labels) if "simulation" in l.lower())
        ax1.legend([handles[idx_sim], handles[idx_anal]], [labels[idx_sim], labels[idx_anal]], loc="upper left", fontsize=8)
    except Exception:
        ax1.legend(loc="upper left", fontsize=8)
    ax1.grid(color="lightgrey", linestyle="--", linewidth=0.5)
    fig1.tight_layout()
    out1 = os.path.join(plots_root, "SlidingRamp_HookePosition.eps")
    plt.savefig(out1, format="eps", bbox_inches="tight")
    plt.close(fig1)
    print(f"Wrote: {out1}")

    fig2, ax2 = plt.subplots(figsize=(5, 5))
    ax2.plot(t_fine, x_mem_anal * 100,
             linewidth=1, linestyle="--", c="k", zorder=1, label="Analytical")
    data = load_position(results_root, "slidingRamp_HookeMemory", "x")
    t, x_num = data[:, 0], data[:, 2]
    dt_data = t[1] - t[0]
    stride = max(1, int(round(marker_dt / dt_data)))
    ax2.plot(t[::stride], (x_num * 100)[::stride],
             linewidth=0, c="tab:blue", marker="o", markersize=5,
             markerfacecolor="tab:blue", zorder=2, label=r"Simulation ($k_t > 0$)")
    ax2.set_xlabel(r"Time, $t$ [s]")
    ax2.set_ylabel(r"Position, $x$ [cm]")
    ax2.set_ylim(-0.05, 1.0)
    ax2.set_xlim(0.0, T_END)
    ax2.set_box_aspect(1)
    handles, labels = ax2.get_legend_handles_labels()
    try:
        idx_anal = labels.index("Analytical")
        idx_sim = next(i for i, l in enumerate(labels) if "simulation" in l.lower())
        ax2.legend([handles[idx_sim], handles[idx_anal]], [labels[idx_sim], labels[idx_anal]], loc="upper left", fontsize=8)
    except Exception:
        ax2.legend(loc="upper left", fontsize=8)
    ax2.grid(color="lightgrey", linestyle="--", linewidth=0.5)
    fig2.tight_layout()
    out2 = os.path.join(plots_root, "SlidingRamp_HookeMemoryPosition.eps")
    plt.savefig(out2, format="eps", bbox_inches="tight")
    plt.close(fig2)
    print(f"Wrote: {out2}")


# ---------------------------------------------------------------------------
def main() -> None:
    parser = argparse.ArgumentParser(
        description="Position plots for the slidingRamp test (Hooke vs HookeMemory)"
    )
    parser.add_argument("--result-dir", default="./results/slidingRamp")
    parser.add_argument("--plot-dir",   default="./plots/slidingRamp")
    args = parser.parse_args()

    print("Sliding Ramp -- physical parameters:")
    print(f"  R         = {R} m  (sphere radius)")
    print(f"  rho       = {RHO} kg/m^3  =>  m = {M:.5f} kg")
    print(f"  theta     = {THETA} deg  =>  tan(theta) = {np.tan(THETA_R):.4f}")
    print(f"  kn        = {KN:.1e} N/m")
    print(f"  kt        = {KT:.1e} N/m  (HookeMemory only)")
    print(f"  muc       = {MUC}")
    print(f"  G_tan     = {G_TAN:.4f} m/s^2  (drives sliding in X)")
    print(f"  G_nor     = {G_NOR:.4f} m/s^2  (presses into floor)")
    print(f"  F_n       = {F_N:.4f} N   (static normal force)")
    print(f"  Coulomb   = muc*F_n = {_coulomb:.4f} N")
    print(f"  x_eq      = {X_EQ*1e3:.4f} mm  (HookeMemory spring equilibrium)")
    print(f"  omega_t   = {OMEGA_T:.4f} rad/s  =>  T_osc = {T_OSC*1e3:.2f} ms")
    print(f"  x_peak    = {2*X_EQ*1e3:.4f} mm  (max spring force = {_ft_peak:.2f} N < {_coulomb:.2f} N)")
    print(f"  Hooke x(T_END) = {0.5*G_TAN*T_END**2*100:.2f} cm  (slides {0.5*G_TAN*T_END**2*100:.1f} cm)")
    print()

    plot_comparison(args.result_dir, args.plot_dir)


if __name__ == "__main__":
    main()
