#!/usr/bin/env python3
"""Plot contact-table benchmark latency vs capacity (CPU/GPU) with speedup."""

import os
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.use("Agg")
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica",
    "font.size": 16,
})

SPEEDUP_COLOR = 'tab:blue'

def plot_benchmarks(csv_file='data/contact_table_bench.csv', output_file='data/contact_table_bench.png'):
    """Plot one figure per load factor: latency (log2) and speedup vs capacity."""
    df = pd.read_csv(csv_file)

    load_factors = sorted(df['loadFactor'].unique())

    all_speedups = []
    for lf in load_factors:
        df_lf_tmp = df[df['loadFactor'] == lf].sort_values('capacity')
        all_speedups.extend((df_lf_tmp['host_ms'] / df_lf_tmp['kernel_ms']).values.tolist())
    all_speedups = np.array(all_speedups)
    all_speedups = all_speedups[np.isfinite(all_speedups) & (all_speedups > 0)]
    if all_speedups.size == 0:
        global_speed_min, global_speed_max = 1.0, 1.0
    else:
        global_speed_min = float(np.min(all_speedups)) - 1
        global_speed_max = float(np.max(all_speedups)) + 1

    for idx, lf in enumerate(load_factors):
        fig, ax1 = plt.subplots(figsize=(6, 5))

        df_lf = df[df['loadFactor'] == lf].sort_values('capacity')

        ax1.set_xscale('log', base=2)
        ax1.set_yscale('log', base=2)

        ax1.plot(df_lf['capacity'], df_lf['host_ms'],
                 marker='o', linestyle='-', label='CPU', color='black', linewidth=1, markersize=6, markerfacecolor='black')
        ax1.plot(df_lf['capacity'], df_lf['kernel_ms'],
                 marker='s', linestyle='--', label='GPU', color='black', linewidth=1, markersize=6, markerfacecolor='black')

        ax1.set_xlabel(r'Capacity [-]')
        ax1.set_ylabel(r'Latency [ms]')
        ax1.set_axisbelow(True)
        ax1.grid(color='lightgrey', linestyle='--', linewidth=0.5, which='both', zorder=0)
        ax1.legend(loc='upper left')

        ax2 = ax1.twinx()
        speedup = df_lf['host_ms'] / df_lf['kernel_ms']
        ax2.plot(df_lf['capacity'], speedup,
            marker='o', linestyle='-', color=SPEEDUP_COLOR, linewidth=1, markersize=4, markerfacecolor=SPEEDUP_COLOR)
        ax2.set_ylabel(r'Speedup [-]', color=SPEEDUP_COLOR)
        ax2.tick_params(axis='y', labelcolor=SPEEDUP_COLOR)
        ax2.set_ylim(global_speed_min, global_speed_max)
        ax2.set_xscale('log', base=2)

        out_file = output_file
        base, ext = os.path.splitext(output_file)
        out_file = f"{base}_LF{lf}{ext}"
        fig.tight_layout()
        plt.savefig(out_file, format='eps', bbox_inches='tight')
        plt.close(fig)
        print(f"Plot saved to: {out_file}")

if __name__ == '__main__':
    import sys
    csv_file = sys.argv[1] if len(sys.argv) > 1 else 'data/contact_table_bench.csv'
    output_file = sys.argv[2] if len(sys.argv) > 2 else 'data/contact_table_bench.eps'

    plot_benchmarks(csv_file, output_file)
