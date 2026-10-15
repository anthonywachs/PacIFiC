#!/usr/bin/env python3
"""Plot collision detection benchmark results from comprehensive CSV data."""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({"text.usetex": True, "font.family": "Helvetica", "font.size": 16})

CPU_COLOR = '#4A90E2'
GPU_COLOR = '#1E3A8A'

COMBO_COLORS = [
    '#E3F2FD',
    '#BBDEFB',
    '#90CAF9',
    '#64B5F6',
    '#42A5F5',
    '#2196F3',
    '#1E88E5',
    '#1565C0'
]

def create_shape_label(row):
    """ShapeType + aspect ratio -> short label (S1, B4, ...)."""
    shape = row['ShapeType']
    aspect = int(row['AspectRatio'])

    if shape == 'Sphere':
        return f'S{aspect}'
    elif shape == 'Box':
        return f'B{aspect}'
    elif shape == 'Superquadric':
        return f'SQ{aspect}'
    return shape

def create_combo_code(row):
    """Abbreviate GJK algorithm / representation / transform flags."""
    algo = 'S' if row['GJKAlgo'] == 'SignedVolume' else 'J'

    rep = 'T' if row['GJKRepresentation'] == 'Transform' else 'Q'

    trans = 'R' if row['UseRelativeTransform'] == 1 else 'A'

    return f'{algo}{rep}{trans}'

def find_best_configuration(df, precision, platform, particle_count, shape_label):
    """Median-fastest combo and its time for one (precision, platform, N, shape)."""
    subset = df[
        (df['Precision'] == precision) &
        (df['Platform'] == platform) &
        (df['ParticleCount'] == particle_count) &
        (df['ShapeLabel'] == shape_label)
    ].copy()

    if len(subset) == 0:
        return None, None, None

    subset['ComboCode'] = subset.apply(create_combo_code, axis=1)

    combo_medians = subset.groupby('ComboCode')['TotalTime_ms'].median()

    best_combo = combo_medians.idxmin()
    best_time = combo_medians.min()

    pair_count = subset['PairCount'].iloc[0]

    return best_combo, best_time, pair_count

def plot_precision(csv_file, precision, output_file):
    """Best CPU vs best GPU time per N, 2x2 panels by shape."""
    df = pd.read_csv(csv_file)

    df['ShapeLabel'] = df.apply(create_shape_label, axis=1)

    shapes = ['S1', 'B1', 'SQ4', 'B4']
    shape_titles = {
        'S1': 'Sphere (AR=1)',
        'B1': 'Box (AR=1)',
        'SQ4': 'Superquadric (AR=4)',
        'B4': 'Box (AR=4)'
    }

    particle_counts = sorted(df['ParticleCount'].unique(), reverse=True)

    fig, axes = plt.subplots(2, 2, figsize=(18, 14))
    fig.suptitle(f'Collision Detection Performance - {precision} Precision', fontsize=18, y=0.995)

    all_times = []
    shape_data = {}

    for shape in shapes:
        cpu_times = []
        gpu_times = []
        cpu_labels = []
        gpu_labels = []
        y_labels = []

        for pc in particle_counts:
            cpu_combo, cpu_time, pair_count = find_best_configuration(
                df, precision, 'CPU', pc, shape
            )

            gpu_combo, gpu_time, _ = find_best_configuration(
                df, precision, 'GPU', pc, shape
            )

            if cpu_time is not None and gpu_time is not None:
                cpu_times.append(cpu_time)
                gpu_times.append(gpu_time)
                cpu_labels.append(cpu_combo)

                speedup = cpu_time / gpu_time
                gpu_labels.append(f'{gpu_combo} ({speedup:.1f}x)')

                y_labels.append(f'{pc}\n({pair_count})')

                all_times.extend([cpu_time, gpu_time])

        shape_data[shape] = {
            'cpu_times': cpu_times,
            'gpu_times': gpu_times,
            'cpu_labels': cpu_labels,
            'gpu_labels': gpu_labels,
            'y_labels': y_labels
        }

    if len(all_times) > 0:
        min_time = min(all_times)
        max_time = max(all_times)
        x_min = min_time * 0.5
        x_max = max_time * 5.0
    else:
        x_min, x_max = 0.001, 100

    for idx, shape in enumerate(shapes):
        ax = axes[idx // 2, idx % 2]

        data = shape_data.get(shape)
        if data is None or len(data['cpu_times']) == 0:
            ax.text(0.5, 0.5, 'No data available',
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_title(shape_titles.get(shape, shape), fontweight='bold')
            continue

        cpu_times = data['cpu_times']
        gpu_times = data['gpu_times']
        cpu_labels = data['cpu_labels']
        gpu_labels = data['gpu_labels']
        y_labels = data['y_labels']

        n = len(cpu_times)
        y_pos = np.arange(n) * 2

        bars_cpu = ax.barh(y_pos + 0.9, cpu_times, height=0.8,
                          color=CPU_COLOR, label='CPU', alpha=0.8)

        bars_gpu = ax.barh(y_pos, gpu_times, height=0.8,
                          color=GPU_COLOR, label='GPU', alpha=0.8)

        for i, (bar_cpu, bar_gpu) in enumerate(zip(bars_cpu, bars_gpu)):
            width_cpu = bar_cpu.get_width()
            ax.text(width_cpu * 1.15, bar_cpu.get_y() + bar_cpu.get_height()/2,
                   cpu_labels[i], ha='left', va='center', fontsize=12, fontweight='bold')

            width_gpu = bar_gpu.get_width()
            ax.text(width_gpu * 1.15, bar_gpu.get_y() + bar_gpu.get_height()/2,
                   gpu_labels[i], ha='left', va='center', fontsize=12, fontweight='bold')

        ax.set_yticks(y_pos + 0.9)
        ax.set_yticklabels(y_labels, rotation=90, va='center', ha='center')
        ax.tick_params(axis='y', which='major', pad=12)

        ax.set_xlabel(r'Clock Time [ms]', fontweight='bold')
        ax.set_ylabel(r'No. Particles' + '\n' + r'(Pairs)', fontweight='bold')
        ax.set_title(shape_titles.get(shape, shape), fontweight='bold')

        ax.set_xscale('log')
        ax.set_xlim(x_min, x_max)

        ax.grid(axis='x', color='lightgrey', linestyle='--', linewidth=0.5, alpha=0.7)
        ax.set_axisbelow(True)

        if idx == 0:
            ax.legend(loc='lower right', fontsize=14)

    plt.tight_layout()
    plt.subplots_adjust(left=0.16)
    plt.savefig(output_file, format='eps', bbox_inches='tight')
    print(f"Plot saved to: {output_file}")
    plt.close()

def plot_all_combinations(csv_file, precision, platform, output_file):
    """All eight GJK combos vs time for a few N, 2x2 panels by shape."""
    df = pd.read_csv(csv_file)

    df['ShapeLabel'] = df.apply(create_shape_label, axis=1)
    df['ComboCode'] = df.apply(create_combo_code, axis=1)

    shapes = ['S1', 'B1', 'SQ4', 'B4']
    shape_titles = {
        'S1': 'Sphere (AR=1)',
        'B1': 'Box (AR=1)',
        'SQ4': 'Superquadric (AR=4)',
        'B4': 'Box (AR=4)'
    }

    all_combos = ['JQA', 'JQR', 'JTA', 'JTR', 'SQA', 'SQR', 'STA', 'STR']

    df_filtered = df[(df['Precision'] == precision) & (df['Platform'] == platform)].copy()

    all_particle_counts = sorted(df_filtered['ParticleCount'].unique())
    if len(all_particle_counts) == 0:
        print(f"No data for {platform}-{precision}")
        return

    if len(all_particle_counts) >= 3:
        selected_counts = [
            all_particle_counts[-1],
            all_particle_counts[len(all_particle_counts)//2],
            all_particle_counts[0]
        ]
    else:
        selected_counts = list(reversed(all_particle_counts))

    fig, axes = plt.subplots(2, 2, figsize=(18, 14))
    fig.suptitle(f'GJK Configuration Comparison - {platform} {precision} Precision', fontsize=18, y=0.995)

    all_times = []

    for shape in shapes:
        for pc in selected_counts:
            subset = df_filtered[
                (df_filtered['ShapeLabel'] == shape) &
                (df_filtered['ParticleCount'] == pc)
            ]
            if len(subset) > 0:
                combo_medians = subset.groupby('ComboCode')['TotalTime_ms'].median()
                all_times.extend(combo_medians.values.tolist())

    if len(all_times) > 0:
        min_time = min(all_times)
        max_time = max(all_times)
        x_min = min_time * 0.5
        x_max = max_time * 3.0
    else:
        x_min, x_max = 0.001, 100

    for idx, shape in enumerate(shapes):
        ax = axes[idx // 2, idx % 2]

        all_data = []
        y_labels = []

        for pc in selected_counts:
            subset = df_filtered[
                (df_filtered['ShapeLabel'] == shape) &
                (df_filtered['ParticleCount'] == pc)
            ]

            if len(subset) == 0:
                continue

            pair_count = subset['PairCount'].iloc[0]
            y_labels.append(f'{pc}\n({pair_count})')

            combo_times = []
            for combo in all_combos:
                combo_subset = subset[subset['ComboCode'] == combo]
                if len(combo_subset) > 0:
                    median_time = combo_subset['TotalTime_ms'].median()
                    combo_times.append(median_time)
                else:
                    combo_times.append(np.nan)

            all_data.append(combo_times)

        if len(all_data) == 0:
            ax.text(0.5, 0.5, 'No data available',
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_title(shape_titles.get(shape, shape), fontweight='bold')
            continue

        n_particles = len(all_data)
        n_combos = len(all_combos)
        bar_height = 0.8 / n_combos
        y_positions = np.arange(n_particles) * (n_combos * bar_height + 0.5)

        for combo_idx, combo in enumerate(all_combos):
            combo_values = [all_data[i][combo_idx] for i in range(n_particles)]
            combo_values_safe = [v if not np.isnan(v) and v > 0 else 0.001 for v in combo_values]
            y_pos = y_positions + combo_idx * bar_height

            ax.barh(y_pos, combo_values_safe, height=bar_height * 0.9,
                   color=COMBO_COLORS[combo_idx], label=combo, alpha=0.9)

        ax.set_yticks(y_positions + (n_combos * bar_height) / 2 - bar_height / 2)
        ax.set_yticklabels(y_labels, rotation=90, va='center', ha='center')
        ax.tick_params(axis='y', which='major', pad=12)

        ax.set_xlabel(r'Clock Time [ms]', fontweight='bold')
        ax.set_ylabel(r'No. Particles' + '\n' + r'(Pairs)', fontweight='bold')
        ax.set_title(shape_titles.get(shape, shape), fontweight='bold')

        ax.set_xscale('log')
        ax.set_xlim(x_min, x_max)

        ax.grid(axis='x', color='lightgrey', linestyle='--', linewidth=0.5, alpha=0.7)
        ax.set_axisbelow(True)

        if idx == 0:
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(reversed(handles), reversed(labels), loc='lower right', fontsize=10, ncol=1)

    plt.tight_layout()
    plt.subplots_adjust(left=0.16)
    plt.savefig(output_file, format='eps', bbox_inches='tight')
    print(f"Plot saved to: {output_file}")
    plt.close()

def plot_benchmarks(csv_file='data/collision_benchmark_comprehensive.csv'):
    """Write optimal-config and all-combo EPS figures for each precision/platform."""
    try:
        df = pd.read_csv(csv_file)
    except FileNotFoundError:
        print(f"Error: Could not find {csv_file}")
        return

    precisions = df["Precision"].unique()

    print("\n" + "="*80)
    print("Generating Optimal Configuration Plots")
    print("="*80)
    for precision in precisions:
        output_file = f'data/collision_benchmark_{precision.lower()}.eps'
        print(f"\nGenerating plot for {precision} precision...")
        plot_precision(csv_file, precision, output_file)

    print("\n" + "="*80)
    print("Generating All Combinations Comparison Plots")
    print("="*80)
    platforms = ['CPU', 'GPU']
    for platform in platforms:
        for precision in precisions:
            output_file = f'data/collision_combos_{platform.lower()}_{precision.lower()}.eps'
            print(f"\nGenerating plot for {platform} {precision} precision...")
            plot_all_combinations(csv_file, precision, platform, output_file)

if __name__ == '__main__':
    import sys
    csv_file = sys.argv[1] if len(sys.argv) > 1 else 'data/collision_benchmark_comprehensive.csv'

    plot_benchmarks(csv_file)
