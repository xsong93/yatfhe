import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import json
import sys
import os
import argparse
import numpy as np
from scipy import stats
from scipy.cluster.vq import kmeans2


def analyze_data(data, method_name):
    """Analyze benchmark data and calculate statistics"""
    times_ms = [item['time_us'] / 1000 for item in data['iterations']]
    data_mean = np.mean(times_ms)

    return {
        'method': method_name,
        'all_times': times_ms,
        'avg': data_mean,
        'pressure_ratio': data.get('pressure_ratio', 1)  # Add default value
    }


def plot_analysis(json_file_paths, output_filename="fig_avg_analysis"):
    """Plot benchmark analysis chart"""

    # Check if files exist
    existing_files = [f for f in json_file_paths if os.path.exists(f)]
    if not existing_files:
        print("No valid files found")
        return

    # Analyze each dataset
    analyses = []
    for file_path in existing_files:
        try:
            with open(file_path, 'r', encoding='utf-8') as file:
                data = json.load(file)

            # Get method name
            method_name = data.get('benchmark_name', 'Unknown').split('/')[-1]

            # Use fixed version analysis
            analysis = analyze_data(data, method_name)
            analyses.append(analysis)

        except Exception as e:
            print(f"Error processing {file_path}: {e}")
            continue

    if not analyses:
        print("No valid analyses to plot")
        return

    # Group data by method
    data_by_method = {}
    for analysis in analyses:
        method = analysis['method']
        if method not in data_by_method:
            data_by_method[method] = {'pressure_ratios': [], 'averages': []}

        data_by_method[method]['pressure_ratios'].append(analysis['pressure_ratio'])
        data_by_method[method]['averages'].append(analysis['avg'])

    # Sort data for each method by pressure ratio
    for method in data_by_method:
        pressure_ratios = data_by_method[method]['pressure_ratios']
        averages = data_by_method[method]['averages']

        # Sort by pressure ratio
        sorted_indices = np.argsort(pressure_ratios)
        data_by_method[method]['pressure_ratios'] = [pressure_ratios[i] for i in sorted_indices]
        data_by_method[method]['averages'] = [averages[i] for i in sorted_indices]

    fig, ax = plt.subplots(figsize=(6, 4))

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

    for i, (method, data) in enumerate(data_by_method.items()):
        pressure_ratios = data['pressure_ratios']
        averages = data['averages']

        # Plot all points with connecting lines
        ax.plot(pressure_ratios, averages,
                 marker='o', linestyle='-', linewidth=1, markersize=2,
                 color=colors[i % len(colors)], alpha=0.7,
                 label=method)

        # Add value labels on each data point
        for j, (x, y) in enumerate(zip(pressure_ratios, averages)):
            ax.annotate(f'{y:.0f}',
                         xy=(x, y),
                         xytext=(0, 10),  # Offset 10 points upward
                         textcoords='offset points',
                         ha='center', va='bottom',
                         fontsize=6,
                         color=colors[i % len(colors)])

    # ax.set_title('Average Execution Time under Different Pressure Ratio',
    #               fontsize=8, fontweight='bold', pad=20)
    ax.set_xlabel('Pressure Ratio', fontsize=12)
    ax.set_ylabel('Average Execution Time (milliseconds)', fontsize=12)
    # ax.grid(True, alpha=0.3, linestyle='--')

    # Set x-axis ticks
    unique_pressure_ratios = sorted(set([ratio for data in data_by_method.values()
                                         for ratio in data['pressure_ratios']]))
    ax.set_xticks(unique_pressure_ratios)
    ax.set_xticklabels([f'{ratio}x' for ratio in unique_pressure_ratios],
                       rotation=60,
                       ha='center',
                       va='top')

    ax.tick_params(axis='both', labelsize=6)

    # Add legend
    ax.legend(loc='best', fontsize=8, frameon=True, fancybox=True, framealpha=0.8)

    plt.tight_layout()
    plt.savefig(output_filename + '.png', dpi=300, bbox_inches='tight')
    # plt.show()

    # Save as PDF (vector format) with tight bounding box and no padding
    plt.savefig(output_filename + '.pdf', format='pdf', dpi=300,
                bbox_inches='tight', pad_inches=0)

    return analyses


# Usage example
if __name__ == "__main__":

    # Define files for all pressure ratios
    pressure_ratios = [1, 5, 10, 25, 100]

    json_files = []
    for ratio in pressure_ratios:
        json_files.append(f'server/tfhe_benchmark_results_{ratio}.json')
        json_files.append(f'server/wwl+24_benchmark_results_{ratio}.json')
        json_files.append(f'server/ours_benchmark_results_{ratio}.json')

    plot_analysis(json_files, 'fig_avg_analysis')