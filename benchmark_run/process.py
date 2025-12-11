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


def plot_analysis(json_file_paths, output_filename="avg_analysis.png"):
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

    # Create figure with two subplots: left for main plot, right for statistics
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8), gridspec_kw={'width_ratios': [3, 1]})

    colors = ['#1f77b4', '#ff7f0e']  # Blue and orange

    # Main plot (left)
    for i, (method, data) in enumerate(data_by_method.items()):
        pressure_ratios = data['pressure_ratios']
        averages = data['averages']

        # Plot all points with connecting lines
        ax1.plot(pressure_ratios, averages,
                 marker='o', linestyle='-', linewidth=2, markersize=8,
                 color=colors[i % len(colors)], alpha=0.7,
                 label=method)

        # Add value labels on each data point
        for j, (x, y) in enumerate(zip(pressure_ratios, averages)):
            ax1.annotate(f'{y:.0f}',
                         xy=(x, y),
                         xytext=(0, 10),  # Offset 10 points upward
                         textcoords='offset points',
                         ha='center', va='bottom',
                         fontsize=9,
                         color=colors[i % len(colors)])

    ax1.set_title('Average Execution Time under Different Pressure Ratio',
                  fontsize=16, fontweight='bold', pad=20)
    ax1.set_xlabel('Pressure Ratio', fontsize=12)
    ax1.set_ylabel('Average Execution Time (milliseconds)', fontsize=12)
    ax1.grid(True, alpha=0.3, linestyle='--')

    # Set x-axis ticks
    unique_pressure_ratios = sorted(set([ratio for data in data_by_method.values()
                                         for ratio in data['pressure_ratios']]))
    ax1.set_xticks(unique_pressure_ratios)
    ax1.set_xticklabels([f'{ratio}x' for ratio in unique_pressure_ratios])

    # Add legend
    ax1.legend(loc='best', fontsize=10, frameon=True, fancybox=True, framealpha=0.8)

    # Statistics panel (right)
    ax2.axis('off')  # Turn off axis for the statistics panel

    # Create statistical information
    stats_text = "PERFORMANCE TREND ANALYSIS\n"
    stats_text += "=" * 30 + "\n\n"
    slopes = []

    for method, data in data_by_method.items():
        pressure_ratios = data['pressure_ratios']
        averages = data['averages']

        if len(averages) > 1:
            # Calculate trend (slope)
            slope, intercept, r_value, p_value, std_err = stats.linregress(pressure_ratios, averages)
            slopes.append(slope)
            trend = "↗ Increasing" if slope > 0 else ("↘ Decreasing" if slope < 0 else "→ Stable")

            stats_text += f"{method}:\n"
            stats_text += f"  Slope: {slope:.2f}\n"
            stats_text += f"  Min: {min(averages):.0f} ms\n"
            stats_text += f"  Max: {max(averages):.0f} ms\n"
        else:
            slopes.append(0)  # Default if not enough data

    # Add comparison summary
    if len(data_by_method) > 1:
        # Calculate which method has better performance
        methods = list(data_by_method.keys())

        # Method with smallest slope (best)
        best_slope_idx = np.argmin(slopes)
        best_method = methods[best_slope_idx]
        best_slope = slopes[best_slope_idx]

        # Method with largest slope (worst)
        worst_slope_idx = np.argmax(slopes)
        worst_method = methods[worst_slope_idx]
        worst_slope = slopes[worst_slope_idx]

        # Calculate slope ratio
        if best_slope != 0:
            slope_ratio = abs(worst_slope / best_slope)
        else:
            slope_ratio = float('inf') if worst_slope != 0 else 1.0

        stats_text += "\nCOMPARISON SUMMARY\n"
        stats_text += "=" * 20 + "\n\n"
        stats_text += f"Best: {best_method}\n"
        stats_text += f"Worst: {worst_method}\n"
        stats_text += f"Slope Ratio: {slope_ratio:.2f}x\n"

    # Add the statistics text to the right panel
    ax2.text(0.0, 0.98, stats_text, transform=ax2.transAxes, fontsize=9,
             verticalalignment='top', horizontalalignment='left',
             bbox=dict(boxstyle="round,pad=0.5", facecolor="lightyellow", alpha=0.9,
                       edgecolor='black', linewidth=1),
             fontfamily='monospace', linespacing=1.2)

    plt.tight_layout()
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    plt.show()

    return analyses


# Usage example
if __name__ == "__main__":

    # Define files for all pressure ratios
    pressure_ratios = [1, 5, 10, 25, 100]

    json_files = []
    for ratio in pressure_ratios:
        json_files.append(f'server/ginx_benchmark_results_{ratio}.json')
        json_files.append(f'server/lazy_benchmark_results_{ratio}.json')

    analyses = plot_analysis(json_files)

    # Print detailed statistical information
    print("="*50)
    print("PERFORMANCE ANALYSIS SUMMARY")
    print("="*50)

    # Group by method
    ginx_data = []
    lazy_data = []

    for analysis in analyses:
        if 'ginx' in analysis['method'].lower():
            ginx_data.append((analysis['pressure_ratio'], analysis['avg']))
        elif 'lazy' in analysis['method'].lower():
            lazy_data.append((analysis['pressure_ratio'], analysis['avg']))

    # Sort
    ginx_data.sort(key=lambda x: x[0])
    lazy_data.sort(key=lambda x: x[0])

    print("\nGINX Performance:")
    print("-" * 30)
    for ratio, avg in ginx_data:
        print(f"  Pressure {ratio}x: {avg:.0f}ms")

    if len(ginx_data) > 1:
        ratios = [x[0] for x in ginx_data]
        avgs = [x[1] for x in ginx_data]
        slope, intercept, r_value, p_value, std_err = stats.linregress(ratios, avgs)
        print(f"\n  Trend: {'Increasing' if slope > 0 else 'Decreasing' if slope < 0 else 'Stable'}")
        print(f"  Slope: {slope:.2f} ms per pressure unit")

    print("\nLAZY Performance:")
    print("-" * 30)
    for ratio, avg in lazy_data:
        print(f"  Pressure {ratio}x: {avg:.0f}ms")

    if len(lazy_data) > 1:
        ratios = [x[0] for x in lazy_data]
        avgs = [x[1] for x in lazy_data]
        slope, intercept, r_value, p_value, std_err = stats.linregress(ratios, avgs)
        print(f"\n  Trend: {'Increasing' if slope > 0 else 'Decreasing' if slope < 0 else 'Stable'}")
        print(f"  Slope: {slope:.2f} ms per pressure unit")