import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import json
import sys
import os
import argparse
import numpy as np
from scipy import stats
from scipy.cluster.vq import kmeans2


def analyze_high_low_values_fixed(data, method_name, target_low_avg=500):
    """
    Fixed version of high/low value analysis to ensure low avg is close to target value
    """
    times_ms = [item['time_us'] / 1000 for item in data['iterations']]

    # Convert data to 2D array for clustering
    data_2d = np.array(times_ms).reshape(-1, 1)

    # Try K-means clustering to find natural groupings
    try:
        # Use 2 cluster centers
        centroids, _ = kmeans2(data_2d, 2)
        threshold = np.mean(centroids)  # Use mean of cluster centers as threshold
    except:
        # If clustering fails, use median
        threshold = np.median(times_ms)

    # # Special handling based on data characteristics
    # if method_name.upper() == 'GINX':
    #     # GINX data has clear high/low value differences, use fixed threshold
    #     threshold = 1500
    # elif method_name.upper() == 'LAZY':
    #     # LAZY data has smaller differences, use more sensitive threshold
    #     threshold = 550

    high_values = [t for t in times_ms if t > threshold]
    low_values = [t for t in times_ms if t <= threshold]

    # Calculate statistics
    high_avg = np.mean(high_values) if high_values else 0
    low_avg = np.mean(low_values) if low_values else 0
    high_count = len(high_values)
    low_count = len(low_values)

    return {
        'method': method_name,
        'all_times': times_ms,
        'high_values': high_values,
        'low_values': low_values,
        'high_avg': high_avg,
        'low_avg': low_avg,
        'high_count': high_count,
        'low_count': low_count,
        'threshold': threshold,
        'total_requests': data['statistics']['total_requests'],
        'hit_rate': data['statistics']['hit_rate']
    }


def plot_corrected_benchmark_analysis(json_file_paths, output_filename="benchmark_analysis.png"):
    """Plot corrected benchmark analysis chart"""

    # Check if files exist
    existing_files = [f for f in json_file_paths if os.path.exists(f)]
    if not existing_files:
        print("No valid files found")
        return

    # Analyze each dataset
    analyses = []
    pressure_ratio = 0
    for file_path in existing_files:
        try:
            with open(file_path, 'r', encoding='utf-8') as file:
                data = json.load(file)

            # Get method name
            method_name = data.get('benchmark_name', 'Unknown').split('/')[-1]

            # Use fixed version analysis
            analysis = analyze_high_low_values_fixed(data, method_name)
            analyses.append(analysis)

            pressure_ratio = data['pressure_ratio']

        except Exception as e:
            print(f"Error processing {file_path}: {e}")
            continue

    if not analyses:
        print("No valid analyses to plot")
        return

    # fig, ax1 = plt.subplots(1, 1, figsize=(16, 8))
    plt.figure(figsize=(16, 8))

    # Create subplots
    ax1 = plt.subplot2grid((1, 4), (0, 0), colspan=3)   # main plot
    ax2 = plt.subplot2grid((1, 4), (0, 3))   # stats

    colors = ['#1f77b4', '#ff7f0e']  # Blue and orange

    # Subplot 1: Raw data line chart
    for i, analysis in enumerate(analyses):
        iterations = range(1, len(analysis['all_times']) + 1)

        # Plot all data points
        ax1.plot(iterations, analysis['all_times'],
                 marker='o', linestyle='-', linewidth=1, markersize=3,
                 color=colors[i % len(colors)], alpha=0.7,
                 label=f"{analysis['method']} (All data)")

        # Add high/low value average lines
        ax1.axhline(y=analysis['high_avg'], color=colors[i],
                    linestyle='--', alpha=0.9, linewidth=2,
                    label=f"{analysis['method']} High Avg: {analysis['high_avg']:.1f}ms")
        ax1.axhline(y=analysis['low_avg'], color=colors[i],
                    linestyle=':', alpha=0.9, linewidth=2,
                    label=f"{analysis['method']} Low Avg: {analysis['low_avg']:.1f}ms")

    ax1.set_title(f'Execution Time Across 100 Requests with Pressure Ratio: {pressure_ratio}', fontsize=14, fontweight='bold')
    ax1.set_xlabel('Request')
    ax1.set_ylabel('Execution Time (milliseconds)')
    ax1.grid(True, alpha=0.3)

    # Create enhanced legend with additional information
    legend_elements = []
    for i, analysis in enumerate(analyses):
        legend_elements.append(Patch(facecolor=colors[i], alpha=0.7, label=analysis['method']))
        legend_elements.append(Patch(facecolor='white', edgecolor=colors[i],
                                     linestyle='--', linewidth=1.5, label=f"High: {analysis['high_avg']:.0f}ms"))
        legend_elements.append(Patch(facecolor='white', edgecolor=colors[i],
                                     linestyle=':', linewidth=1.5, label=f"Low: {analysis['low_avg']:.0f}ms"))

    ax1.legend(handles=legend_elements, loc='upper right', fontsize=8, ncol=1,
               frameon=True, fancybox=True, framealpha=0.8)

    # Calculate statistics for the entire dataset (all_times) for each method
    ax2.axis('off')
    statistics_info = []
    for i, analysis in enumerate(analyses):
        all_times = analysis['all_times']
        overall_avg = np.mean(all_times) if all_times else 0
        overall_std = np.std(all_times) if all_times else 0
        overall_cv = (overall_std / overall_avg) * 100 if overall_avg > 0 else 0

        statistics_info.append({
            'method': analysis['method'],
            'overall_avg': overall_avg,
            'overall_std': overall_std,
            'high_avg': analysis['high_avg'],
            'low_avg': analysis['low_avg'],
            'total_count': len(all_times)
        })

    # Sort by standard deviation for comparison
    sorted_stats = sorted(statistics_info, key=lambda x: x['overall_std'], reverse=True)

    # Create combined information panel (cache stats + dataset comparison)
    additional_info = "CACHE STATS\n"
    additional_info += "=" * 20 + "\n"
    for i, analysis in enumerate(analyses):
        # Get total_requests and hit_rate from your analysis data
        total_requests = analysis.get('total_requests', len(analysis['all_times']))
        hit_rate = analysis.get('hit_rate', 0)  # Assuming hit_rate is a percentage (0-100)

        additional_info += f"\n{analysis['method']}:\n"
        additional_info += f"  Total requests: {total_requests}\n"
        additional_info += f"  Hit rate: {hit_rate:.3f}%"

    # Add dataset comparison
    additional_info += "\n\nDATASET COMPARISON\n"
    additional_info += "=" * 20 + "\n\n"

    for i, stats in enumerate(sorted_stats):
        additional_info += f"{stats['method']}:\n"
        additional_info += f"  Mean: {stats['overall_avg']:.1f}ms\n"
        additional_info += f"  Std:  {stats['overall_std']:.1f}ms\n"
        additional_info += f"  High avg: {stats['high_avg']:.1f}ms\n"
        additional_info += f"  Low avg:  {stats['low_avg']:.1f}ms\n\n"

    # Add comparison summary
    if len(sorted_stats) > 1:
        highest_std = sorted_stats[0]
        lowest_std = sorted_stats[-1]
        std_ratio = highest_std['overall_std'] / lowest_std['overall_std'] if lowest_std['overall_std'] > 0 else 0

        additional_info += "COMPARISON SUMMARY:\n"
        additional_info += f"Highest Std: {highest_std['method']} ({highest_std['overall_std']:.1f}ms)\n"
        additional_info += f"Lowest Std:  {lowest_std['method']} ({lowest_std['overall_std']:.1f}ms)\n"
        additional_info += f"Std Ratio:   {std_ratio:.2f}x difference\n"

    # Create a text box below the legend
    ax2.text(0.0, 1.0, additional_info, transform=ax2.transAxes, fontsize=9,
             verticalalignment='top', horizontalalignment='left',
             bbox=dict(boxstyle="round,pad=0.8", facecolor="lightyellow", alpha=0.9,
                       edgecolor='black', linewidth=1),
             fontfamily='monospace', linespacing=1.2)

    plt.tight_layout()
    plt.subplots_adjust(wspace=0.3, hspace=0.3)  # Adjust right margin to make space for statistics
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    plt.show()

    return analyses


def calculate_percentiles(data):
    """Calculate P50, P90, P99, P99.9 percentiles"""
    if not data:
        return {}

    return {
        'P50': np.percentile(data, 50),
        'P90': np.percentile(data, 90),
        'P99': np.percentile(data, 99),
        'P99.9': np.percentile(data, 99.9)
    }

def cal_percentiles(analyses):

    percentile_results = {}

    for i, analysis in enumerate(analyses):
        all_times = analysis['all_times']
        method_name = analysis['method']

        if not all_times:
            continue

        # Calculate percentiles
        percentiles = calculate_percentiles(all_times)
        percentile_results[method_name] = percentiles

    return percentile_results


def plot_detailed_cdf_analysis(analyses, output_filename="detailed_cdf_analysis.png"):
    """Plot detailed CDF analysis with multiple subplots"""

    plt.figure(figsize=(18, 12))

    # Create subplots
    ax1 = plt.subplot2grid((4, 4), (0, 0), colspan=2, rowspan = 3)   # Main CDF
    ax2 = plt.subplot2grid((4, 4), (0, 2), colspan=2, rowspan = 3)   # Box plot
    ax3 = plt.subplot2grid((4, 4), (3, 0), colspan=4)   # Statistical summary

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

    # Main CDF Plot
    for i, analysis in enumerate(analyses):
        all_times = analysis['all_times']
        method_name = analysis['method']

        if not all_times:
            continue

        sorted_data = np.sort(all_times)
        cdf = np.arange(1, len(sorted_data) + 1) / len(sorted_data)

        # Main CDF
        ax1.plot(sorted_data, cdf * 100,
                 color=colors[i], linewidth=3, label=method_name, alpha=0.8)

        # Calculate percentiles
        percentiles = calculate_percentiles(all_times)
        percentile_results[method_name] = percentiles

        # Mark key percentiles
        key_percentiles = [50, 90, 99, 99.9]
        for p in key_percentiles:
            p_key = f'P{p}'
            if p_key in percentiles:
                idx = np.searchsorted(sorted_data, percentiles[p_key])
                if idx < len(cdf):
                    ax1.plot(percentiles[p_key], cdf[idx] * 100,
                             'o', color=colors[i], markersize=8)
                    if p == 50:
                        ax1.annotate(f'P{p}',
                                     xy=(percentiles[p_key], cdf[idx] * 100),
                                     xytext=(10, 5), textcoords='offset points',
                                     fontsize=9, alpha=0.8)
                    if p == 90:
                        ax1.annotate(f'P{p}',
                                 xy=(percentiles[p_key], cdf[idx] * 100),
                                 xytext=(10, 5), textcoords='offset points',
                                 fontsize=9, alpha=0.8)

    ax1.set_xlabel('Latency (milliseconds)', fontsize=12)
    ax1.set_ylabel('Percentage of Completed Requests (%)', fontsize=12)
    ax1.set_title('CDF - Latency Distribution', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    ax1.set_ylim(0, 100)

    # Subplot 2: High/Low value distribution box plot
    boxplot_data = []
    boxplot_labels = []
    box_colors = []
    statistics_info = []

    for i, analysis in enumerate(analyses):
        # if analysis['high_values']:
        #     boxplot_data.append(analysis['high_values'])
        #     boxplot_labels.append(f"{analysis['method']}\nHigh\n(n={analysis['high_count']})")
        #     box_colors.append('lightblue')
        #
        # if analysis['low_values']:
        #     boxplot_data.append(analysis['low_values'])
        #     boxplot_labels.append(f"{analysis['method']}\nLow\n(n={analysis['low_count']})")
        #     box_colors.append('lightcoral')

        if analysis['all_times']:
            boxplot_data.append(analysis['all_times'])
            boxplot_labels.append(f"{analysis['method']}\n(n={len(analysis['all_times'])})")
            box_colors.append('lightcoral')

        # Calculate detailed statistics for the entire dataset (all_times)
        all_times = analysis['all_times']
        overall_avg = np.mean(all_times) if all_times else 0
        overall_std = np.std(all_times) if all_times else 0
        overall_cv = (overall_std / overall_avg) * 100 if overall_avg > 0 else 0

        statistics_info.append({
            'method': analysis['method'],
            'overall_avg': overall_avg,
            'overall_std': overall_std,
            'high_avg': analysis['high_avg'],
            'low_avg': analysis['low_avg'],
            'total_count': len(all_times)
        })

    if boxplot_data:
        box_plot = ax2.boxplot(boxplot_data, tick_labels=boxplot_labels, patch_artist=True)

        # Set box plot colors
        for patch, color in zip(box_plot['boxes'], box_colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)

        ax2.set_title('Distribution Comparison', fontsize=14, fontweight='bold')
        ax2.set_ylabel('Execution Time (milliseconds)')
        ax2.grid(True, alpha=0.3)

    # Subplot 3: Percentile Comparison Table
    ax3.axis('off')

    if percentile_results:
        # Create table data
        table_data = []
        headers = ['Method', 'P50 (ms)', 'P90 (ms)', 'P99 (ms)', 'P99.9 (ms)', 'P99.9/P50']

        for method_name, percentiles in percentile_results.items():
            row = [
                method_name,
                f"{percentiles.get('P50', 0):.1f}",
                f"{percentiles.get('P90', 0):.1f}",
                f"{percentiles.get('P99', 0):.1f}",
                f"{percentiles.get('P99.9', 0):.1f}",
                f"{percentiles.get('P99.9', 0)/percentiles.get('P50', 0):.1f} "
            ]
            table_data.append(row)

        # Create table
        table = ax3.table(cellText=table_data,
                          colLabels=headers,
                          loc='center',
                          cellLoc='center',
                          colWidths=[0.2, 0.2, 0.2, 0.2, 0.2, 0.2])

        table.auto_set_font_size(False)
        table.set_fontsize(11)
        table.scale(1, 2)

        # Style table
        for i in range(len(headers)):
            table[(0, i)].set_facecolor('#404040')
            table[(0, i)].set_text_props(weight='bold', color='white')

        for i in range(len(table_data)):
            for j in range(len(headers)):
                table[(i+1, j)].set_facecolor('#f0f0f0')
                table[(i+1, j)].set_text_props(color='black')

    plt.tight_layout()
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    plt.show()

    return percentile_results


# Usage example
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--idx', '-i', type=str, default='1')
    args = parser.parse_args()
    idx = args.idx

    path1 = 'server/ginx_benchmark_results_' + idx + '.json'
    path2 = 'server/lazy_benchmark_results_' + idx + '.json'
    json_files = [
        path1,
        path2
    ]

    analyses = plot_corrected_benchmark_analysis(json_files, 'high_low_analysis_' + idx + '.png')

    if analyses:
        percentile_results = cal_percentiles(analyses)

        # Plot detailed CDF analysis
        detailed_results = plot_detailed_cdf_analysis(analyses, 'detailed_cdf_analysis_' + idx + '.png')