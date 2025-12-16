import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import json
import sys
import os
import argparse
import numpy as np
from scipy import stats
from scipy.cluster.vq import kmeans2

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


def plot_combined_analysis(json_file_paths, output_filename="combined_analysis.png"):
    """Plot combined analysis with 3 charts on first row and 1 table on second row"""

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

            pressure_ratio = data.get('pressure_ratio', 0)

        except Exception as e:
            print(f"Error processing {file_path}: {e}")
            continue

    if not analyses:
        print("No valid analyses to plot")
        return

    # Create figure with 2 rows: 3 charts on first row, 1 table on second row
    fig = plt.figure(figsize=(16, 6))

    gs = fig.add_gridspec(2, 6, height_ratios=[4, 1])

    # Row 1: Three charts
    ax1 = fig.add_subplot(gs[0, 0:2])  # Raw data line chart
    ax2 = fig.add_subplot(gs[0, 2:4])  # Main CDF Plot
    ax3 = fig.add_subplot(gs[0, 4:6])  # High/Low value distribution box plot

    # Row 2: One compact table spanning all columns
    ax4 = fig.add_subplot(gs[1, :])  # Statistics table

    colors = ['#1f77b4', '#ff7f0e']  # Blue and orange

    # 1. Raw data line chart
    for i, analysis in enumerate(analyses):
        iterations = range(1, len(analysis['all_times']) + 1)

        # Plot all data points
        ax1.plot(iterations, analysis['all_times'],
                 marker='o', linestyle='-', linewidth=1, markersize=3,
                 color=colors[i], alpha=0.7,
                 label=f"{analysis['method']}")

        # Add high/low value average lines
        ax1.axhline(y=analysis['high_avg'], color=colors[i],
                    linestyle='--', alpha=0.8, linewidth=1.5)
        ax1.axhline(y=analysis['low_avg'], color=colors[i],
                    linestyle=':', alpha=0.8, linewidth=1.5)

    ax1.set_title(f'Raw Data Line Chart\nPressure Ratio: {pressure_ratio}', fontsize=12, fontweight='bold')
    ax1.set_xlabel('Request', fontsize=10)
    ax1.set_ylabel('Execution Time (ms)', fontsize=10)
    ax1.grid(True, alpha=0.3)
    ax1.tick_params(labelsize=9)
    ax1.legend(fontsize=9, loc='upper right')

    # 2. Main CDF Plot
    for i, analysis in enumerate(analyses):
        all_times = analysis['all_times']
        method_name = analysis['method']

        if not all_times:
            continue

        sorted_data = np.sort(all_times)
        cdf = np.arange(1, len(sorted_data) + 1) / len(sorted_data)

        # Main CDF
        ax2.plot(sorted_data, cdf * 100,
                 color=colors[i], linewidth=2, label=method_name, alpha=0.8)

        # Calculate and mark key percentiles
        percentiles = calculate_percentiles(all_times)
        key_percentiles = [50, 90, 99, 99.9]

        for p in key_percentiles:
            p_key = f'P{p}'
            if p_key in percentiles:
                idx = np.searchsorted(sorted_data, percentiles[p_key])
                if idx < len(cdf):
                    ax2.plot(percentiles[p_key], cdf[idx] * 100,
                             'o', color=colors[i], markersize=6, alpha=0.8)

    ax2.set_xlabel('Latency (ms)', fontsize=10)
    ax2.set_ylabel('Percentage (%)', fontsize=10)
    ax2.set_title('CDF - Latency Distribution', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=9)
    ax2.set_ylim(0, 100)
    ax2.tick_params(labelsize=9)

    # 3. Distribution box plot
    boxplot_data = []
    boxplot_labels = []
    box_colors = []

    for i, analysis in enumerate(analyses):
        if analysis['all_times']:
            boxplot_data.append(analysis['all_times'])
            boxplot_labels.append(f"{analysis['method']}\n(n={len(analysis['all_times'])})")
            box_colors.append('lightcoral')

    if boxplot_data:
        box_plot = ax3.boxplot(boxplot_data, tick_labels=boxplot_labels, patch_artist=True)

        # Set box plot colors
        for patch, color in zip(box_plot['boxes'], box_colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)

        ax3.set_title('Distribution Comparison', fontsize=12, fontweight='bold')
        ax3.set_ylabel('Execution Time (ms)', fontsize=10)
        ax3.grid(True, alpha=0.3)
        ax3.tick_params(labelsize=9, rotation=0)

    # 4. Statistics table (second row, spanning all columns)
    ax4.axis('off')

    # Calculate statistics for the table
    statistics_info = []
    percentile_results = {}

    for i, analysis in enumerate(analyses):
        all_times = analysis['all_times']
        method_name = analysis['method']

        # Basic statistics
        overall_avg = np.mean(all_times) if all_times else 0
        overall_std = np.std(all_times) if all_times else 0
        overall_cv = (overall_std / overall_avg) * 100 if overall_avg > 0 else 0

        # Percentiles
        percentiles = calculate_percentiles(all_times)
        percentile_results[method_name] = percentiles

        # Cache stats
        total_requests = analysis.get('total_requests', len(all_times))
        hit_rate = analysis.get('hit_rate', 0)

        statistics_info.append({
            'method': method_name,
            'overall_avg': overall_avg,
            'overall_std': overall_std,
            'overall_cv': overall_cv,
            'high_avg': analysis['high_avg'],
            'low_avg': analysis['low_avg'],
            'total_count': len(all_times),
            'total_requests': total_requests,
            'hit_rate': hit_rate,
            'percentiles': percentiles
        })

    # Create table data
    table_data = []

    # Headers
    headers = [
        'Method',
        'Cache HR (%)',
        'Avg (ms)',
        'Std (ms)',
        'High Avg (ms)',
        'Low Avg (ms)',
        'P50 (ms)',
        'P90 (ms)',
        'P99 (ms)',
        'P99.9 (ms)',
        'P99.9/P50'
    ]

    # Add data for each method
    for stats in statistics_info:
        method = stats['method']
        percentiles = stats['percentiles']
        p50 = percentiles.get('P50', 0)
        p999 = percentiles.get('P99.9', 0)
        p999_p50_ratio = p999 / p50 if p50 > 0 else 0

        row = [
            method,
            f"{stats['hit_rate']:.1f}",
            f"{stats['overall_avg']:.1f}",
            f"{stats['overall_std']:.1f}",
            f"{stats['high_avg']:.1f}",
            f"{stats['low_avg']:.1f}",
            f"{p50:.1f}",
            f"{percentiles.get('P90', 0):.1f}",
            f"{percentiles.get('P99', 0):.1f}",
            f"{p999:.1f}",
            f"{p999_p50_ratio:.1f}"
        ]
        table_data.append(row)

    # Create table
    table = ax4.table(cellText=table_data,
                      colLabels=headers,
                      loc='center',
                      cellLoc='center')

    # Style table
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 2)

    # Style header row
    for i in range(len(headers)):
        table[(0, i)].set_facecolor('#404040')
        table[(0, i)].set_text_props(weight='bold', color='white', fontsize=10)

    # Style data rows
    for i in range(len(table_data)):
        for j in range(len(headers)):
            cell = table[(i+1, j)]
            cell.set_facecolor('#f8f9fa' if i%2==0 else 'white')
            cell.set_text_props(color='black', fontsize=9)

    # Highlight important columns
    highlight_cols = ['Avg (ms)', 'Std (ms)', 'P99.9/P50']
    col_indices = [headers.index(col) for col in highlight_cols if col in headers]

    for i in range(len(table_data)):
        for col_idx in col_indices:
            cell = table[(i+1, col_idx)]
            if headers[col_idx] in ['P99.9/P50']:
                cell.set_text_props(weight='bold', fontsize=9)
            else:
                cell.set_text_props(weight='bold', fontsize=9)

    # Adjust layout
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.4, wspace=0.4)  # Adjust spacing between subplots

    # Save figure
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    plt.show()

    return analyses, statistics_info, percentile_results


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

    analyses, statistics_info, percentile_results = (
        plot_combined_analysis(json_files, 'fig_combined_benchmark_analysis_' + idx + '.png'))
