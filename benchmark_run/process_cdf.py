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


def plot_combined_analysis(json_file_paths, output_filename):
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
    fig, ax1 = plt.subplots(figsize=(3, 4))

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

    # 1. Main CDF Plot
    f_size = 8
    marker_size = 0
    line_width = 1
    for i, analysis in enumerate(analyses):
        all_times = analysis['all_times']
        method_name = analysis['method']

        if not all_times:
            continue

        sorted_data = np.sort(all_times)
        cdf = np.arange(1, len(sorted_data) + 1) / len(sorted_data)

        # Main CDF
        ax1.plot(sorted_data, cdf * 100,
                 color=colors[i], linewidth=line_width, label=method_name, alpha=0.8)

        # Calculate and mark key percentiles
        percentiles = calculate_percentiles(all_times)
        key_percentiles = [50, 90, 99, 99.9]

        for p in key_percentiles:
            p_key = f'P{p}'
            if p_key in percentiles:
                idx = np.searchsorted(sorted_data, percentiles[p_key])
                if idx < len(cdf):
                    ax1.plot(percentiles[p_key], cdf[idx] * 100,
                             'o', color=colors[i], markersize=marker_size, alpha=0.8)

    ax1.set_xlabel('Latency (ms)', fontsize=f_size)
    ax1.set_ylabel('Percentage (%)', fontsize=f_size)
    ax1.set_title('CDF - Latency Distribution', fontsize=f_size, fontweight='bold')
    # ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=f_size)
    ax1.set_ylim(0, 100)
    ax1.tick_params(labelsize=f_size)

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.4, wspace=0.4)  # Adjust spacing between subplots

    # Save figure
    plt.savefig(output_filename + '.png', dpi=300, bbox_inches='tight')
    # plt.show()

    # Save as PDF (vector format) with tight bounding box and no padding
    plt.savefig(output_filename + '.pdf', format='pdf', dpi=300,
                bbox_inches='tight', pad_inches=0)

    return analyses


# Usage example
if __name__ == "__main__":

    pressure_ratios = ['1', '5', '10', '25', '50']

    for idx in pressure_ratios:
        path1 = 'server/tfhe_benchmark_results_' + idx + '.json'
        path2 = 'server/wwl+24_benchmark_results_' + idx + '.json'
        path3 = 'server/ours_benchmark_results_' + idx + '.json'
        json_files = [
            path1,
            path2,
            path3
        ]

        plot_combined_analysis(json_files, 'fig_cdf_' + idx)
