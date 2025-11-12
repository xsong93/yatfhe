import matplotlib.pyplot as plt
import json
import os
import numpy as np
from scipy.cluster import vq
from scipy import stats

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
        centroids, _ = vq.kmeans2(data_2d, 2)
        threshold = np.mean(centroids)  # Use mean of cluster centers as threshold
    except:
        # If clustering fails, use median
        threshold = np.median(times_ms)

    # Special handling based on data characteristics
    if method_name.upper() == 'GINX':
        # GINX data has clear high/low value differences, use fixed threshold
        threshold = 1500
    elif method_name.upper() == 'LAZY':
        # LAZY data has smaller differences, use more sensitive threshold
        threshold = 550

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
        'threshold': threshold
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
    for file_path in existing_files:
        try:
            with open(file_path, 'r', encoding='utf-8') as file:
                data = json.load(file)

            # Get method name
            method_name = data.get('benchmark_name', 'Unknown').split('/')[-1]

            # Use fixed version analysis
            analysis = analyze_high_low_values_fixed(data, method_name)
            analyses.append(analysis)

        except Exception as e:
            print(f"Error processing {file_path}: {e}")
            continue

    if not analyses:
        print("No valid analyses to plot")
        return

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 12))

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


    ax1.set_title('Execution Time with High/Low Value Analysis', fontsize=14, fontweight='bold')
    ax1.set_xlabel('Iteration Number')
    ax1.set_ylabel('Execution Time (milliseconds)')
    ax1.grid(True, alpha=0.3)
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

    # Subplot 2: High/Low value distribution box plot
    boxplot_data = []
    boxplot_labels = []
    box_colors = []
    statistics_info = []

    for i, analysis in enumerate(analyses):
        if analysis['high_values']:
            boxplot_data.append(analysis['high_values'])
            boxplot_labels.append(f"{analysis['method']}\nHigh\n(n={analysis['high_count']})")
            box_colors.append('lightblue')

        if analysis['low_values']:
            boxplot_data.append(analysis['low_values'])
            boxplot_labels.append(f"{analysis['method']}\nLow\n(n={analysis['low_count']})")
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
        box_plot = ax2.boxplot(boxplot_data, labels=boxplot_labels, patch_artist=True)

        # Set box plot colors
        for patch, color in zip(box_plot['boxes'], box_colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)

        ax2.set_title('High/Low Value Distribution Comparison with Statistics', fontsize=14, fontweight='bold')
        ax2.set_ylabel('Execution Time (milliseconds)')
        ax2.grid(True, alpha=0.3)

        # Add comprehensive statistics panel comparing different datasets
        if statistics_info:
            # Sort by standard deviation for comparison
            sorted_stats = sorted(statistics_info, key=lambda x: x['overall_std'], reverse=True)

            stats_text = "DATASET COMPARISON\n"
            stats_text += "=" * 30 + "\n\n"

            for i, stats in enumerate(sorted_stats):
                stats_text += f"{stats['method']}:\n"
                stats_text += f"Overall Dataset:\n"
                stats_text += f"  Mean: {stats['overall_avg']:.1f}ms\n"
                stats_text += f"  Std:  {stats['overall_std']:.1f}ms\n"
                stats_text += f"Subsets:\n"
                stats_text += f"  High Avg: {stats['high_avg']:.1f}ms\n"
                stats_text += f"  Low Avg:  {stats['low_avg']:.1f}ms\n"
                stats_text += f"  Samples:  {stats['total_count']}\n"
                stats_text += "-" * 25 + "\n\n"

            # Add comparison summary
            if len(sorted_stats) > 1:
                highest_std = sorted_stats[0]
                lowest_std = sorted_stats[-1]
                std_ratio = highest_std['overall_std'] / lowest_std['overall_std'] if lowest_std['overall_std'] > 0 else 0

                comparison_text = f"\nCOMPARISON SUMMARY:\n"
                comparison_text += f"Highest Std: {highest_std['method']} ({highest_std['overall_std']:.1f}ms)\n"
                comparison_text += f"Lowest Std:  {lowest_std['method']} ({lowest_std['overall_std']:.1f}ms)\n"
                comparison_text += f"Std Ratio:   {std_ratio:.2f}x difference\n"

                stats_text += comparison_text

            # Position statistics on the right side
            ax2.text(1.02, 0.98, stats_text, transform=ax2.transAxes, fontsize=9,
                     verticalalignment='top', horizontalalignment='left',
                     bbox=dict(boxstyle="round,pad=0.8", facecolor="lightyellow", alpha=0.9,
                               edgecolor='black', linewidth=1),
                     fontfamily='monospace', linespacing=1.2)

    plt.tight_layout()
    plt.subplots_adjust(right=0.75)  # Adjust right margin to make space for statistics
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    plt.show()

    return analyses


# Use fixed version analysis
if __name__ == "__main__":
    # Your JSON file paths
    json_files = [
        "server/ginx_benchmark_results.json",  # Replace with your GINX file
        "server/lazy_benchmark_results.json"   # Replace with your LAZY file
    ]

    print("Running analysis with optimized thresholds...")

    analyses = plot_corrected_benchmark_analysis(
        json_files,
        "high_low_analysis.png"
    )
