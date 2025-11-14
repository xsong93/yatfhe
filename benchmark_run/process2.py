import matplotlib.pyplot as plt
import json
import os
import numpy as np
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

    ax1.set_title('Execution Time Across 100 Requests', fontsize=14, fontweight='bold')
    ax1.set_xlabel('Iteration Number')
    ax1.set_ylabel('Execution Time (milliseconds)')
    ax1.grid(True, alpha=0.3)

    # Create enhanced legend with additional information
    legend_handles, legend_labels = ax1.get_legend_handles_labels()

    # Add additional information below the legend
    additional_info = ""
    for i, analysis in enumerate(analyses):
        # Get total_requests and hit_rate from your analysis data
        total_requests = analysis.get('total_requests', len(analysis['all_times']))
        hit_rate = analysis.get('hit_rate', 0)  # Assuming hit_rate is a percentage (0-100)

        additional_info += f"\n{analysis['method']}:\n"
        additional_info += f"  Total requests: {total_requests}\n"
        additional_info += f"  Hit rate: {hit_rate:.3f}%"

    # Create a text box below the legend
    ax1.text(1.02, 0.2, additional_info, transform=ax1.transAxes, fontsize=9,
             verticalalignment='top', horizontalalignment='left',
             bbox=dict(boxstyle="round,pad=0.8", facecolor="lightyellow", alpha=0.9,
                       edgecolor='black', linewidth=1),
             fontfamily='monospace', linespacing=1.2)

    # Adjust legend position to make space for additional info
    ax1.legend(handles=legend_handles, labels=legend_labels,
               bbox_to_anchor=(1.02, 1), loc='upper left')

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

def plot_cdf_with_percentiles(analyses, output_filename="cdf_analysis.png"):
    """Plot CDF graph with percentile markers"""

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']  # Different colors for each dataset

    # Subplot 1: CDF Graph
    percentile_results = {}

    for i, analysis in enumerate(analyses):
        all_times = analysis['all_times']
        method_name = analysis['method']

        if not all_times:
            continue

        # Sort data for CDF
        sorted_data = np.sort(all_times)

        # Calculate CDF values (0% to 100%)
        cdf = np.arange(1, len(sorted_data) + 1) / len(sorted_data)

        # Plot CDF
        ax1.plot(sorted_data, cdf * 100,
                 color=colors[i % len(colors)],
                 linewidth=2,
                 label=f'{method_name}',
                 alpha=0.8)

        # Calculate percentiles
        percentiles = calculate_percentiles(all_times)
        percentile_results[method_name] = percentiles

        # Mark percentile points on CDF
        percentile_values = [50, 90, 99, 99.9]
        for p in percentile_values:
            if p in [50, 90, 99, 99.9]:
                p_key = f'P{p}'
                if p_key in percentiles:
                    # Find the closest CDF value for this percentile
                    idx = np.searchsorted(sorted_data, percentiles[p_key])
                    if idx < len(cdf):
                        ax1.plot(percentiles[p_key], cdf[idx] * 100,
                                 'o', color=colors[i % len(colors)],
                                 markersize=6, alpha=0.8)
                        ax1.annotate(f'P{p}',
                                     xy=(percentiles[p_key], cdf[idx] * 100),
                                     xytext=(10, 5), textcoords='offset points',
                                     fontsize=9, alpha=0.8)

    ax1.set_xlabel('Latency (milliseconds)', fontsize=12)
    ax1.set_ylabel('Percentage of Completed Requests (%)', fontsize=12)
    ax1.set_title('Cumulative Distribution Function (CDF) of Latency',
                  fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    ax1.set_ylim(0, 100)

    # Subplot 2: Percentile Comparison Table
    ax2.axis('off')

    if percentile_results:
        # Create table data
        table_data = []
        headers = ['Method', 'P50 (ms)', 'P90 (ms)', 'P99 (ms)', 'P99.9 (ms)']

        for method_name, percentiles in percentile_results.items():
            row = [
                method_name,
                f"{percentiles.get('P50', 0):.1f}",
                f"{percentiles.get('P90', 0):.1f}",
                f"{percentiles.get('P99', 0):.1f}",
                f"{percentiles.get('P99.9', 0):.1f}"
            ]
            table_data.append(row)

        # Create table
        table = ax2.table(cellText=table_data,
                          colLabels=headers,
                          loc='center',
                          cellLoc='center',
                          colWidths=[0.2, 0.2, 0.2, 0.2, 0.2])

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

        ax2.set_title('Latency Percentile Comparison',
                      fontsize=14, fontweight='bold', pad=20)

    plt.tight_layout()
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    plt.show()

    return percentile_results


def plot_detailed_cdf_analysis(analyses, output_filename="detailed_cdf_analysis.png"):
    """Plot detailed CDF analysis with multiple subplots"""

    fig = plt.figure(figsize=(18, 12))

    # Create subplots
    ax1 = plt.subplot2grid((2, 4), (0, 0), colspan=3)  # Main CDF
    ax2 = plt.subplot2grid((2, 4), (1, 2), colspan=2)   # P90-P99.9 zoom
    ax3 = plt.subplot2grid((2, 4), (0, 3))              # Statistical summary
    ax4 = plt.subplot2grid((2, 4), (1, 0), colspan=2)   # P50-P90 zoom


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

    # -P90 Zoom
    for i, analysis in enumerate(analyses):
        all_times = analysis['all_times']
        if not all_times:
            continue

        sorted_data = np.sort(all_times)
        cdf = np.arange(1, len(sorted_data) + 1) / len(sorted_data)

        # Filter data for P50-90 range
        mask = cdf <= 0.9
        if np.any(mask):
            x_data = sorted_data[mask]
            y_data = cdf[mask] * 100

            # Plot original data
            ax4.plot(x_data, y_data, color=colors[i], linewidth=2,
                     label=analysis['method'], alpha=0.7)

    ax4.set_xlabel('Latency (ms)')
    ax4.set_ylabel('Percentage (%)')
    ax4.set_title('P0 - P90 Range', fontsize=12)
    ax4.grid(True, alpha=0.3)

    # P90-P99.9 Zoom
    for i, analysis in enumerate(analyses):
        all_times = analysis['all_times']
        if not all_times:
            continue

        sorted_data = np.sort(all_times)
        cdf = np.arange(1, len(sorted_data) + 1) / len(sorted_data)

        # Filter data for P90-P99.9 range
        mask = (cdf >= 0.9) & (cdf <= 0.999)
        if np.any(mask):
            x_data = sorted_data[mask]
            y_data = cdf[mask] * 100

            # Plot original data
            ax2.plot(x_data, y_data, color=colors[i], linewidth=2,
                     label=analysis['method'], alpha=0.7)

            # Add fitted line (linear regression)
            if len(x_data) > 1:
                # Perform linear regression
                slope, intercept, r_value, p_value, std_err = stats.linregress(x_data, y_data)

                # Create fitted line
                x_fit = np.linspace(x_data.min(), x_data.max(), 100)
                y_fit = slope * x_fit + intercept

                # Plot fitted line
                ax2.plot(x_fit, y_fit, color=colors[i], linestyle='--',
                         linewidth=1.5, alpha=0.8,
                         label=f"{analysis['method']} fit (k={slope:.4f})")

                # Add slope information as text
                ax2.text(0.05, 0.95 - i*0.1, f"{analysis['method']}: k={slope:.4f}",
                         transform=ax2.transAxes, fontsize=9, color=colors[i],
                         bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

    ax2.set_xlabel('Latency (ms)')
    ax2.set_ylabel('Percentage (%)')
    ax2.set_title('P90 - P99.9 Range (Tail Latency) with Fitted Lines', fontsize=12)
    ax2.grid(True, alpha=0.3)
    # ax2.legend(fontsize=8)

    # Statistical Summary
    ax3.axis('off')
    if percentile_results:
        stats_text = "STATISTICAL SUMMARY\n"
        stats_text += "=" * 20 + "\n\n"

        for method_name, percentiles in percentile_results.items():
            stats_text += f"{method_name}:\n"
            stats_text += f"P50:  {percentiles.get('P50', 0):.1f}ms\n"
            stats_text += f"P90:  {percentiles.get('P90', 0):.1f}ms\n"
            stats_text += f"P99:  {percentiles.get('P99', 0):.1f}ms\n"
            stats_text += f"P99.9: {percentiles.get('P99.9', 0):.1f}ms\n"
            stats_text += "\n"

        ax3.text(0.1, 0.95, stats_text, transform=ax3.transAxes, fontsize=10,
                 verticalalignment='top', fontfamily='monospace',
                 bbox=dict(boxstyle="round,pad=0.5", facecolor="lightgray"))
        ax3.set_title('Statistical Summary', fontsize=12, fontweight='bold')

    plt.tight_layout()
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    plt.show()

    return percentile_results


# Usage example
if __name__ == "__main__":

    idx = '5'

    path1 = 'server/ginx_benchmark_results_' + idx + '.json'
    path2 = 'server/lazy_benchmark_results_' + idx + '.json'
    json_files = [
        path1,
        path2
    ]

    analyses = plot_corrected_benchmark_analysis(json_files, 'high_low_analysis_' + idx + '.png')

    if analyses:
        # Plot CDF with percentiles
        percentile_results = plot_cdf_with_percentiles(analyses, 'cdf_analysis_' + idx + '.png')

        # Plot detailed CDF analysis
        detailed_results = plot_detailed_cdf_analysis(analyses, 'detailed_cdf_analysis_' + idx + '.png')