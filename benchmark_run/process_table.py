import matplotlib.pyplot as plt
import json
import os
import numpy as np


def calculate_percentiles(data):
    """Calculate P50, P90, P99, P99.9 percentiles"""
    if not data:
        return {'P50': 0, 'P90': 0, 'P99': 0, 'P99.9': 0}

    return {
        'P50': np.percentile(data, 50),
        'P90': np.percentile(data, 90),
        'P99': np.percentile(data, 99),
        'P99.9': np.percentile(data, 99.9)
    }


def analyze_data(data, method_name):
    """Analyze benchmark data and calculate comprehensive statistics"""
    times_ms = [item['time_us'] / 1000 for item in data['iterations']]
    data_mean = np.mean(times_ms) if times_ms else 0
    data_std = np.std(times_ms) if times_ms else 0

    # Calculate percentiles
    percentiles = calculate_percentiles(times_ms)

    return {
        'method': method_name,
        'all_times': times_ms,
        'avg': data_mean,
        'std': data_std,
        'p50': percentiles['P50'],
        'p99': percentiles['P99'],
        'p99_9': percentiles['P99.9'],
        'pressure_ratio': data.get('pressure_ratio', 1),
        'hit_rate': data['statistics']['hit_rate']
    }


def create_performance_table(json_file_paths, output_filename="performance_table.pdf"):
    """Create a PDF with a table showing performance metrics similar to Excel"""

    # Check if files exist
    existing_files = [f for f in json_file_paths if os.path.exists(f)]
    if not existing_files:
        print("No valid files found")
        return None, None

    # Analyze each dataset
    analyses = []
    for file_path in existing_files:
        try:
            with open(file_path, 'r', encoding='utf-8') as file:
                data = json.load(file)

            # Get method name
            method_name = data.get('benchmark_name', 'Unknown').split('/')[-1]

            # Use analysis function
            analysis = analyze_data(data, method_name)
            analyses.append(analysis)

        except Exception as e:
            print(f"Error processing {file_path}: {e}")
            continue

    if not analyses:
        print("No valid analyses to plot")
        return None, None

    # Group data by pressure ratio and method
    data_by_pressure_ratio = {}
    for analysis in analyses:
        pressure_ratio = analysis['pressure_ratio']
        method = analysis['method']

        if pressure_ratio not in data_by_pressure_ratio:
            data_by_pressure_ratio[pressure_ratio] = {}

        data_by_pressure_ratio[pressure_ratio][method] = analysis

    # Get all pressure ratios and methods
    pressure_ratios = sorted(data_by_pressure_ratio.keys())

    # Determine all methods
    all_methods = set()
    for ratio_data in data_by_pressure_ratio.values():
        all_methods.update(ratio_data.keys())
    all_methods = sorted(list(all_methods))

    # Create figure for the table
    # Adjust figure size based on number of rows (more compact)
    n_rows = len(pressure_ratios) * len(all_methods) + 1  # +1 for header
    # Reduced row height and figure margins
    fig_height = max(3, 0.3 * n_rows)  # Reduced minimum height and row height
    fig_width = 10  # Slightly reduced width for compactness

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    ax.axis('off')

    # Create table data
    table_data = []

    # Define column headers (matching the Excel table)
    headers = [
        'Pressure Ratio',
        'Method',
        'Cache HR(%)',
        'Avg(ms)',
        'Std(ms)',
        'P50(ms)',
        'P99(ms)',
        'P99.9(ms)'
    ]

    # Add rows for each pressure ratio and method
    for pressure_ratio in pressure_ratios:
        for method in all_methods:
            if method in data_by_pressure_ratio[pressure_ratio]:
                data = data_by_pressure_ratio[pressure_ratio][method]
                std_ratio = data['std'] / data_by_pressure_ratio[pressure_ratio]['OURS']['std']
                row = [
                    f"{pressure_ratio}",
                    method,
                    f"{data.get('hit_rate', 0):.2f}",
                    f"{data['avg']:.1f}",
                    f"{data['std']:.1f} ({std_ratio:.1f}X)",
                    f"{data['p50']:.1f}",
                    f"{data['p99']:.1f}",
                    f"{data['p99_9']:.1f}"
                ]
            else:
                # If data is missing for this pressure ratio and method, fill with empty values
                row = [f"{pressure_ratio}", method, "-", "-", "-", "-", "-", "-"]

            table_data.append(row)

    # Create table
    table = ax.table(
        cellText=table_data,
        colLabels=headers,
        loc='center',
        cellLoc='center',
        colWidths=[0.12, 0.1, 0.12, 0.1, 0.1, 0.1, 0.1, 0.1]  # Adjusted column widths
    )

    # Style the table
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.2)

    # Style header row
    for i in range(len(headers)):
        cell = table[(0, i)]
        cell.set_text_props(weight='bold', color='black', fontsize=9)  # Slightly smaller header font
        cell.set_height(0.15)  # Set header row height

    # Style data rows
    for i in range(len(table_data)):
        row_idx = i + 1
        for j in range(len(headers)):
            cell = table[(row_idx, j)]
            cell.set_height(0.12)

            if table_data[i][1] == 'TFHE':
                cell.set_facecolor('#E6F3FF')  # Light blue for TFHE
            elif table_data[i][1] == 'WWL+24':
                cell.set_facecolor('#FFF2E6')  # Light orange for WWL+24
            elif table_data[i][1] == 'OURS':
                cell.set_facecolor('#98FB98')  # Light green for OURS
            else:
                # Alternate row colors for other methods
                row_color = '#FFFFFF' if i % 2 == 0 else '#F2F2F2'
                cell.set_facecolor(row_color)

    # Highlight important columns
    highlight_cols = ['Std(ms)']
    col_indices = [headers.index(col) for col in highlight_cols if col in headers]
    for i in range(len(table_data)):
        for col_idx in col_indices:
            cell = table[(i+1, col_idx)]
            cell.set_text_props(weight='bold')

    # Remove all borders to make it cleaner
    for key, cell in table.get_celld().items():
        cell.set_edgecolor('black')
        cell.set_linewidth(0.5)

    # Adjust layout - REMOVE ALL PADDING
    plt.tight_layout(pad=0)  # Remove all padding

    # Set all margins to zero
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0, wspace=0, hspace=0)

    # Save as PDF (vector format) with tight bounding box and no padding
    plt.savefig(output_filename, format='pdf', dpi=300,
                bbox_inches='tight', pad_inches=0)
    print(f"Table saved as: {output_filename}")

    # Also save as PNG for quick viewing
    png_filename = output_filename.replace('.pdf', '.png')
    plt.savefig(png_filename, dpi=300,
                bbox_inches='tight', pad_inches=0)
    print(f"Table also saved as: {png_filename}")

    plt.show()

    return analyses, data_by_pressure_ratio


# Usage example
if __name__ == "__main__":

    pressure_ratios = [1, 5, 10, 20, 50]
    methods = ['TFHE', 'WWL+24', 'OURS']

    json_files = []
    for ratio in pressure_ratios:
        json_files.append(f'server/tfhe_benchmark_results_{ratio}.json')
        json_files.append(f'server/wwl+24_benchmark_results_{ratio}.json')
        json_files.append(f'server/ours_benchmark_results_{ratio}.json')

    existing_files = [f for f in json_files if os.path.exists(f)]

    if existing_files:
        print(f"Found {len(existing_files)} data files")

        analyses, data_by_pressure_ratio = create_performance_table(
            existing_files,
            "fig_performance_summary.pdf"
        )
    else:
        print("No data files found.")