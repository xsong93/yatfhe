import os
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

def get_file_size(file_path, unit='MB'):
    """
    Get file size and convert to specified unit

    Parameters:
    -----------
    file_path : str
        File path
    unit : str
        Unit: 'B', 'KB', 'MB', 'GB'

    Returns:
    --------
    float
        File size in specified unit
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")

    # Get file size in bytes
    size_bytes = os.path.getsize(file_path)

    # Convert to specified unit
    units = {
        'B': 1,
        'KB': 1024,
        'MB': 1024**2,
        'GB': 1024**3
    }

    if unit not in units:
        raise ValueError(f"Unsupported unit: {unit}, please use 'B', 'KB', 'MB', or 'GB'")

    return size_bytes / units[unit]


def get_best_unit(size_bytes):
    """
    Automatically select the best unit based on file size

    Parameters:
    -----------
    size_bytes : int
        File size in bytes

    Returns:
    --------
    tuple
        (size, unit)
    """
    if size_bytes < 1024:
        return size_bytes, 'B'
    elif size_bytes < 1024**2:
        return size_bytes / 1024, 'KB'
    elif size_bytes < 1024**3:
        return size_bytes / (1024**2), 'MB'
    else:
        return size_bytes / (1024**3), 'GB'


def plot_file_size_comparison(file_paths, output_filename):
    """
    Plot bar chart comparing file sizes, with multiplier labels relative to 'WWL+24'.

    Parameters:
    -----------
    file_paths : list
        List of file paths
    output_filename : str
        Output filename
    """
    if len(file_paths) < 2:
        raise ValueError("At least two file paths are required")

    # Check if files exist
    missing_files = [f for f in file_paths if not os.path.exists(f)]
    if missing_files:
        raise FileNotFoundError(f"The following files do not exist: {missing_files}")

    # Get file names and sizes
    file_names = []
    file_sizes_bytes = []
    file_sizes_display = []

    for file_path in file_paths:
        # Get file name
        file_name = os.path.basename(file_path).split('_')[1]
        if file_name == 'LAZY':
            file_name = 'OURS'
        if file_name == 'GINX':
            file_name = 'TFHE'
        file_names.append(file_name)

        # Get file size in bytes
        size_bytes = os.path.getsize(file_path)
        file_sizes_bytes.append(size_bytes)

        # Convert to best unit for display
        size_value, unit = get_best_unit(size_bytes)
        file_sizes_display.append(f"{size_value:.2f} {unit}")

    # --- New: calculate multiplier relative to baseline (WWL+24) ---
    base_label = 'WWL+24'
    try:
        base_idx = file_names.index(base_label)
    except ValueError:
        print(f"Warning: Baseline label '{base_label}' not found. Using the first file as baseline (multipliers may be meaningless).")
        base_idx = 0

    base_size = file_sizes_bytes[base_idx]
    multipliers = [size / base_size for size in file_sizes_bytes]
    # ----------------------------------------------------------------

    # Create figure
    fig, ax = plt.subplots(figsize=(8, 6))

    # Set colors
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

    # Plot bar chart
    bars = ax.bar(file_names, file_sizes_bytes,
                  color=colors[:len(file_names)],
                  edgecolor='black', linewidth=1.5,
                  alpha=0.8)

    # Add value labels on bars (size and multiplier)
    for bar, size_display, mult in zip(bars, file_sizes_display, multipliers):
        height = bar.get_height()
        # Show file size (original label)
        ax.text(bar.get_x() + bar.get_width()/2., height + 14.17,
                size_display, ha='center', va='bottom',
                fontsize=10, fontweight='bold')
        # Show multiplier (new, gray small text)
        ax.text(bar.get_x() + bar.get_width()/2., height*0.5,
                f"({mult:.2f}x)", ha='center', va='bottom',
                fontsize=9, color='black', fontweight='bold')

    # Adjust y-axis upper limit to reserve space for multiplier labels
    max_height = max(file_sizes_bytes)
    ax.set_ylim(0, max_height * 1.15)

    # Set title and labels
    ax.set_title('Key Size Comparison', fontsize=16, fontweight='bold', pad=20)
    ax.set_xlabel('Method', fontsize=12)
    ax.set_ylabel('Key Size (Bytes)', fontsize=12)

    # Set y-axis format
    def format_bytes(x, pos):
        if x >= 1e9:  # GB
            return f'{x/1e9:.2f}G'
        elif x >= 1e6:  # MB
            return f'{x/1e6:.2f}M'
        elif x >= 1e3:  # KB
            return f'{x/1e3:.2f}K'
        else:
            return f'{x:.0f}'

    ax.yaxis.set_major_formatter(ticker.FuncFormatter(format_bytes))

    # Add grid
    ax.grid(True, axis='y', alpha=0.3, linestyle='--')

    # Add border
    for spine in ax.spines.values():
        spine.set_edgecolor('gray')
        spine.set_linewidth(1)

    plt.tight_layout()

    # Save figure
    plt.savefig(output_filename + '.png', dpi=300, bbox_inches='tight')
    print(f"Chart saved as: {output_filename}")

    # plt.show()

    plt.savefig(output_filename + '.pdf', format='pdf', dpi=300,
                bbox_inches='tight', pad_inches=0)

    return file_names, file_sizes_bytes, file_sizes_display


# Usage example
if __name__ == "__main__":
    file1 = "server/BSK_GINX_1.bin"
    file3 = "server/BSK_LAZY_1.bin"
    file2 = "server/BSK_WWL+24_1.bin"

    file_paths = [file1, file2, file3]

    try:
        # Create basic comparison chart
        print("="*50)
        print("Basic File Size Comparison Chart")
        print("="*50)

        names, sizes_bytes, sizes_display = plot_file_size_comparison(
            file_paths,
            "fig_key_size_comparison"
        )

    except FileNotFoundError as e:
        print(f"Error: {e}")