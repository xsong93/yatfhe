import json
import os

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats
from scipy.interpolate import griddata
from scipy.spatial import ConvexHull


def load_and_process_data(file_path):
    """Load and preprocess benchmark data"""
    with open(file_path) as f:
        data = json.load(f)

    benchmarks = []
    for run in data['benchmarks']:
        benchmarks.append({
            'batchSize': int(run['batchSize']),
            'tasksPerThread': int(run['tasksPerThread']),
            'real_time_ms': run['real_time'] / 1000,  # Convert to milliseconds
            'cpu_time_ms': run['cpu_time'] / 1000,
            'iterations': run['iterations'],
            'throughput': run['iterations'] / (run['real_time'] / 1e6)  # ops/sec
        })

    return pd.DataFrame(benchmarks)


def load_and_process_data2(file_path):
    """Load and preprocess benchmark data"""
    with open(file_path) as f:
        data = json.load(f)

    benchmarks = []
    for run in data['benchmarks']:
        benchmarks.append({
            'method': run['name'].split('/')[1],
            'real_time_ms': run['real_time'] / 1000,
        })

    return pd.DataFrame(benchmarks)


def create_heatmap(df, out_path: str):
    """Create heatmap"""
    plt.figure(figsize=(14, 10))

    # Create pivot table and convert to microseconds
    pivot = df.pivot_table(values='real_time_ms',
                           index='tasksPerThread',
                           columns='batchSize')

    # Custom colormap (yellow -> green)
    cmap = LinearSegmentedColormap.from_list('perf', ['#FFD700', '#2ECC71'])

    # Plot heatmap
    ax = sns.heatmap(pivot, annot=True, fmt=".1f", cmap=cmap,
                     cbar_kws={'label': 'Execution Time (ms)'},
                     annot_kws={'fontsize': 9})

    # Calculate performance metrics
    min_time = pivot.min().min()
    max_time = pivot.max().max()
    std_dev = pivot.values.std()

    # Highlight optimal regions
    optimal_mask = pivot <= (min_time + 0.2 * std_dev)
    for (j, i), val in np.ndenumerate(pivot):
        if optimal_mask.iloc[j, i]:
            ax.add_patch(plt.Rectangle((i, j), 1, 1, fill=False,
                                       edgecolor='red', lw=3))

    # Add analysis textbox
    analysis_text = f"""
    Performance Analysis:
    - Fastest config: {pivot.stack().idxmin()} @ {min_time:.1f}ms
    - Slowest config: {pivot.stack().idxmax()} @ {max_time:.1f}ms
    - Optimal zone (red): Within 20% of best time
    - Color range: {min_time:.1f}-{max_time:.1f}ms (Δ={max_time - min_time:.1f}ms)
    """
    plt.text(0.5, -0.3, analysis_text, transform=ax.transAxes,
             fontsize=11, bbox=dict(facecolor='white', alpha=0.8))

    # Formatting
    plt.title("Execution Time Heatmap with Performance Zones\n"
              "(Lower values are better)", pad=20)
    plt.xlabel("Batch Size", labelpad=10)
    plt.ylabel("Tasks Per Thread", labelpad=10)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()


def create_3d_surface(df, out_path: str):
    """Create 3D performance surface with dense plateau marking"""
    fig = plt.figure(figsize=(18, 14))
    ax = fig.add_subplot(111, projection='3d')

    # Prepare data
    X = df['batchSize'].values
    Y = df['tasksPerThread'].values
    Z = df['real_time_ms'].values

    z_min = Z.min() - 6
    z_max = Z.max() + 5

    # Create grid for surface plot
    xi = np.linspace(X.min(), X.max(), 100)
    yi = np.linspace(Y.min(), Y.max(), 100)
    zi = griddata((X, Y), Z, (xi[None, :], yi[:, None]), method='cubic')

    # Create surface
    surf = ax.plot_surface(xi[None, :], yi[:, None], zi,
                           cmap='berlin', edgecolor='none', alpha=0.7)

    # Calculate performance metrics
    min_idx = np.argmin(Z)
    max_idx = np.argmax(Z)
    min_point = (X[min_idx], Y[min_idx], Z[min_idx])
    max_point = (X[max_idx], Y[max_idx], Z[max_idx])
    mean_time = np.mean(Z)
    std_dev = np.std(Z)

    # Define plateau threshold
    plateau_threshold = mean_time - 0.5 * std_dev
    plateau_mask = Z <= plateau_threshold

    if plateau_mask.sum() > 1:
        plateau_points = np.column_stack([X[plateau_mask], Y[plateau_mask], Z[plateau_mask]])

        # 2. cal point density
        def calculate_local_density(points, radius=2.0):
            densities = []
            for i, point in enumerate(points):
                distances = np.sqrt(np.sum((points - point) ** 2, axis=1))
                density = np.sum(distances < radius)
                densities.append(density)
            return np.array(densities)

        densities = calculate_local_density(plateau_points[:, :2])  # XY plane only

        # 3. filter dense area
        density_threshold = np.percentile(densities, 75)
        dense_mask = densities >= density_threshold
        dense_points = plateau_points[dense_mask]

        # 4. marking & plotting
        if len(dense_points) > 0:
            ax.scatter(dense_points[:, 0], dense_points[:, 1], dense_points[:, 2],
                       color='gold', s=150, alpha=0.9, label='Dense Plateau Core')

            ax.scatter(plateau_points[~dense_mask, 0],
                       plateau_points[~dense_mask, 1],
                       plateau_points[~dense_mask, 2],
                       color='purple', s=80, alpha=0.6, label='Plateau Region')

            if len(dense_points) >= 3:
                try:
                    hull = ConvexHull(dense_points[:, :2])
                    boundary_label_added = False

                    for simplex in hull.simplices:
                        hull_points = dense_points[simplex, :2]
                        hull_z = np.full(2, Z.min() - 5)
                        label = 'Dense Region Boundary' if not boundary_label_added else ""
                        ax.plot(hull_points[:, 0], hull_points[:, 1], hull_z,
                                color='red', linewidth=3, linestyle='--',
                                label=label)
                        boundary_label_added = True  # 后续绘制不再添加标签
                except:
                    pass

        # 6. project 2D heatmap
        if len(dense_points) > 0:
            x_range = np.linspace(X.min(), X.max(), 50)
            y_range = np.linspace(Y.min(), Y.max(), 50)
            xx, yy = np.meshgrid(x_range, y_range)

            positions = np.vstack([xx.ravel(), yy.ravel()])
            kernel = stats.gaussian_kde(plateau_points[:, :2].T)
            density_grid = np.reshape(kernel(positions).T, xx.shape)

            ax.contourf(xx, yy, density_grid, zdir='z', offset=Z.min() - 10,
                        levels=10, alpha=0.3, cmap='hot')

    ax.scatter(X, Y, Z, color='black', s=30, alpha=0.4, label='Data Points')
    ax.scatter(*min_point, color='red', s=300, marker='*', label=f'Optimal: {Z[min_idx]:.1f}ms')

    dense_analysis = ""
    if plateau_mask.sum() > 1 and 'dense_points' in locals():
        dense_analysis = f"""
    Dense Plateau Analysis:
    - Plateau Points: {plateau_mask.sum()}/{len(Z)} total
    - Dense Core Size: {len(dense_points)} points
    - Density Threshold: {density_threshold:.1f} (top 25%)
    - Core Center: ({np.mean(dense_points[:, 0]):.1f}, {np.mean(dense_points[:, 1]):.1f})
    - Performance Range in Core: {dense_points[:, 2].min():.1f}-{dense_points[:, 2].max():.1f}ms
        """

    analysis_text = f"""
    Performance Characteristics:
    - Optimal Config: batch={min_point[0]:.0f}, tasks={min_point[1]:.0f}
    - Performance Range: {Z.min():.1f}ms to {Z.max():.1f}ms 
    - Plateau Threshold: {plateau_threshold:.1f}ms
    - Standard Deviation: {std_dev:.1f}ms
    {dense_analysis}
    """

    ax.text2D(0.05, 0.95, analysis_text, transform=ax.transAxes,
              bbox=dict(facecolor='white', alpha=0.9), fontsize=10)

    ax.set_zlim(z_min, z_max)
    ax.set_xlabel('Batch Size', labelpad=15, fontsize=12)
    ax.set_ylabel('Tasks Per Thread', labelpad=15, fontsize=12)
    ax.set_zlabel('Execution Time (ms)', labelpad=15, fontsize=12)
    ax.set_title('3D Performance Surface with Dense Plateau Analysis\n'
                 '(Lower is better)', pad=25, fontsize=20, loc="right")

    ax.legend(loc='upper right', bbox_to_anchor=(0.95, 0.85), fontsize=10)

    cbar = fig.colorbar(surf, shrink=0.5, aspect=5, pad=0.1)
    cbar.set_label('Execution Time (ms)', fontsize=12)

    ax.view_init(elev=25, azim=45)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()

    if plateau_mask.sum() > 1:
        return {
            'total_points': len(Z),
            'plateau_points': plateau_mask.sum(),
            'dense_core_size': len(dense_points),
            'density_threshold': density_threshold,
            'core_center': (np.mean(dense_points[:, 0]), np.mean(dense_points[:, 1])),
            'core_performance_range': (dense_points[:, 2].min(), dense_points[:, 2].max())
        }
    return None


def create_line_plot(df, out_path: str):
    """Create single line plot of time vs batch size"""
    plt.figure(figsize=(10, 6))

    # Plot settings
    markers = ['o', 's', '^', 'D', 'v', 'p', '*', 'X', '<', '>']
    colors = plt.cm.berlin(np.linspace(0, 1, len(df['tasksPerThread'].unique())))
    plt.style.use('seaborn-v0_8')

    # Plot each tasksPerThread as separate line
    for i, tpt in enumerate(sorted(df['tasksPerThread'].unique())):
        subset = df[df['tasksPerThread'] == tpt]
        plt.plot(subset['batchSize'],
                 subset['real_time_ms'],
                 marker=markers[i % len(markers)],
                 linestyle='-',
                 linewidth=2,
                 markersize=6,
                 color=colors[i],
                 label=f'Tasks/Thread={tpt}')

    # Customize plot
    plt.xlabel('Batch Size', fontsize=12)
    plt.ylabel('Execution Time (ms)', fontsize=12)
    plt.title('Execution Time by Batch Size', fontsize=14, pad=20)
    plt.grid(True, alpha=0.3)

    # 2-column legend
    plt.legend(ncol=2,
               loc='upper right',
               frameon=True,
               framealpha=0.9,
               fancybox=True,
               shadow=True,
               fontsize=9)

    # Save and close
    plt.tight_layout()
    plt.subplots_adjust(right=0.85)
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()


def plot_benchmark_barchart(df, out_path: str, method_name: str):
    # Calculate performance metrics
    fastest = df.loc[df['real_time_ms'].idxmin()]
    slowest = df.loc[df['real_time_ms'].idxmax()]

    plt.figure(figsize=(14, 8))
    ax = sns.barplot(data=df.sort_values('real_time_ms'),
                     x='method',
                     y='real_time_ms',
                     hue="method", palette="dark", legend=False)

    # Highlight fastest method
    fastest_idx = df['real_time_ms'].idxmin()
    ax.patches[fastest_idx].set_edgecolor('red')
    ax.patches[fastest_idx].set_linewidth(3)

    # Add performance ratio annotations
    for i, (_, row) in enumerate(df.sort_values('real_time_ms').iterrows()):
        ratio = (slowest['real_time_ms'] / row['real_time_ms'])
        ax.text(i, row['real_time_ms'] / 2, f'{ratio:.1f}x',
                ha='center', va='center', color='white', fontweight='bold')

    # Add value labels on top of bars
    for p in ax.patches:
        ax.text(p.get_x() + p.get_width() / 2.,  # x-position
                p.get_height() + 0.02 * max(df['real_time_ms']),  # y-position (slightly above bar)
                f'{p.get_height():.1f} ms',  # text
                ha='center', va='bottom',  # alignment
                fontsize=10)

    # Formatting
    plt.title(method_name + ' Performance Comparison\n'
              'Benchmark Iteration = 500',
              pad=20)
    plt.xlabel('Method', labelpad=15)
    plt.ylabel('Execution Time (ms)', labelpad=15)
    plt.xticks(rotation=45, ha='right')

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()


def generate_report(df, out_path: str):
    """Generate comprehensive analysis report with UTF-8 encoding"""
    report = []

    # Find optimal configurations
    min_time = df.loc[df['real_time_ms'].idxmin()]

    report.append("=== Performance Analysis Report ===")
    report.append(f"\nFastest Configuration:")
    report.append(f"- Batch Size: {min_time['batchSize']}")
    report.append(f"- Tasks/Thread: {min_time['tasksPerThread']}")
    report.append(f"- Time: {min_time['real_time_ms']:.1f} ms")

    # Scaling analysis
    report.append("\n=== Scaling Analysis ===")
    for bs in sorted(df['batchSize'].unique()):
        subset = df[df['batchSize'] == bs]
        min_time = subset['real_time_ms'].min()
        report.append(f"Batch Size {bs}: Min time = {min_time:.1f} ms")

    with open(out_path, 'w', encoding='utf-8') as f:
        f.write("\n".join(report))


if __name__ == "__main__":
    # out_dir = 'results_i5'
    # data1 = 'result_batchsize_i5.json'
    # data2 = 'result_i5_bench_blindrotate_all.json'

    out_dir = 'results_9950x3d'
    data2 = 'result_9950x3d_bench_blindrotate_basemethod.json'
    data3 = 'result_9950x3d_bench_keygen.json'
    data4 = 'result_9950x3d_bench_readkey.json'
    data5 = 'result_9950x3d_bench_writekey.json'

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    # Load and process data
    # df = load_and_process_data(data1)
    df2 = load_and_process_data2(data2)
    df3 = load_and_process_data2(data3)
    df4 = load_and_process_data2(data4)
    df5 = load_and_process_data2(data5)

    # Generate visualizations
    # create_heatmap(df, out_dir+'/heatmap.png')
    # create_3d_surface(df, out_dir+'/3d_surface_dense_plateau.png')
    # create_line_plot(df, out_dir+'/execution_time.png')
    plot_benchmark_barchart(df2, out_dir+'/comparison_barchart2.png', "Blind Rotation")

    plot_benchmark_barchart(df3, out_dir+'/comparison_barchart3.png', "Bootstrapping Key Gen")
    plot_benchmark_barchart(df4, out_dir+'/comparison_barchart4.png', "Initial Run")
    plot_benchmark_barchart(df5, out_dir+'/comparison_barchart5.png', "Write Key")


    # Generate report
    # generate_report(df, out_dir+'/performance_report.txt')

    print("Analysis complete. Created:")
    print("- heatmap.png")
    print("- 3d_surface.png")
    print("- line_plot.png")
    print("- bar_chart.png")
    print("- performance_report.txt")
