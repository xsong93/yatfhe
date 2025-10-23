import subprocess
import json
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

class BenchmarkOptimizer:
    def __init__(self, binary_path):
        self.binary_path = binary_path
        self.history = []
        self.best_params = None
        self.best_time = float('inf')  # Initialize with worst possible time
        self.results_cache = {}  # Cache to avoid redundant benchmarks

    def run_benchmark(self, batch_size, tasks_per_thread):
        """Run the benchmark binary with given parameters and return real time"""
        cache_key = (batch_size, tasks_per_thread)
        if cache_key in self.results_cache:
            return self.results_cache[cache_key]

        cmd = [
            self.binary_path,
            "--benchmark_filter=LAZY_MULTITHREAD",
            f"--batchSize={batch_size}",
            f"--tasksPerThread={tasks_per_thread}",
            "--benchmark_format=json"
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        try:
            data = json.loads(result.stdout)
            # Extract the real time (shorter is better)
            for benchmark in data["benchmarks"]:
                if "LAZY_MULTITHREAD" in benchmark["name"]:
                    real_time = benchmark["real_time"]
                    self.results_cache[cache_key] = real_time
                    return real_time
        except json.JSONDecodeError:
            print(f"Error parsing benchmark output for {batch_size},{tasks_per_thread}")
            return float('inf')  # Return worst possible time on error

        return float('inf')

    def optimize_parameter(self, fixed_param, fixed_value, param_to_optimize, range_values):
        """Optimize one parameter while keeping another fixed"""
        results = []
        for value in range_values:
            if param_to_optimize == "batch_size":
                time = self.run_benchmark(value, fixed_value)
            else:
                time = self.run_benchmark(fixed_value, value)
            results.append((value, time))

        # Find parameter with minimum time
        best_value, best_time = min(results, key=lambda x: x[1])
        return best_value, best_time, results

    def coordinate_descent(self, initial_batch_size, initial_tasks_per_thread, max_iter=5):
        """Perform coordinate descent optimization"""
        batch_size = initial_batch_size
        tasks_per_thread = initial_tasks_per_thread
        best_time = float('inf')
        best_params = (batch_size, tasks_per_thread)

        for iteration in range(max_iter):
            print(f"\n=== Iteration {iteration + 1} ===")

            # Optimize tasks_per_thread with fixed batch_size
            print(f"Optimizing tasks_per_thread with batch_size={batch_size}")
            best_tpt, tpt_time, tpt_results = self.optimize_parameter(
                "batch_size", batch_size, "tasks_per_thread", range(1, 21)
            )
            print(f"Best tasks_per_thread={best_tpt} with {tpt_time:.2f} μs")

            # Optimize batch_size with fixed tasks_per_thread
            print(f"Optimizing batch_size with tasks_per_thread={best_tpt}")
            best_bs, bs_time, bs_results = self.optimize_parameter(
                "tasks_per_thread", best_tpt, "batch_size", range(1, 21)
            )
            print(f"Best batch_size={best_bs} with {bs_time:.2f} μs")

            # Check for improvement
            current_time = min(tpt_time, bs_time)
            current_params = (best_bs, best_tpt)

            if current_time < best_time:
                best_time = current_time
                best_params = current_params
                batch_size, tasks_per_thread = best_params
                print(f"New best time: {best_time:.2f} μs")
            else:
                print("No improvement - stopping optimization")
                break

            self.history.append({
                'iteration': iteration + 1,
                'best_params': best_params,
                'best_time': best_time,
                'tasks_results': tpt_results,
                'batch_results': bs_results
            })

        return best_params, best_time

    def local_grid_search(self, center_bs, center_tpt, radius=2):
        """Perform a local grid search around found optimum"""
        bs_range = range(max(1, center_bs - radius), min(21, center_bs + radius + 1))
        tpt_range = range(max(1, center_tpt - radius), min(21, center_tpt + radius + 1))

        best_time = float('inf')
        best_params = (center_bs, center_tpt)

        print(f"\n=== Local Grid Search around {center_bs},{center_tpt} ===")
        for bs in bs_range:
            for tpt in tpt_range:
                time = self.run_benchmark(bs, tpt)
                print(f"  batch_size={bs}, tasks_per_thread={tpt}: {time:.2f} μs")
                if time < best_time:
                    best_time = time
                    best_params = (bs, tpt)

        return best_params, best_time

    def visualize_results(self):
        """Visualize the optimization results"""
        if not self.history:
            print("No results to visualize")
            return

        # Create heatmap of all tested combinations
        all_bs = set()
        all_tpt = set()
        time_map = defaultdict(dict)

        for cache_key, time in self.results_cache.items():
            bs, tpt = cache_key
            all_bs.add(bs)
            all_tpt.add(tpt)
            time_map[bs][tpt] = time

        # Convert to arrays for plotting
        bs_list = sorted(all_bs)
        tpt_list = sorted(all_tpt)
        time_array = np.zeros((len(bs_list), len(tpt_list)))

        for i, bs in enumerate(bs_list):
            for j, tpt in enumerate(tpt_list):
                time_array[i,j] = time_map.get(bs, {}).get(tpt, float('inf'))

        # Create heatmap
        plt.figure(figsize=(10, 8))
        plt.imshow(time_array, cmap='viridis_r', origin='lower')  # Reversed colormap for time (shorter=better)
        plt.colorbar(label='Real Time (μs)')
        plt.xticks(np.arange(len(tpt_list)), tpt_list)
        plt.yticks(np.arange(len(bs_list)), bs_list)
        plt.xlabel('Tasks per Thread')
        plt.ylabel('Batch Size')
        plt.title('Parameter Optimization Heatmap (Lower is Better)')

        # Mark the best point
        best_bs, best_tpt = self.best_params
        bs_idx = bs_list.index(best_bs)
        tpt_idx = tpt_list.index(best_tpt)
        plt.scatter(tpt_idx, bs_idx, color='red', marker='x', s=100)

        plt.show()

    def optimize(self):
        """Main optimization routine with multiple starting points"""
        starting_points = [
            (10, 10), # Middle ground
        ]

        overall_best_params = None
        overall_best_time = float('inf')

        for start_bs, start_tpt in starting_points:
            print(f"\n=== Starting optimization from {start_bs},{start_tpt} ===")
            best_params, best_time = self.coordinate_descent(start_bs, start_tpt)

            # Perform local grid search around found optimum
            best_params, best_time = self.local_grid_search(*best_params)

            if best_time < overall_best_time:
                overall_best_time = best_time
                overall_best_params = best_params

        self.best_params = overall_best_params
        self.best_time = overall_best_time

        print("\n=== Final Results ===")
        print(f"Best parameters: batch_size={self.best_params[0]}, tasks_per_thread={self.best_params[1]}")
        print(f"Best real time: {self.best_time:.2f} μs")

        # Visualize results
        self.visualize_results()

        return self.best_params, self.best_time

if __name__ == "__main__":
    # Replace with your actual benchmark binary path
    optimizer = BenchmarkOptimizer("./bench_batchsize_lazy")
    best_params, best_time = optimizer.optimize()