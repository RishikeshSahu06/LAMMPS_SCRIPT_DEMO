import os
import subprocess
import multiprocessing
import time
import datetime
import csv
from contextlib import contextmanager

@contextmanager
def timing_context(description):
    start_time = time.time()
    yield
    elapsed_time = time.time() - start_time
    return elapsed_time

def run_simulation(params):
    N, damp, seed = params
    output_dir = f"dumps/N{N}_seed{seed}"
    os.makedirs(output_dir, exist_ok=True)

    with open("in.main", "r") as f:
        template = f.read()

    # Replace variables in the template
    input_script = (
        template
        .replace("${n}", str(N))
        .replace("${seed}", str(seed))
        .replace("${damp}", str(damp))
    )

    # Write customized input file
    input_path = f"{output_dir}/in.input"
    with open(input_path, "w") as f:
        f.write(input_script)

    # Run LAMMPS simulation with timing
    print(f"Starting simulation for N={N}, seed={seed}")
    start_time = time.time()
    result = subprocess.run(["lmp", "-in", input_path], capture_output=True, text=True)
    elapsed_time = time.time() - start_time
    print(f"Finished simulation for N={N}, seed={seed} in {elapsed_time:.2f} seconds")
    
    # Save individual timing data to a file in the output directory
    timing_file = f"{output_dir}/timing.txt"
    with open(timing_file, "w") as f:
        f.write(f"N={N}, seed={seed}, damp={damp}\n")
        f.write(f"Execution time: {elapsed_time:.2f} seconds\n")
        f.write(f"Return code: {result.returncode}\n")
        f.write(f"Timestamp: {datetime.datetime.now()}\n")
    
    return (result.returncode, N, seed, elapsed_time)

if __name__ == "__main__":
    overall_start_time = time.time()
    
    # Create directory for timing results
    os.makedirs("timing_results", exist_ok=True)
    
    # Define N and corresponding damp
    N_values = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500]
    damp_values = [100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0]
    seeds = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50]

    # Create all parameter combinations
    all_params = [(N, damp, seed) for N, damp in zip(N_values, damp_values) for seed in seeds]
    
    # Determine the number of CPU cores to use
    num_cores = max(1, multiprocessing.cpu_count() - 1)
    print(f"Running simulations using {num_cores} CPU cores")
    print(f"Total simulations to run: {len(all_params)}")
    
    # Create a pool of workers and distribute the tasks
    with multiprocessing.Pool(processes=num_cores) as pool:
        results = pool.map(run_simulation, all_params)
    
    # Calculate overall execution time
    overall_elapsed_time = time.time() - overall_start_time
    
    # Process and save results
    with open("timing_results/simulation_times.csv", "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["N", "Seed", "Execution Time (s)", "Status"])
        
        for result in results:
            returncode, N, seed, elapsed_time = result
            status = "Success" if returncode == 0 else "Failed"
            writer.writerow([N, seed, f"{elapsed_time:.2f}", status])
    
    # Calculate statistics by N value
    time_by_N = {}
    for result in results:
        returncode, N, seed, elapsed_time = result
        if returncode == 0:  # Only count successful runs
            if N not in time_by_N:
                time_by_N[N] = []
            time_by_N[N].append(elapsed_time)
    
    # Write summary statistics
    with open("timing_results/timing_summary.txt", "w") as f:
        f.write(f"Overall execution time: {overall_elapsed_time:.2f} seconds ({overall_elapsed_time/60:.2f} minutes)\n")
        f.write(f"Started at: {datetime.datetime.fromtimestamp(overall_start_time)}\n")
        f.write(f"Finished at: {datetime.datetime.now()}\n")
        f.write(f"Number of processes: {num_cores}\n")
        f.write(f"Total simulations: {len(all_params)}\n\n")
        
        f.write("Average execution time by N value:\n")
        for N in sorted(time_by_N.keys()):
            avg_time = sum(time_by_N[N]) / len(time_by_N[N])
            f.write(f"N={N}: {avg_time:.2f} seconds (avg of {len(time_by_N[N])} successful runs)\n")
    
    # Check if any simulations failed
    failed_simulations = [(N, seed) for (returncode, N, seed, _) in results if returncode != 0]
    if failed_simulations:
        print(f"Warning: {len(failed_simulations)} simulations failed:")
        with open("timing_results/failed_simulations.txt", "w") as f:
            for N, seed in failed_simulations:
                failure_msg = f"N={N}, seed={seed}"
                print(f"  {failure_msg}")
                f.write(f"{failure_msg}\n")
    else:
        print("All simulations completed successfully")
    
    print(f"Overall execution time: {overall_elapsed_time:.2f} seconds ({overall_elapsed_time/60:.2f} minutes)")
    print(f"Timing results saved to timing_results/")