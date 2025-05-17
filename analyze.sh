#!/bin/bash

# Parallel processing with configurable settings
MAX_PARALLEL_JOBS=11  # Adjust based on your system's capabilities

# Function to process a single N and seed combination
process_analysis() {
    local N=$1
    local seed=$2
    
    local output_dir="analysis/N${N}_seed${seed}"
    local input_file="dumps/N${N}_seed${seed}/production.lammpstrj"
    local output_prefix="${output_dir}/N${N}_seed${seed}"
    
    # Create output directory if it doesn't exist
    mkdir -p "$output_dir"
    
    # Only run analysis if input file exists
    if [ -f "$input_file" ]; then
        echo "Processing N=${N}, seed=${seed}"
        ./analyze "$input_file" "$output_prefix"
    else
        echo "Warning: Input file not found for N=${N}, seed=${seed}"
    fi
}

# Export the function so it can be used with parallel
export -f process_analysis

# Create a temporary file to store all jobs
job_list=$(mktemp)

# Generate the job list
for N in 50 100 150 200 250 300 350 400 450 500; do
    for seed in {1..50}; do
        echo "$N $seed" >> "$job_list"
    done
done

# Run the jobs in parallel
if command -v parallel &> /dev/null; then
    # Use GNU parallel if available
    cat "$job_list" | parallel -j "$MAX_PARALLEL_JOBS" process_analysis {1} {2}
else
    # Fallback to xargs if parallel is not available
    cat "$job_list" | xargs -P "$MAX_PARALLEL_JOBS" -n 2 bash -c 'process_analysis "$0" "$1"'
fi

# Clean up
rm "$job_list"

echo "All analysis jobs completed"