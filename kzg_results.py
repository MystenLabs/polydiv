import json
import os
import pandas as pd
import matplotlib.pyplot as plt

# Directory containing the Criterion JSON files
criterion_dir = 'target/criterion'

# Function to recursively find all relevant JSON files in a directory
def find_benchmark_files(directory):
    benchmark_files = []
    for root, _, files in os.walk(directory):
        if 'benchmark.json' in files and 'estimates.json' in files:
            benchmark_files.append({
                'benchmark': os.path.join(root, 'benchmark.json'),
                'estimates': os.path.join(root, 'estimates.json')
            })
    return benchmark_files

# Collect all benchmark file sets
benchmark_files = find_benchmark_files(criterion_dir)

# Dictionary to hold benchmark results grouped by function
results_by_function = {}

# Read and parse each set of benchmark files
for files in benchmark_files:
    with open(files['benchmark']) as f:
        benchmark_data = json.load(f)
    with open(files['estimates']) as f:
        estimates_data = json.load(f)
    
    # Extract relevant information
    function_id = benchmark_data.get('function_id', '')
    function_parts = function_id.split('/')
    if len(function_parts) >= 3:
        kzg_variant = function_parts[0]
        function_name = function_parts[1]
        size = int(function_parts[2])  # Convert size to integer
        short_function_id = f'{kzg_variant}/{function_name}/{size}'
    else:
        short_function_id = function_id

    # Skip KZGFK functions
    if kzg_variant == 'KZGFK':
        continue
    
    result = {
        'kzg_variant': kzg_variant,
        'function_name': function_name,
        'size': size,
        'mean': estimates_data['mean']['point_estimate']
    }
    
    if function_name not in results_by_function:
        results_by_function[function_name] = []
    results_by_function[function_name].append(result)

# Write results to separate CSV files and create plots
for function_name, results in results_by_function.items():
    # Convert results to DataFrame for easier sorting and manipulation
    df = pd.DataFrame(results)
    
    # Sort the DataFrame by size
    df.sort_values(by='size', inplace=True)
    
    # Write the sorted DataFrame to a CSV file
    csv_filename = f'{function_name}_results.csv'
    df.to_csv(csv_filename, index=False)
    
    # Plot the benchmark results
    plt.figure(figsize=(12, 8))

    colors = {'KZG': 'r', 'KZGDeriv': 'b', 'KZGTabDFK': 'm'}
    for kzg_variant in df['kzg_variant'].unique():
        subset = df[df['kzg_variant'] == kzg_variant]
        plt.plot(subset['size'], subset['mean'], label=kzg_variant, color=colors.get(kzg_variant, 'k'), marker='o', linewidth=4, markersize=6)

    plt.xlabel('Size (n)', fontsize=16, fontweight='bold')
    plt.ylabel('Mean Time in Nanosecends', fontsize=16, fontweight='bold')
    plt.yscale('log')  # Use logarithmic scale for the y-axis
    plt.title(f'{function_name.capitalize()} Benchmark Results', fontsize=16, fontweight='bold')
    plt.legend(fontsize=16)
    plt.grid(False)
    plt.xticks([16, 512, 1024, 2048, 4096, 8192], fontsize=12, fontweight='bold')
    plt.yticks(fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'{function_name}_results.png')
    plt.close()
