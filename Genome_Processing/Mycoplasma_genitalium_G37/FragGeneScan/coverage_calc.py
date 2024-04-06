import numpy as np
from collections import defaultdict

# Function to calculate coverage metrics for each contig
def calculate_coverage_metrics(coverage_file):
    coverage_metrics = defaultdict(list)

    with open(coverage_file, 'r') as file:
        for line in file:
            contig, position, coverage = line.strip().split('\t')
            coverage = int(coverage)

            # Store coverage for each contig
            coverage_metrics[contig].append(coverage)

    # Calculate metrics for each contig
    for contig, coverages in coverage_metrics.items():
        avg_coverage = np.mean(coverages)
        min_coverage = min(coverages)
        max_coverage = max(coverages)
        median_coverage = np.median(coverages)

        # Update coverage metrics
        coverage_metrics[contig] = {
            'Average_Coverage': avg_coverage,
            'Min_Coverage': min_coverage,
            'Max_Coverage': max_coverage,
            'Median_Coverage': median_coverage
        }

    return coverage_metrics

# Function to write coverage metrics to a new file
def write_coverage_metrics(coverage_metrics, output_file):
    with open(output_file, 'w') as file:
        file.write("Contig\tAverage_Coverage\tMin_Coverage\tMax_Coverage\tMedian_Coverage\n")
        for contig, metrics in coverage_metrics.items():
            line = f"{contig}\t{metrics['Average_Coverage']:.2f}\t{metrics['Min_Coverage']}\t{metrics['Max_Coverage']}\t{metrics['Median_Coverage']}\n"
            file.write(line)

# Example usage
coverage_file = 'coverage.txt'
output_file = 'coverage_metrics.txt'

coverage_metrics = calculate_coverage_metrics(coverage_file)
write_coverage_metrics(coverage_metrics, output_file)
