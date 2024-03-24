from collections import Counter
import numpy as np
import scipy.stats as stats

# Function to parse Eggnog mapper output file and extract COG assignments
def parse_eggnog_output(file_path):
    cog_assignments = []
    with open(file_path, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                cogs = line.split('\t')[6]
                for cog in cogs:
                    if cog != 'Z':
                        cog_assignments.append(cog)
    return cog_assignments

# File paths for Eggnog mapper output files
gene_level_file = '/Users/nicholas/Documents/FragGeneScan_ART_errFree_Combined_Emapper.emapper.annotations'
read_level_files = ['/Users/nicholas/Documents/Pyrodigal_ART_errFree_Combined_Emapper.emapper.annotations']

# Parse gene-level Eggnog mapper output file
gene_level_cogs = parse_eggnog_output(gene_level_file)

# Parse read-level Eggnog mapper output files
read_level_cogs = []
for file in read_level_files:
    read_level_cogs.extend(parse_eggnog_output(file))

# Count occurrences of COGs
gene_level_counts = Counter(gene_level_cogs)
read_level_counts = Counter(read_level_cogs)

# Compute normalized counts
total_gene_level = sum(gene_level_counts.values())
total_read_level = sum(read_level_counts.values())
norm_read_level = {cog: count / total_read_level for cog, count in read_level_counts.items()}

# Calculate expected frequencies based on normalized read-level counts
expected = [total_gene_level * norm_read_level.get(cog, 0) for cog in gene_level_counts]

# Convert counts to proportions
observed_proportions = {cog: count / total_read_level for cog, count in read_level_counts.items()}
expected_proportions = {cog: count / total_gene_level for cog, count in zip(gene_level_counts.keys(), expected)}

# Perform chi-squared test
chi2, p = stats.chisquare(list(observed_proportions.values()), list(expected_proportions.values()))

# Print results
print("Gene-level COG counts:", gene_level_counts)
print("Read-level COG counts:", read_level_counts)
print("Normalized read-level COG proportions:", norm_read_level)
print("Observed COG proportions:", observed_proportions)
print("Expected COG proportions:", expected_proportions)
print("Chi-squared test statistic:", chi2)
print("P-value:", p)
