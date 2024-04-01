from collections import Counter
import numpy as np
import matplotlib.pyplot as plt
import gzip
import scipy.stats as stats

# Function to parse Eggnog mapper output file and extract COG assignments
def parse_eggnog_output(file_path):
    cog_assignments = []
    #Lazy
    try:
        with open(file_path, 'r') as file:
            for line in file:
                if not line.startswith('#'):
                    cogs = line.split('\t')[6]
                    for cog in cogs:
                        if cog != 'Z' and cog != '-':
                            cog_assignments.append(cog)
    except:
        with gzip.open(file_path, 'rt') as file:
            for line in file:
                if not line.startswith('#'):
                    cogs = line.split('\t')[6]
                    for cog in cogs:
                        if cog != 'Z' and cog != '-':
                            cog_assignments.append(cog)
    return cog_assignments


# Function to ensure both datasets have the same COGs
def ensure_same_cogs(gene_level_counts, read_level_counts):
    all_cogs = set(gene_level_counts.keys()).union(set(read_level_counts.keys()))
    missing_cogs = []

    for cog in all_cogs:
        if cog not in gene_level_counts:
            missing_cogs.append(cog)
        if cog not in read_level_counts:
            missing_cogs.append(cog)

    for cog in missing_cogs:
        if cog in gene_level_counts:
            del gene_level_counts[cog]
        if cog in read_level_counts:
            del read_level_counts[cog]

    if missing_cogs:
        print(f"Missing COG categories: {missing_cogs}")

    return gene_level_counts, read_level_counts

# File paths for Eggnog mapper output files
# gene_level_file = '../Genome_Processing/Mycoplasma_genitalium_G37/Myco_pep_Emapper.emapper.annotations'
# read_level_files = ['../Genome_Processing/Mycoplasma_genitalium_G37/Processing/ART_Simulated_Reads/Myco_ART_errFree_Combined_Emapper.emapper.annotations',
#                     '../Genome_Processing/Mycoplasma_genitalium_G37/FragGeneScan/FragGeneScan_ART_errFree_Combined_Emapper.emapper.annotations',
#                     '../Genome_Processing/Mycoplasma_genitalium_G37/Pyrodigal/Pyrodigal_ART_errFree_Combined_Emapper.emapper.annotations',
#                     '../Genome_Processing/Staphylococcus_aureus_502A/FragGeneScan/FragGeneScan_ART_errFree_Combined_Emapper.emapper.annotations',
#                     '../Genome_Processing/Staphylococcus_aureus_502A/Pyrodigal/Pyrodigal_ART_errFree_Combined_Emapper.emapper.annotations',
#                     '../Genome_Processing/Staphylococcus_aureus_502A/Staph_pep_Emapper.emapper.annotations']
gene_level_file = '../Genome_Processing/Staphylococcus_aureus_502A/Staph_pep_Emapper.emapper.annotations'
read_level_files = [#'../Genome_Processing/Staphylococcus_aureus_502A/Processing/ART_Simulated_Reads/Staph_ART_errFree_Combined_Emapper.emapper.annotations.gz',
                    '../Genome_Processing/Staphylococcus_aureus_502A/FragGeneScan/FragGeneScan_ART_errFree_Combined_Emapper.emapper.annotations.gz',
                    '../Genome_Processing/Staphylococcus_aureus_502A/Pyrodigal/Pyrodigal_ART_errFree_Combined_Emapper.emapper.annotations.gz',
                    '../Genome_Processing/Mycoplasma_genitalium_G37/Myco_pep_Emapper.emapper.annotations']


# Parse gene-level Eggnog mapper output file
gene_level_cogs = parse_eggnog_output(gene_level_file)

# Count occurrences of COGs at gene level
gene_level_counts = Counter(gene_level_cogs)

# Individual figures for each comparison
for read_level_file in read_level_files:
    # Parse read-level Eggnog mapper output file
    read_level_cogs = parse_eggnog_output(read_level_file)

    # Count occurrences of COGs at read level
    read_level_counts = Counter(read_level_cogs)

    # Ensure both datasets have the same COGs
    gene_level_counts, read_level_counts = ensure_same_cogs(gene_level_counts, read_level_counts)

    # Compute total counts
    total_gene_level = sum(gene_level_counts.values())
    total_read_level = sum(read_level_counts.values())

    # Calculate observed COG proportions
    observed_proportions = {cog: count / total_read_level for cog, count in read_level_counts.items()}

    # Calculate expected COG proportions
    expected_proportions = {cog: count / total_gene_level for cog, count in gene_level_counts.items()}

    # Perform chi-squared test
    chi2, p = stats.chisquare(list(observed_proportions.values()), list(expected_proportions.values()))

    # Plot differences between observed and expected proportions
    cogs = sorted(gene_level_counts.keys())
    observed_values = [observed_proportions[cog] for cog in cogs]
    expected_values = [expected_proportions[cog] for cog in cogs]
    plt.figure(figsize=(10, 6))
    plt.bar(cogs, np.array(observed_values) - np.array(expected_values), color='blue')
    plt.xlabel('COG Category')
    plt.ylabel('Difference (Observed - Expected)')
    plt.title(f'Differences between Observed and Expected COG Proportions for {read_level_file}')
    plt.xticks(rotation=45, ha='right')
    plt.grid(axis='y')
    plt.tight_layout()
    plt.show()

    # Print results
    print("Read-level file:", read_level_file)
    print("Gene-level COG counts:", gene_level_counts)
    print("Read-level COG counts:", read_level_counts)
    print("Observed COG proportions:", {cog: observed_proportions[cog] for cog in cogs})
    print("Expected COG proportions:", {cog: expected_proportions[cog] for cog in cogs})
    print("Chi-squared test statistic:", chi2)
    print("P-value:", p)


# Custom legend names for each file
legend_names = {
    '../Genome_Processing/Staphylococcus_aureus_502A/FragGeneScan/FragGeneScan_ART_errFree_Combined_Emapper.emapper.annotations.gz': 'FragGeneScan',
    '../Genome_Processing/Staphylococcus_aureus_502A/Pyrodigal/Pyrodigal_ART_errFree_Combined_Emapper.emapper.annotations.gz': 'Pyrodigal',
    '../Genome_Processing/Mycoplasma_genitalium_G37/Myco_pep_Emapper.emapper.annotations': "Mycoplasma genitalium 'Full Proteins'"
}

# Combined "meta" figure
plt.figure(figsize=(12, 8))
for read_level_file in read_level_files:
    # Parse read-level Eggnog mapper output file
    read_level_cogs = parse_eggnog_output(read_level_file)

    # Count occurrences of COGs at read level
    read_level_counts = Counter(read_level_cogs)

    # Ensure both datasets have the same COGs
    gene_level_counts, read_level_counts = ensure_same_cogs(gene_level_counts, read_level_counts)

    # Compute total counts
    total_gene_level = sum(gene_level_counts.values())
    total_read_level = sum(read_level_counts.values())

    # Calculate observed COG proportions
    observed_proportions = {cog: count / total_read_level for cog, count in read_level_counts.items()}

    # Calculate expected COG proportions
    expected_proportions = {cog: count / total_gene_level for cog, count in gene_level_counts.items()}

    # Plot differences between observed and expected proportions
    cogs = sorted(gene_level_counts.keys())
    diffs = [observed_proportions[cog] - expected_proportions[cog] for cog in cogs]
    plt.plot(cogs, diffs, label=legend_names.get(read_level_file, read_level_file))  # Use custom legend name if available

plt.xlabel('COG Category', fontsize = 24)
plt.ylabel('Difference (Observed - Expected)', fontsize = 24)
plt.title('Differences between Observed and Expected COG Proportions', fontsize = 24)
plt.xticks(rotation=45, ha='right', fontsize = 14)
plt.yticks(fontsize = 14)
plt.legend(fontsize = 24)
plt.ylim(-0.20, 0.20)
plt.grid(axis='y')
plt.tight_layout()
plt.show()

