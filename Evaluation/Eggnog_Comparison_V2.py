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




# File paths for Eggnog mapper output files
# gene_level_file = '../Genome_Processing/Mycoplasma_genitalium_G37/Myco_pep_Emapper.emapper.annotations'
# read_level_files = ['../Genome_Processing/Mycoplasma_genitalium_G37/Processing/ART_Simulated_Reads/Myco_ART_errFree_Combined_Emapper.emapper.annotations',
#                     '../Genome_Processing/Mycoplasma_genitalium_G37/FragGeneScan/FragGeneScan_ART_errFree_Combined_Emapper.emapper.annotations',
#                     '../Genome_Processing/Mycoplasma_genitalium_G37/Pyrodigal/Pyrodigal_ART_errFree_Combined_Emapper.emapper.annotations',
#                     '../Genome_Processing/Staphylococcus_aureus_502A/FragGeneScan/FragGeneScan_ART_errFree_Combined_Emapper.emapper.annotations',
#                     '../Genome_Processing/Staphylococcus_aureus_502A/Pyrodigal/Pyrodigal_ART_errFree_Combined_Emapper.emapper.annotations',
#                     '../Genome_Processing/Staphylococcus_aureus_502A/Staph_pep_Emapper.emapper.annotations']
gene_level_file = '../Genome_Processing/Mycoplasma_genitalium_G37/Myco_pep_Emapper.emapper.annotations'
read_level_files = [#'../Genome_Processing/Staphylococcus_aureus_502A/Processing/ART_Simulated_Reads/Staph_ART_errFree_Combined_Emapper.emapper.annotations.gz',
                    '../Genome_Processing/Mycoplasma_genitalium_G37/FragGeneScan/FragGeneScan_ART_errFree_Combined_Emapper.emapper.annotations',
                    '../Genome_Processing/Mycoplasma_genitalium_G37/Pyrodigal/Pyrodigal_ART_errFree_Combined_Emapper.emapper.annotations',
                    '../Genome_Processing/Staphylococcus_aureus_502A/Staph_pep_Emapper.emapper.annotations']


# Parse gene-level Eggnog mapper output file
gene_level_cogs = parse_eggnog_output(gene_level_file)

# Count occurrences of COGs at gene level
gene_level_counts = Counter(gene_level_cogs)




# Custom legend names for each file
legend_names = {
    '../Genome_Processing/Mycoplasma_genitalium_G37/FragGeneScan/FragGeneScan_ART_errFree_Combined_Emapper.emapper.annotations': 'M. genitalium - FragGeneScan',
    '../Genome_Processing/Mycoplasma_genitalium_G37/Pyrodigal/Pyrodigal_ART_errFree_Combined_Emapper.emapper.annotations': 'M. genitalium - Pyrodigal',
    '../Genome_Processing/Staphylococcus_aureus_502A/Staph_pep_Emapper.emapper.annotations': 'S. aureus - Full Proteins'
}

# Combined "meta" figure
plt.figure(figsize=(12, 8))

# Collect all unique COGs from both gene-level and read-level datasets
all_cogs = set(gene_level_counts.keys())
for read_level_file in read_level_files:
    read_level_cogs = parse_eggnog_output(read_level_file)
    all_cogs.update(read_level_cogs)

# Initialize dictionaries to store counts for each COG
gene_level_counts_all = {cog: gene_level_counts.get(cog, 0) for cog in all_cogs}
observed_proportions_all = {}
expected_proportions_all = {}

# Compute total counts at gene level
total_gene_level = sum(gene_level_counts.values())

# Calculate observed and expected COG proportions for each COG
for cog in all_cogs:
    # Observed proportion at gene level
    observed_proportions_all[cog] = 0
    if cog in gene_level_counts:
        observed_proportions_all[cog] += gene_level_counts[cog] / total_gene_level

    # Initialize expected proportion at gene level
    expected_proportions_all[cog] = 0

# Plot differences between observed and expected proportions for each read-level file
for read_level_file in read_level_files:
    # Parse read-level Eggnog mapper output file
    read_level_cogs = parse_eggnog_output(read_level_file)

    # Count occurrences of COGs at read level
    read_level_counts = Counter(read_level_cogs)

    # Ensure both datasets have the same COGs
    # gene_level_counts, read_level_counts = ensure_same_cogs(gene_level_counts, read_level_counts)

    # Compute total counts at read level
    total_read_level = sum(read_level_counts.values())

    # Calculate observed COG proportions for each COG
    observed_proportions = {cog: count / total_read_level for cog, count in read_level_counts.items()}

    # Update observed proportions for all COGs
    for cog, proportion in observed_proportions.items():
        observed_proportions_all[cog] = proportion

    # Calculate expected COG proportions for each COG
    for cog, count in gene_level_counts_all.items():
        expected_proportions_all[cog] = count / total_gene_level

    # Sort COGs
    sorted_cogs = sorted(all_cogs)

    # Calculate differences between observed and expected proportions for each COG
    diffs = [observed_proportions_all[cog] - expected_proportions_all[cog] for cog in sorted_cogs]

    # Plot differences
    plt.plot(sorted_cogs, diffs, label=legend_names.get(read_level_file, read_level_file))  # Use custom legend name if available

plt.xlabel('COG Category', fontsize=24)
plt.ylabel('Difference (Observed - Expected)', fontsize=24)
plt.title('Differences between Observed and Expected COG Proportions', fontsize=24)
plt.xticks(rotation=45, ha='right', fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=18)
plt.ylim(-0.20, 0.20)
plt.grid(axis='y')
plt.tight_layout()
plt.show()


