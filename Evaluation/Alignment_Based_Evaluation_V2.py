from collections import defaultdict
import statistics

def add_coverage(gene_coverage, start, end):
    for pos in range(start, end + 1):
        gene_coverage[pos] += 1

# Initialize dictionaries for gene coverage and gene length
gene_coverage = defaultdict(lambda: defaultdict(int))
gene_length = {}

# Initialize variables for counting predictions, reads, and hits
total_predictions = 0
reads_with_predictions = set()
reads_with_multiple_predictions = 0
predictions_with_hit = 0
total_genes_with_prediction = 0

# Read the DIAMOND BLAST TSV file
with open('../Genome_Processing/Mycoplasma_genitalium_G37/'
           'Pyrodigal/Pyrodigal_ART_errFree_Combined_DP_80', 'r') as blast_file:
# with open('../Genome_Processing/Mycoplasma_genitalium_G37/'
#           'Myco_Reads_Diamond_Peps', 'r') as blast_file:
    for line in blast_file:
        fields = line.strip().split('\t')
        qseqid, sseqid, pident, qlen, slen, length, sstart, send, score = line.strip().split('\t')

        # Track the number of predictions with a hit
        if qseqid != sseqid:
            predictions_with_hit += 1

        # Track the number of genes or hits with a prediction
        if sseqid not in gene_length:
            gene_length[sseqid] = int(slen)
            total_genes_with_prediction += 1

        # Add the coverage for this hit to the gene's coverage
        add_coverage(gene_coverage[sseqid], int(sstart), int(send))

# Read the FASTA file containing predictions
with open('../Genome_Processing/Mycoplasma_genitalium_G37/'
          'Pyrodigal/Pyrodigal_ART_errFree_Combined.faa', 'r') as fasta_file:
    for line in fasta_file:
        if line.startswith('>'):
            total_predictions += 1
            read_id = line.strip().split('/')[0]
            reads_with_predictions.add(read_id)
            if '_' in read_id:
                reads_with_multiple_predictions += 1

# Calculate average number of predictions per base across all genes

import statistics

calculated_coverage = []
missing_position_counts = []
genes_with_single_mapped_position = 0
positions_with_single_mapping = 0
total_positions = 0

for gene, length in gene_length.items():
    coverage_numbers = gene_coverage.get(gene, {})  # Get coverage for the gene; default to empty dict if not present
    if len(coverage_numbers) == length:  # Check if coverage information is available for all positions
        median_coverage = statistics.mean(coverage_numbers.values())  # Calculate mean coverage
        calculated_coverage.append(median_coverage)

        # Count positions with single mapping
        positions_with_single_mapping += sum(1 for cov in coverage_numbers.values() if cov == 1)
        total_positions += length

        if 1 in coverage_numbers.values():  # Check if any position has only one read mapped
            genes_with_single_mapped_position += 1
    else:
        # If coverage information is incomplete, calculate mean coverage considering zero coverage for missing positions
        total_positions += length
        covered_positions = set(coverage_numbers.keys())  # Positions with coverage information
        missing_positions = length - len(covered_positions)  # Count missing positions
        missing_position_counts.append(missing_positions)
        zero_coverage_positions = {pos: 0 for pos in range(1, length + 1) if
                                   pos not in covered_positions}  # Fill missing positions with zero coverage
        all_positions_coverage = {**coverage_numbers, **zero_coverage_positions}  # Combine coverage data
        median_coverage = statistics.mean(all_positions_coverage.values())  # Calculate mean coverage
        calculated_coverage.append(median_coverage)

mean_coverage = statistics.mean(calculated_coverage)
mean_missing_positions = statistics.mean(missing_position_counts)
proportion_genes_with_single_mapped_position = genes_with_single_mapped_position / len(gene_length)
proportion_positions_with_single_mapping = positions_with_single_mapping / total_positions

print(f"Mean Coverage: {mean_coverage}")
print(f"Mean Missing Positions per Gene: {mean_missing_positions}")
print(f"Number of Genes with Missing Positions: {missing_positions}")

print(
    f"Proportion of Genes with At Least One Position with Only One Mapping: {proportion_genes_with_single_mapped_position:.2%}")
if proportion_genes_with_single_mapped_position > 0:
    proportion_positions_with_single_mapping_filtered = positions_with_single_mapping / (total_positions * proportion_genes_with_single_mapped_position)
    print(f"Proportion of Positions with Only One Mapping Across Genes with Only One Mapping: {proportion_positions_with_single_mapping_filtered:.2%}")
else:
    print("Proportion of Positions with Only One Mapping Across Genes with Only One Mapping: N/A")



total_bases = sum(gene_length.values())
total_predictions_bases = sum(sum(gene_coverage[gene].values()) for gene in gene_coverage)
average_predictions_per_base = total_predictions_bases / total_bases

# Print results
print(f"Number of predictions: {total_predictions}")
print(f"Number of reads with multiple predictions: {reads_with_multiple_predictions}")
print(f"Number of predictions with a hit: {predictions_with_hit}")
print(f"Number of genes or hits with a prediction: {total_genes_with_prediction}")
print(f"Average number of predictions per base across all genes: {average_predictions_per_base:.2f}")
