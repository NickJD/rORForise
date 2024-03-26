def add_coverage(gene_coverage, start, end):
    for pos in range(start, end + 1):
        gene_coverage.add(pos)

# Initialize dictionaries for gene coverage and gene length
gene_coverage = {}
gene_length = {}

# Initialize variables for counting predictions, reads, and hits
total_predictions = 0
reads_with_predictions = set()
reads_with_multiple_predictions = 0
predictions_with_hit = 0
total_genes_with_prediction = 0

# Initialize counters for coverage percentages
coverage_counters = {25: 0, 50: 0, 75: 0, 90: 0, 99: 0, 100: 0}

# Read the DIAMOND BLAST TSV file
with open('../Genome_Processing/Staphylococcus_aureus_502A/'
          'Naive-StORF-V3/Naive-StORF-V3_ART_errFree_Combined_DP_80', 'r') as blast_file:
    for line in blast_file:
        fields = line.strip().split('\t')
        qseqid, sseqid, pident, qlen, slen, length, sstart, send, score = line.strip().split('\t')

        # Track the number of predictions with a hit
        if qseqid != sseqid:
            predictions_with_hit += 1

        # Track the number of genes or hits with a prediction
        if sseqid not in gene_coverage:
            gene_coverage[sseqid] = set()
            total_genes_with_prediction += 1

        # Add the coverage for this hit to the gene's coverage
        add_coverage(gene_coverage[sseqid], int(sstart), int(send))

# Read the FASTA file containing predictions
with open('../Genome_Processing/Staphylococcus_aureus_502A/'
          'Naive-StORF-V3/Naive-StORF-V3_ART_errFree_Combined.faa', 'r') as fasta_file:
    for line in fasta_file:
        if line.startswith('>'):
            total_predictions += 1
            read_id = line.strip().split('/')[0]
            reads_with_predictions.add(read_id)
            if '_' in read_id:
                reads_with_multiple_predictions += 1

# Calculate coverage statistics for each gene
lowest_coverage = float('inf')
highest_coverage = 0
total_coverage_percent = 0

for gene_id, positions in gene_coverage.items():
    gene_length[gene_id] = max(positions)
    gene_cov = len(positions)
    gene_cov_percent = (gene_cov / gene_length[gene_id]) * 100

    # Update coverage counters
    for threshold in coverage_counters.keys():
        if gene_cov_percent >= threshold:
            coverage_counters[threshold] += 1

    # Update lowest and highest coverage
    lowest_coverage = min(lowest_coverage, gene_cov_percent)
    highest_coverage = max(highest_coverage, gene_cov_percent)
    total_coverage_percent += gene_cov_percent

# Calculate average coverage across all genes
average_coverage_percent = total_coverage_percent / len(gene_coverage) if len(gene_coverage) > 0 else 0

# Calculate the proportion of times a read with multiple predictions had those extra predictions with a hit
if reads_with_multiple_predictions > 0:
    proportion_extra_predictions_with_hit = (predictions_with_hit - reads_with_multiple_predictions) / reads_with_multiple_predictions
else:
    proportion_extra_predictions_with_hit = 0

# Print results
print(f"Number of predictions: {total_predictions}")
print(f"Number of reads with multiple predictions: {reads_with_multiple_predictions}")
print(f"Number of predictions with a hit: {predictions_with_hit}")
print(f"Number of genes or hits with a prediction: {total_genes_with_prediction}")
print("Coverage statistics:")
print(f"Average coverage across all genes: {average_coverage_percent:.2f}%")
print(f"Lowest coverage: {lowest_coverage:.2f}%")
print(f"Highest coverage: {highest_coverage:.2f}%")
for threshold, count in coverage_counters.items():
    print(f"Number of genes with at least {threshold}% coverage: {count}")
print(f"Proportion of times a read with multiple predictions had those extra predictions with a hit: {proportion_extra_predictions_with_hit:.2f}")
