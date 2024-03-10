def add_coverage(gene_coverage, start, end):
    for pos in range(start, end + 1):
        gene_coverage.add(pos)


# Initialize dictionaries for gene coverage and gene length
gene_coverage = {}
gene_length = {}

# Read the DIAMOND BLAST TSV file
with open('../Genome_Processing/Staphylococcus_aureus_502A/FragGeneScan/FGS_ART_errFree_paired_DIAMOND_80_99', 'r') as file:
    for line in file:
        qseqid, sseqid, pident, qlen, slen, length, sstart, send = line.strip().split('\t')

        # Convert string to integers for relevant columns
        slen, sstart, send = map(int, [slen, sstart, send])

        # Initialize the gene in the dictionaries if it's not already there
        if sseqid not in gene_coverage:
            gene_coverage[sseqid] = set()
            gene_length[sseqid] = slen

        # Add the coverage for this fragment to the gene's coverage
        add_coverage(gene_coverage[sseqid], min(sstart, send), max(sstart, send))

# Initialize variables for lowest and highest coverage
lowest_coverage = float('inf')
highest_coverage = 0

# Calculate coverage for each gene and accumulate total coverage percentage
total_coverage_percent = 0
for gene_id, positions in gene_coverage.items():
    gene_cov = len(positions)
    gene_cov_percent = (gene_cov / gene_length[gene_id]) * 100

    # Print coverage for each gene
    print(f"Coverage for {gene_id}: {gene_cov_percent:.2f}%")

    # Update lowest and highest coverage
    lowest_coverage = min(lowest_coverage, gene_cov_percent)
    highest_coverage = max(highest_coverage, gene_cov_percent)

    total_coverage_percent += gene_cov_percent

# Calculate average coverage across all genes
average_coverage_percent = total_coverage_percent / len(gene_coverage)
print(f"Average coverage across all genes: {average_coverage_percent:.2f}%")

# Print lowest and highest coverage
print(f"Lowest coverage: {lowest_coverage:.2f}%")
print(f"Highest coverage: {highest_coverage:.2f}%")
