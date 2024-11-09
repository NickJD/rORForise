def parse_gff_line(line):
    fields = line.strip().split('\t')
    return (fields[0], int(fields[3]), int(fields[4]), fields[6])

def parse_bed_line(line):
    fields = line.strip().split('\t')
    return (fields[0], int(fields[1]), int(fields[2]), fields[5])

# Parse GFF file (from step 4)
with open('small.gff', 'r') as gff_file:
    gff_data = [parse_gff_line(line) for line in gff_file]

# Parse BED file (from step 3)
with open('small.bed', 'r') as bed_file:
    bed_data = [parse_bed_line(line) for line in bed_file]

# Find intersections
intersections = [gff for gff in gff_data for bed in bed_data if gff[0] == bed[0] and gff[1] == bed[1] and gff[2] == bed[2] and gff[3] == bed[3]]

# Print each FragGeneScan prediction and its correctness
for gene in intersections:
    print(f"FragGeneScan predicted a gene on sequence {gene[0]} with start {gene[1]} and end {gene[2]}. This prediction was confirmed by read alignments.")
