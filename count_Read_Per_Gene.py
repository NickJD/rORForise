import collections

predicted_reads_dict = collections.defaultdict()

# Specify the path to the processed GFF file
gff_file = "Myco_joined_readmode.gff"

with open(gff_file, "r") as file:
    for line in file:
        fields = line.strip().split("\t")
        read = fields[0]
        start = fields[3]
        stop = fields[4]
        frame = fields[6]
        predicted_reads_dict[read] = [start,stop,frame]


bed_reads_dict = collections.defaultdict(dict)

# Specify the path to the processed BED file
bed_file = "mapped_reads_05.bed"

with open(bed_file, "r") as file:
    for line in file:
        fields = line.strip().split("\t")
        read_id = fields[3]  # Assuming the read ID is in the 4th column
        attribute_type = fields[8]  # Assuming the CDS gene name is in the 1st column

        if attribute_type == 'CDS':
            read_start = int(fields[1])  # Start position of the read
            read_end = int(fields[2])  # End position of the read
            cds_start = int(fields[9])  # Start position of the CDS gene
            cds_end = int(fields[10])  # End position of the CDS gene
            gene_id = fields[14].split(';')[0].split(':')[1]
            frame = fields[5]
            #reads_dict[gene_id].append(read_id)

            # Determine alignment position within the gene
            if read_start >= cds_start and read_end <= cds_end:  # Read is completely within the CDS gene
                alignment_position = "Middle"
                position_within_gene = read_start - cds_start
                alignment_start = max(read_start, cds_start)
                alignment_end = min(read_end, cds_end)
            elif read_start < cds_start:  # Read is on the left of the CDS gene
                alignment_position = "Left Edge"
                alignment_start = max(read_start, cds_start)
                alignment_end = min(read_end, cds_end)
            elif read_end > cds_end:  # Read is on the right of the CDS gene
                alignment_position = "Right Edge"
                alignment_start = max(read_start, cds_start)
                alignment_end = min(read_end, cds_end)
            else:  # Read overlaps both ends of the CDS gene
                alignment_position = "Spanning"
                position_within_gene = "N/A"
                alignment_start = max(read_start, cds_start)
                alignment_end = min(read_end, cds_end)

            # Store the read information in the dictionary
            bed_reads_dict[gene_id].update({read_id:[cds_start,cds_end,alignment_position,read_start,read_end,alignment_start,alignment_end,frame]})

# Printing the dictionary
for cds_gene, reads in bed_reads_dict.items():
    print("CDS gene:", cds_gene)
    for read in reads:
        print("Read ID:", read)
        print("Alignment Position:", reads[read][2])
        print("Alignment Start:", reads[read][5])
        print("Alignment End:", reads[read][6])
        print("----")
    print("Number of Reads: ",len(reads))


NoP = 0
Correct = 0
NoCorrect = 0

for gene, reads in bed_reads_dict.items():
    print(gene)
    for read in reads:
        print(read)
        mapped_read_info = reads[read]
        try:
            predicted_read = predicted_reads_dict[read]
            #
            if mapped_read_info[2] == 'Left Edge':
                gene_start = mapped_read_info[0]
                mapped_read_start = mapped_read_info[3]
                if mapped_read_info[-1] == '-':
                    mapped_read_length = mapped_read_info[4] - mapped_read_info[3]
                    read_alignment_start = mapped_read_length+1 - int(predicted_read[1])
                else:
                    read_alignment_start = int(predicted_read[0])
                if mapped_read_start + read_alignment_start == gene_start:
                    print("Found correct start")
                    Correct +=1
                else:
                    print("Not correct start")
                    NoCorrect +=1
        except KeyError:
            print("Read Not Predicted")
            NoP +=1

print(Correct)
print(NoCorrect)
print(NoP)









