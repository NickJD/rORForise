from collections import defaultdict

predicted_reads_dict = defaultdict()


# Specify the path to the processed GFF file
gff_file = "Myco_joined_readmode.gff"

with open(gff_file, "r") as file:
    for line in file:
        fields = line.strip().split("\t")
        gff_read = fields[0]
        gff_start = int(fields[3])
        gff_stop = int(fields[4])
        gff_read_frame = fields[6]
        predicted_reads_dict[gff_read] = [gff_start, gff_stop, gff_read_frame]

bed_reads_dict = defaultdict(dict)



# Specify the path to the processed BED file
bed_file = "mapped_reads_05.bed"

with open(bed_file, "r") as file:
    for line in file:
        fields = line.strip().split("\t")
        bed_read_id = fields[3]
        attribute_type = fields[8]

        if attribute_type == 'CDS':
            bed_read_start = int(fields[1]) + 1  # Start position of the read
            bed_read_end = int(fields[2]) + 1  # End position of the read
            bed_cds_start = int(fields[9])  # Start position of the CDS gene
            bed_cds_end = int(fields[10])  # End position of the CDS gene
            bed_gene_id = fields[14].split(';')[0].split(':')[1]
            bed_mapped_read_frame = fields[5]
            bed_CDS_frame = fields[12]

            # Determine alignment position within the gene
            if bed_read_start >= bed_cds_start and bed_read_end <= bed_cds_end:  # Read is completely within the CDS gene
                bed_alignment_position = "Middle"
            elif bed_read_start < bed_cds_start and bed_read_end > bed_cds_end:  # Read overlaps both ends of the CDS gene
                bed_alignment_position = "Spanning"
            elif bed_read_start < bed_cds_start and bed_read_end > bed_cds_start:  # Read is on the left of the CDS gene
                bed_alignment_position = "Left Edge"
            elif bed_read_start < bed_cds_end and bed_read_end > bed_cds_end:  # Read is on the right of the CDS gene
                bed_alignment_position = "Right Edge"
            else:
                print("WAT")
            bed_alignment_start = max(bed_read_start, bed_cds_start)
            bed_alignment_end = min(bed_read_end, bed_cds_end)
            # Store the read information in the dictionary
            bed_reads_dict[bed_gene_id][bed_read_id] = (bed_cds_start, bed_cds_end, bed_CDS_frame, bed_alignment_position,
                                                        bed_read_start, bed_read_end, bed_alignment_start, bed_alignment_end)











