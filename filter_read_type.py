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
        predicted_reads_dict[gff_read] = [gff_start, gff_stop, gff_read_frame] # does this allow for more than one ORF per read?

bed_reads_dict = defaultdict(dict)

prediction_category = defaultdict(list)

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

            try:
                read_ORF = predicted_reads_dict[bed_read_id]
                if bed_CDS_frame == '+' and bed_mapped_read_frame == '+' and read_ORF[2] == '+':
                    prediction_category['Category_1'].append([bed_gene_id,bed_read_id,read_ORF])

                elif bed_CDS_frame == '+' and bed_mapped_read_frame == '+' and read_ORF[2] == '-':
                    prediction_category['Category_2'].append([bed_gene_id,bed_read_id,read_ORF])

                elif bed_CDS_frame == '+' and bed_mapped_read_frame == '-' and read_ORF[2] == '-':
                    prediction_category['Category_3'].append([bed_gene_id,bed_read_id,read_ORF])

                elif bed_CDS_frame == '+' and bed_mapped_read_frame == '-' and read_ORF[2] == '+':
                    prediction_category['Category_4'].append([bed_gene_id,bed_read_id,read_ORF])

                elif bed_CDS_frame == '-' and bed_mapped_read_frame == '+' and read_ORF[2] == '+':
                    prediction_category['Category_5'].append([bed_gene_id,bed_read_id,read_ORF])

                elif bed_CDS_frame == '-' and bed_mapped_read_frame == '+' and read_ORF[2] == '-':
                    prediction_category['Category_6'].append([bed_gene_id,bed_read_id,read_ORF])

                elif bed_CDS_frame == '-' and bed_mapped_read_frame == '-' and read_ORF[2] == '-':
                    prediction_category['Category_7'].append([bed_gene_id,bed_read_id,read_ORF])

                elif bed_CDS_frame == '-' and bed_mapped_read_frame == '-' and read_ORF[2] == '+':
                    prediction_category['Category_8'].append([bed_gene_id,bed_read_id,read_ORF])
                else:
                    print("HMM")
            except KeyError:
                continue


print("And Here")











