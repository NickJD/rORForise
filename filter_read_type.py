from collections import defaultdict

predicted_reads_dict = defaultdict(list)


# Specify the path to the processed GFF file
gff_file = "./Genome_Processing/Mycoplasma/Processing/FragGeneScan/Myco_ART_errFree_paired.gff"

with open(gff_file, "r") as file:
    for line in file:
        fields = line.strip().split("\t")
        if len(fields) == 9:
            gff_read = fields[0].replace('@','')
            gff_start = int(fields[3])
            gff_stop = int(fields[4])
            gff_read_frame = fields[6]
            prediction = fields[8].split('=')[1].split(';')[0]
            predicted_reads_dict[gff_read].append([prediction, gff_start, gff_stop, gff_read_frame])

bed_reads_dict = defaultdict(dict)

prediction_category = defaultdict(list)

# Creating a nested defaultdict with lists as values
predictions = defaultdict(lambda: defaultdict(list))

# Specify the path to the processed BED file
bed_file = "./Genome_Processing/Mycoplasma/Processing/Myco_Reads_Intersect.bed"

catch = []

with open(bed_file, "r") as file:
    for line in file:
        fields = line.strip().split("\t")
        bed_read_id = fields[3]
        attribute_type = fields[14]

        if attribute_type == 'CDS':
            bed_read_start = int(fields[1]) + 1  # Start position of the read
            bed_read_end = int(fields[2]) + 1  # End position of the read
            bed_cds_start = int(fields[15])  # Start position of the CDS gene
            bed_cds_end = int(fields[16])  # End position of the CDS gene
            bed_gene_id = fields[20].split(';')[0].split(':')[1]
            bed_mapped_read_frame = fields[5]
            bed_CDS_frame = fields[18]

            try:
                prediction = predicted_reads_dict[bed_read_id]
                for pred in prediction:
                    if pred[0] not in predictions:
                        print("")
                    else:
                        print("D$$K")
                    pred.extend([bed_read_start,bed_read_end])
                    if bed_CDS_frame == '+' and bed_mapped_read_frame == '+' and pred[3] == '+':
                        prediction_category['Category_1'].append([bed_gene_id,bed_read_id,prediction])
                        predictions[pred[0]][bed_gene_id].append(['Category_1', pred])

                    elif bed_CDS_frame == '+' and bed_mapped_read_frame == '+' and pred[3] == '-':
                        prediction_category['Category_2'].append([bed_gene_id,bed_read_id,prediction])
                        predictions[pred[0]][bed_gene_id].append(['Category_2', pred])

                    elif bed_CDS_frame == '+' and bed_mapped_read_frame == '-' and pred[3] == '-':
                        prediction_category['Category_3'].append([bed_gene_id,bed_read_id,pred])

                        predictions[pred[0]][bed_gene_id].append(['Category_3', pred])

                    elif bed_CDS_frame == '+' and bed_mapped_read_frame == '-' and pred[3] == '+':
                        prediction_category['Category_4'].append([bed_gene_id,bed_read_id,pred])

                        predictions[pred[0]][bed_gene_id].append(['Category_4', pred])

                    elif bed_CDS_frame == '-' and bed_mapped_read_frame == '+' and pred[3] == '+':
                        prediction_category['Category_5'].append([bed_gene_id,bed_read_id,pred])

                        predictions[pred[0]][bed_gene_id].append(['Category_5', pred])

                    elif bed_CDS_frame == '-' and bed_mapped_read_frame == '+' and pred[3] == '-':
                        prediction_category['Category_6'].append([bed_gene_id,bed_read_id,pred])

                        predictions[pred[0]][bed_gene_id].append(['Category_6', pred])

                    elif bed_CDS_frame == '-' and bed_mapped_read_frame == '-' and pred[3] == '-':
                        prediction_category['Category_7'].append([bed_gene_id,bed_read_id,pred])

                        predictions[pred[0]][bed_gene_id].append(['Category_7', pred])

                    elif bed_CDS_frame == '-' and bed_mapped_read_frame == '-' and pred[3] == '+':
                        prediction_category['Category_8'].append([bed_gene_id,bed_read_id,pred])

                        predictions[pred[0]][bed_gene_id].append(['Category_8', pred])
                    else:
                        print("HMM")
                    print("Another")
            except KeyError:
                continue


print("And Here")

count = 0
# Assuming predictions is a dictionary of dictionaries
for key, sub_dict in predictions.items():
    if sum(len(value) for value in sub_dict.values()) > 1:
        print(f"Sub-dictionary at key {key} has more than one value:")
        count += 1
        for sub_key, sub_value in sub_dict.items():
            print(f"    {sub_key}: {sub_value}")


print(count)







