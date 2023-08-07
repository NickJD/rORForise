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
        read_frame = fields[6]
        predicted_reads_dict[read] = [start,stop,read_frame]


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
            mapped_read_frame = fields[5]
            CDS_frame = fields[12]
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
            bed_reads_dict[gene_id].update({read_id:[cds_start,cds_end,CDS_frame,alignment_position,read_start,read_end,alignment_start,alignment_end]})

# Printing the dictionary
# for cds_gene, reads in bed_reads_dict.items():
#     print("CDS gene:", cds_gene)
#     for read in reads:
#         print("Read ID:", read)
#         print("Alignment Position:", reads[read][2])
#         print("Alignment Start:", reads[read][5])
#         print("Alignment End:", reads[read][6])
#         print("----")
#     print("Number of Reads: ",len(reads))


NoP = 0
L_Correct = 0
L_NoCorrect = 0

R_Correct = 0
R_NoCorrect = 0

M_Correct = 0
M_NoCorrect = 0

for gene, alignment_information in bed_reads_dict.items():
    print(f'Gene: {gene}')
    for compared_read in alignment_information:
        print(f'Read: {compared_read}')  # cds_start,cds_end,alignment_position,read_start,read_end,alignment_start,alignment_end,frame
        mapped_read_info = alignment_information[compared_read]
        print(f'Read start: {mapped_read_info[4]}, Read end: {mapped_read_info[5]}, Alignment position: {mapped_read_info[3]}, CDS Frame: {mapped_read_info[2]}')
        try:
            predicted_read = predicted_reads_dict[compared_read]
            print(f'Predicted ORF start: {predicted_read[0]}, Predicted ORF end: {predicted_read[1]}, Predicted ORF mapping frame to CDS Gene: {predicted_read[2]}')

            if mapped_read_info[3] == 'Left Edge':
                gene_start = mapped_read_info[0]
                mapped_read_start = mapped_read_info[4]
                if mapped_read_info[2] == '-':
                    mapped_read_length = mapped_read_info[5] - mapped_read_info[4]
                    read_prediction_start = mapped_read_length + 1 - int(predicted_read[1])
                else:
                    read_prediction_start = int(predicted_read[0])
                print(
                    f'Expected gene start: {gene_start}, Read prediction start: {read_prediction_start+mapped_read_start}, Mapped read start: {mapped_read_start}')
                if mapped_read_start + read_prediction_start == gene_start:
                    print("Found correct start")
                    L_Correct += 1
                else:
                    print("Not correct start")
                    L_NoCorrect += 1

            elif mapped_read_info[3] == 'Right Edge':
                gene_start = mapped_read_info[0]
                gene_stop = mapped_read_info[1]
                mapped_read_stop = mapped_read_info[5]
                mapped_read_start = mapped_read_info[4]
                if mapped_read_info[2] == '-':
                    mapped_read_length = mapped_read_info[5] - mapped_read_info[4]
                    read_prediction_stop = mapped_read_length + 1 - int(predicted_read[0])
                else:
                    read_prediction_stop = int(predicted_read[1])
                print(
                    f'Expected gene stop: {gene_stop}, Read prediction stop: {read_prediction_stop+mapped_read_start}, Mapped read stop: {mapped_read_stop}')
                if mapped_read_start + read_prediction_stop == gene_stop:  # and start is either 1/2/3? (+2 is to account for the condon size of 3)
                    print("Found correct stop")
                    R_Correct += 1
                else:
                    print("Not correct stop")
                    R_NoCorrect += 1

            elif mapped_read_info[3] == 'Middle':
                mapped_read_start = mapped_read_info[4]
                mapped_read_stop = mapped_read_info[5]
                read_prediction_start = int(predicted_read[0])
                read_prediction_stop = int(predicted_read[1])
                read_prediction_frame = predicted_read[2]
                CDS_frame = mapped_read_info[2]
                mapped_read_length = mapped_read_info[5] - mapped_read_info[4]
                #if CDS_frame == read_prediction_frame:
                if read_prediction_start in [1,2,3] and read_prediction_stop in [mapped_read_length,mapped_read_length-1,mapped_read_length-2]:
                    print("Found correct position and frame")
                    M_Correct += 1
                else:
                    print("Not correct position and/or frame")
                    M_NoCorrect += 1
                #else:
                print(f'Mapped read start: {mapped_read_start}, Mapped read stop: {mapped_read_stop}, Read prediction start: {read_prediction_start}, Read prediction stop: {read_prediction_stop}, Correct frame: {CDS_frame}')

        except KeyError:
            print("Read Not Predicted")
            NoP += 1



print(L_Correct)
print(R_Correct)
print(M_Correct)
print(L_NoCorrect)
print(R_NoCorrect)
print(M_NoCorrect)
print(NoP)









