import collections

predicted_reads_dict = collections.defaultdict()

# Specify the path to the processed GFF file
gff_file = "Myco_joined_readmode.gff"

with open(gff_file, "r") as file:
    for line in file:
        fields = line.strip().split("\t")
        gff_read = fields[0]
        gff_start = fields[3]
        gff_stop = fields[4]
        gff_read_frame = fields[6]
        predicted_reads_dict[gff_read] = [gff_start,gff_stop,gff_read_frame]


bed_reads_dict = collections.defaultdict(dict)

# Specify the path to the processed BED file
bed_file = "mapped_reads_05.bed"

with open(bed_file, "r") as file:
    for line in file:
        fields = line.strip().split("\t")
        bed_read_id = fields[3]
        attribute_type = fields[8]

        if attribute_type == 'CDS':
            bed_read_start = int(fields[1])  # Start position of the read
            bed_read_end = int(fields[2])  # End position of the read
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
            bed_reads_dict[bed_gene_id][bed_read_id] = [bed_cds_start, bed_cds_end, bed_CDS_frame, bed_alignment_position, bed_read_start, bed_read_end, bed_alignment_start, bed_alignment_end]

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
        (cds_start, cds_end, cds_frame, alignment_position, read_start, read_end, alignment_start, alignment_end) = mapped_read_info
        print(f'Read start: {read_start}, Read end: {read_end}, Alignment position: {alignment_position}, CDS Frame: {cds_frame}')
        try:
            predicted_read = predicted_reads_dict[compared_read]
            (read_orf_start, read_orf_end, read_orf_frame) = predicted_read
            print(f'Predicted ORF start: {read_orf_start}, Predicted ORF end: {read_orf_end}, Predicted ORF mapping frame to CDS Gene: {read_orf_frame}')

            if alignment_position == 'Left Edge':
                # if cds_frame == '-':
                #     mapped_read_length = read_end - read_start
                #     read_prediction_start = mapped_read_length + 1 - int(read_orf_end)
                # else:
                #     read_prediction_start = int(read_orf_start)
                read_prediction_start = int(read_orf_start) # Seems like we do not need to recalculate start/stop
                print(f'Expected gene start: {cds_start}, Read prediction start: {read_prediction_start+read_start}, Mapped read start: {read_start}')
                if read_start + int(read_prediction_start) == cds_start:
                    print("Found correct start")
                    L_Correct += 1
                else:
                    print("Not correct start")
                    L_NoCorrect += 1

            elif alignment_position == 'Right Edge':
                # if cds_frame == '-':
                #     mapped_read_length = read_end - read_start
                #     read_prediction_stop = mapped_read_length + 1 - int(read_orf_start)
                # else:
                #     read_prediction_stop = int(read_orf_end)
                read_prediction_stop = int(read_orf_end) # Seems like we do not need to recalculate start/stop
                print(f'Expected gene stop: {cds_end}, Read prediction stop: {read_prediction_stop+read_start}, Mapped read stop: {read_end}')
                if alignment_start + read_prediction_stop == cds_end:  # and start is either 1/2/3? (+2 is to account for the condon size of 3)
                    print("Found correct stop")
                    R_Correct += 1
                else:
                    print("Not correct stop")
                    R_NoCorrect += 1

            elif alignment_position == 'Middle': # needs work
                mapped_read_length = read_end - read_start
                #if CDS_frame == read_prediction_frame:
                if read_orf_start in [1,2,3] and read_orf_end in [mapped_read_length,mapped_read_length-1,mapped_read_length-2]: # Need to check if it needs to be more than -2 because the stop codon starts -3?
                    print("Found correct position and frame")
                    M_Correct += 1
                else:
                    print("Not correct position and/or frame")
                    M_NoCorrect += 1
                #else:
                print(f'Mapped read start: {alignment_start}, Mapped read stop: {alignment_end}, Read prediction start: {read_orf_start}, Read prediction stop: {read_orf_end}, Correct frame: {cds_frame}') #The frame bit here is not finished

        except KeyError:
            print("Read Not Predicted")
            NoP += 1


print("\n\nBelow is correct and not correct not proportions") # At the moment we are getting very low accuracy so I am a bit unsure what is going on...
print("\nLeft: " + str(L_Correct) + "/" + str(L_NoCorrect) + "\nRight: " + str(R_Correct) + "/" + str(R_NoCorrect) + "\nMiddle: " + str(M_Correct) + "/" + str(M_NoCorrect) + "\nNoPrediction: " + str(NoP))








