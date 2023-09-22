import collections

predicted_reads_dict = collections.defaultdict()

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

bed_reads_dict = collections.defaultdict(dict)

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
            bed_reads_dict[bed_gene_id][bed_read_id] = (
            bed_cds_start, bed_cds_end, bed_CDS_frame, bed_alignment_position, bed_read_start, bed_read_end,
            bed_alignment_start, bed_alignment_end)

NoP = 0
L_Correct = 0
L_Alternative_Starts = collections.defaultdict(
    str)  # Could just point per alt start codon but might expand on this later
L_NoCorrect = 0

R_Correct = 0
R_NoCorrect = 0

M_Correct = 0
M_NoCorrect = 0

S_Correct = 0
S_NoCorrect = 0
S_Alternative_Starts = collections.defaultdict(str)

for gene, alignment_information in bed_reads_dict.items():
    print(f'Gene: {gene}')
    for compared_read in alignment_information:
        if 'Chromosome-49188/1' in compared_read:
            print(2)
        print(
            f'Read: {compared_read}')  # cds_start,cds_end,cd_frame, alignment_position,read_start,read_end,alignment_start,alignment_end
        mapped_read_info = alignment_information[compared_read]
        (cds_start, cds_end, cds_frame, alignment_position, read_start, read_end, alignment_start,
         alignment_end) = mapped_read_info
        print(
            f'Read start: {read_start}, Read end: {read_end}, Alignment position: {alignment_position}, CDS Frame: {cds_frame}')
        try:
            predicted_read = predicted_reads_dict[compared_read]
            (read_orf_start, read_orf_end, read_orf_frame) = predicted_read
            print(
                f'Predicted ORF start: {read_orf_start}, Predicted ORF end: {read_orf_end}, Predicted ORF mapping frame to CDS Gene: {read_orf_frame}')

            if alignment_position == 'Left Edge':
                print("Left Edge")
                print(
                    f'Expected gene start: {cds_start}, Read prediction start: {read_orf_start + read_start}, Mapped read start: {read_start}')
                Left_read_start = read_start + read_orf_start - 1
                diff = cds_start - Left_read_start
                if Left_read_start == cds_start:  # if diff == 0
                    print("Found correct start")
                    L_Correct += 1
                elif diff % 3 == 0:
                    print("Correct Frame predicted but wrong start codon selected")
                    L_Alternative_Starts[compared_read] = diff  # Left_read_start
                else:
                    print("Not correct start")
                    L_NoCorrect += 1

            elif alignment_position == 'Right Edge':
                print("Right Edge")
                print(
                    f'Expected gene stop: {cds_end}, Read prediction stop: {read_orf_end + read_start}, Mapped read stop: {read_end}')
                if alignment_start + read_orf_end - 1 == cds_end:
                    print("Found correct stop")
                    R_Correct += 1
                else:
                    print("Not correct stop")
                    R_NoCorrect += 1

            elif alignment_position == 'Middle':
                print("Middle")
                mapped_read_length = read_end - read_start + 1
                if read_orf_start in [1, 2, 3] and read_orf_end > mapped_read_length - 3:
                    ## Checking if predicted 'frame' is divisible and therefore 'in-frame' with genome cds_start
                    middle_read_start = read_start + read_orf_start - 1
                    diff = middle_read_start - cds_start
                    if diff % 3 == 0:
                        print("We think this is the correct Frame")
                        M_Correct += 1
                    else:  # afc
                        print("We do not think this is the correct Frame")
                        M_NoCorrect += 1
                else:
                    print("We do not think this is the correct Frame")
                    M_NoCorrect += 1
                print(
                    f'Mapped read start: {alignment_start}, Mapped read stop: {alignment_end}, Read prediction start: {read_orf_start}, Read prediction stop: {read_orf_end}, Correct frame: {cds_frame}')  # The frame bit here is not finished

            else:  # alignment_position == 'Spanning'
                print('Spanning')
                if read_start + read_orf_end - 1 == cds_end:
                    Left_read_start = read_start + read_orf_start - 1
                    start_diff = cds_start - Left_read_start
                    if start_diff == 0:
                        print("We think this is the correct Frame")
                        S_Correct += 1
                    elif start_diff % 3 == 0:
                        print("Correct Frame predicted but wrong start codon selected")
                        S_Alternative_Starts[compared_read] = start_diff  # Left_read_start
                    else:
                        print("We do not think this is the correct Frame")
                        S_NoCorrect += 1
                else:
                    print("We do not think this is the correct Frame")
                    S_NoCorrect += 1

        except KeyError:
            print("Read Not Predicted")
            NoP += 1

# At the moment we are getting very low accuracy so I am a bit unsure what is going on...

print("\n\nBelow is correct number of predictions out of the total for each category")
print("\nLeft: " + str(L_Correct) + "/" + str(L_NoCorrect + L_Correct))
print("Right: " + str(R_Correct) + "/" + str(R_NoCorrect + R_Correct))
print("Middle: " + str(M_Correct) + "/" + str(M_NoCorrect + M_Correct))
print("Spanning: " + str(S_Correct) + "/" + str(S_NoCorrect + S_Correct))
print("NoPrediction: " + str(NoP))

l_diffs = L_Alternative_Starts.values()
print(f"Left correct frame but wrong start codon: {len(l_diffs)}")
print("Left Start differences:", sorted(l_diffs))

s_diffs = S_Alternative_Starts.values()
print("Spanning Start differences:", sorted(s_diffs))





