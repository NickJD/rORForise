import sys
from collections import defaultdict

def revCompIterative(watson): #Gets Reverse Complement
    complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
                   'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M',
                   'M': 'K', 'V': 'B', 'B': 'V', 'H': 'D', 'D': 'H',
                   'U': 'U', 'i': 'i', '-':'-'}
    watson = watson.upper()
    watson_rev = watson[::-1]
    try:
        crick = ''.join([complements[nt] for nt in watson_rev])
    except KeyError as e:
        missing_key = e.args[0]
        print("FASTA sequence has unexected character: " + missing_key, file=sys.stderr)
        sys.exit()
    return crick

with open("Myco.fa", "r") as fasta_file:
    sequence = "".join(line.strip() for line in fasta_file if not line.startswith(">"))
rev_sequence = revCompIterative(sequence)
seq_length = len(sequence)

def getCodons(frame,seq_length,sequence,rev_sequence,start,stop, codon_type, anno_type):
    startCodon, stopCodon = None,None
    if '-' in frame:  # Reverse Compliment starts and stops adjusted
        if "Start_Codon" in codon_type:
            #if "CDS" in anno_type:
            r_start = seq_length - stop
            startCodon = rev_sequence[r_start - 2:r_start + 1]
            #else: # Yes this needs to be cleaned
            #    startCodon = rev_sequence[start - 1:start + 2]
            i = rev_sequence[r_start-3:r_start+9]
        else:
            r_stop = seq_length - start
            stopCodon = rev_sequence[r_stop - 2:r_stop + 1]
    elif '+' in frame:
        if "Start_Codon" in codon_type:
            startCodon = sequence[start - 1:start + 2]
        else:
            stopCodon = sequence[stop - 3:stop]
    return startCodon, stopCodon

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

NoP = 0
L_Correct = 0
L_Alternative_Starts, L_Alternative_Starts_runoff  = defaultdict(int), defaultdict(int)
L_NoCorrect = 0

R_Correct = 0
R_Alternative_Stops, R_Alternative_Stops_runoff  = defaultdict(int), defaultdict(int)
R_NoCorrect = 0

M_Correct = 0
M_NoCorrect = 0

S_Correct = 0
S_NoCorrect = 0
S_Alternative_Starts = defaultdict(int)

for gene, alignment_information in bed_reads_dict.items():
    print(f'Gene: {gene}')
    for compared_read in alignment_information:
        print(f'Read: {compared_read}')  # cds_start,cds_end,cd_frame, alignment_position,read_start,read_end,alignment_start,alignment_end
        mapped_read_info = alignment_information[compared_read]
        (cds_start, cds_end, cds_frame, alignment_position, read_start, read_end, alignment_start,
         alignment_end) = mapped_read_info
        print(f'Read start: {read_start}, Read end: {read_end}, Alignment position: {alignment_position}, CDS Frame: {cds_frame}')
        if 'Chromosome-32366/1' in compared_read: # This is the funky read
            print("$$")
        try:
            predicted_orf = predicted_reads_dict[compared_read]
            (read_orf_start, read_orf_end, read_orf_frame) = predicted_orf
            mapped_read_length = read_end - read_start + 1
            if '-' in read_orf_frame: # This needs to be up here for the print statements
                r_start = mapped_read_length - read_orf_end + 1 # maybe we don't need the +1 / -1
                r_stop = mapped_read_length - read_orf_start + 1
                Left_read_orf_start = read_start + r_start - 1
                Left_read_orf_stop = read_start + r_stop - 1
            else:
                Left_read_orf_start = read_start + read_orf_start - 1
                Left_read_orf_stop = read_start + read_orf_end - 1
            orf_start_codon = sequence[Left_read_orf_start - 1:Left_read_orf_start + 2]  # this is working - mostly

            print(f'Predicted ORF start: {Left_read_orf_start}, Predicted ORF end: {Left_read_orf_stop}, Predicted ORF mapping frame to CDS Gene: {read_orf_frame}')

            if alignment_position == 'Left Edge':
                print("Left Edge")
                print(f'Expected gene start: {cds_start}, Read prediction start: {Left_read_orf_start}, Mapped read start: {read_start}')
                diff = cds_start - Left_read_orf_start
                if diff == 0:
                    print("Found correct start")
                    L_Correct += 1
                elif diff % 3 == 0:
                    # Function to append 1 to subkey
                    print("Correct Frame predicted but wrong start codon selected")
                    cds_start_codon = getCodons(cds_frame, seq_length, sequence, rev_sequence, cds_start, cds_end, "Start_Codon", "CDS")
                    #getCodons(read_orf_frame, seq_length, sequence, rev_sequence, Left_read_orf_start, Left_read_orf_stop, "Start_Codon", "ORF")
                    alt_codon = cds_start_codon[0] + '_' + orf_start_codon
                    if read_orf_start in [1, 2, 3] or read_orf_end >= mapped_read_length - 3: # the or bit accounts for negative ORFs
                        L_Alternative_Starts_runoff[alt_codon] += 1
                    else:
                        L_Alternative_Starts[alt_codon] +=1
                    L_NoCorrect += 1
                else:
                    print("Not correct start")
                    L_NoCorrect += 1

            elif alignment_position == 'Right Edge':
                print("Right Edge")
                print(f'Expected gene stop: {cds_end}, Read prediction stop: {read_orf_end + read_start}, Mapped read stop: {read_end}')
                Right_read_orf_stop = read_start + read_orf_end - 1
                diff = cds_end - Right_read_orf_stop
                if Right_read_orf_stop == cds_end:  # if diff == 0
                    print("Found correct stop")
                    R_Correct += 1
                elif diff % 3 == 0:
                    # Function to append 1 to subkey
                    print("Correct Frame predicted but wrong stop codon selected")
                    # cds_start_codon = getCodons(cds_frame, seq_length, sequence, rev_sequence, None, cds_end)
                    # orf_start_codon = getCodons(read_orf_frame, seq_length, sequence, rev_sequence, None, Right_read_orf_stop)
                    # alt_codon = cds_start_codon[1] + '_' + orf_start_codon[1]
                    # if read_orf_end > mapped_read_length - 3: # If a runoff prediction, technically not a stop codon
                    #     R_Alternative_Stops_runoff[alt_codon] += 1
                    # else:
                    #     R_Alternative_Stops[alt_codon] +=1
                    R_NoCorrect += 1
                else:
                    print("Not correct stop")
                    R_NoCorrect += 1

            elif alignment_position == 'Middle':
                print("Middle")
                if read_orf_start in [1, 2, 3] and read_orf_end >= mapped_read_length - 3:
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
                print(f'Mapped read start: {alignment_start}, Mapped read stop: {alignment_end}, Read prediction start: {read_orf_start}, Read prediction stop: {read_orf_end}, Correct frame: {cds_frame}')  # The frame bit here is not finished

            else:  # alignment_position == 'Spanning'
                print('Spanning')
                if read_start + read_orf_end - 1 == cds_end:
                    Left_read_orf_start = read_start + read_orf_start - 1
                    start_diff = cds_start - Left_read_orf_start
                    if start_diff == 0:
                        print("We think this is the correct Frame")
                        S_Correct += 1
                    elif start_diff % 3 == 0:
                        print("Correct Frame predicted but wrong start codon selected")
                        S_Alternative_Starts[compared_read] = start_diff  # Left_read_orf_start
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

l_diffs = len(L_Alternative_Starts.values()) + len(L_Alternative_Starts_runoff.values())
print(f"Left correct frame but wrong start codon: " + str(l_diffs))
#print("Left Start differences:", sorted(l_diffs))

s_diffs = S_Alternative_Starts.values()
print("Spanning Start differences:", sorted(s_diffs))





