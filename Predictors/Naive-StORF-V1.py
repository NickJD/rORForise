import re,sys
import collections


reads_in = open('../Genome_Processing/Staphylococcus_aureus_502A/Processing/ART_Simulated_Reads/Staph_ART_errFree_Combined.fasta', 'r')

predictions_fasta = open('../Genome_Processing/Staphylococcus_aureus_502A/Naive-StORF-V1/Naive-StORF-V1_ART_errFree_Combined.faa','w')
predictions_gff = open('../Genome_Processing/Staphylococcus_aureus_502A/Naive-StORF-V1/Naive-StORF-V1_ART_errFree_Combined.gff','w')


###################
gencode = {
      'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
      'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
      'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
      'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
      'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
      'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
      'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
      'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
      'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
      'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
      'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
      'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
      'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
      'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
      'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
      'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}

def translate_frame(sequence):
    translate = ''.join([gencode.get(sequence[3 * i:3 * i + 3], 'X') for i in range(len(sequence) // 3)])
    return translate

############################

def trim_sequence(seq,direction): # I don't think direction is needed. I have left this in but the run is the same
    remainder = len(seq) % 3
    if remainder == 1:
        if direction == '+':
            return seq[:-1],remainder
        elif direction == '-':
            return seq[:-1], remainder
    elif remainder == 2:
        if direction == '+':
            return seq[:-2],remainder
        elif direction == '-':
            return seq[:-2], remainder
    else:
        return seq,remainder

############################
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

def get_stops(sequence):
    stops = []
    for stop_codon in ['TGA','TAG','TAA']:  # Find all Stops in seq
        stops += [match.start() for match in re.finditer(re.escape(stop_codon), sequence)]
    stops.sort()
    #stops = [stop + 1 for stop in stops]
    return stops

def fasta_load(fasta_in):
    dna_regions = collections.OrderedDict()
    first = True
    #### Default for when presented with standard fasta file
    for line in fasta_in:
        line = line.strip()
        if line.startswith('>') and first == False:  # Check if first seq in file
            dna_regions.update({dna_region_id:seq})
            seq = ''
            dna_region_id = line.split()[0].replace('>', '')
        elif line.startswith('>'):
            seq = ''
            dna_region_id = line.split()[0].replace('>', '')
        else:
            seq += str(line)
            first = False
    dna_regions.update({dna_region_id:seq})



    return dna_regions

sequences = fasta_load(reads_in)


def find_longest_interval(sequence, stop_codon_positions, direction):
    seq_length = len(seq)
    stop_codons_per_frame = {0: set(), 1: set(), 2: set()}
    all_intervals = {0: [], 1: [], 2: []}

    for pos in stop_codon_positions:
        frame = pos % 3
        stop_codons_per_frame[frame].add(pos)

    for frame in range(3):
        stops = sorted(list(stop_codons_per_frame[frame]))
        if not stops:
            all_intervals[frame].append((1, seq_length))
            continue

        if stops[0] != 0:
            all_intervals[frame].append((1, stops[0]))

        for i in range(len(stops) - 1):
            all_intervals[frame].append((stops[i], stops[i + 1]))

        if stops[-1] != seq_length - 1:
            all_intervals[frame].append((stops[-1], seq_length))

    longest_intervals = {}
    for frame, intervals in all_intervals.items():
        longest_interval = max(intervals, key=lambda x: x[1] - x[0])
        longest_interval = list(longest_interval)
        seq_pos = [longest_interval[0], longest_interval[1]]

        if direction == '+':
            if longest_interval[0] == 1 or seq_pos[0] == 1:
                seq_pos[0] = seq_pos[0] - 1
                # longest_interval[0] -= 1
                # current_seq = sequence[longest_interval[0] - 1:longest_interval[1]]
            if longest_interval[1] not in [148, 149, 150]:  # hard coded for now
                # if 150 in list(longest_interval[1]):
                seq_pos[1] = seq_pos[1] + 3
                longest_interval[1] += 3
            #    current_seq = sequence[longest_interval[0]:longest_interval[1]]

            current_seq = sequence[seq_pos[0]:seq_pos[1]]

            current_seq, remainder = trim_sequence(current_seq,direction)
            new_end = longest_interval[1] - remainder
            new_longest_interval = (longest_interval[0], new_end)

        elif direction == '-':
            start, end = longest_interval
            new_start = max(seq_length - end -2, 1)
            if start == 1:
                new_end = seq_length
            else:
                new_end = max(seq_length - start, 1)
            longest_interval = [new_start, new_end]
            #seq_pos = longest_interval

            if longest_interval[1] == 150:# or seq_pos[0] == 1:
                seq_pos[0] = seq_pos[0] - 1
                #longest_interval[0] -= 1
                #current_seq = sequence[longest_interval[0] - 1:longest_interval[1]]
            if longest_interval[0] not in [1, 2, 3]: # hard coded for now
            #if 150 in list(longest_interval[1]):
                seq_pos[1] = seq_pos[1] + 3
                #longest_interval[1] += 3
            #    current_seq = sequence[longest_interval[0]:longest_interval[1]]

            current_seq = sequence[seq_pos[0]:seq_pos[1]]

            current_seq, remainder = trim_sequence(current_seq,direction)
            new_end = longest_interval[1] - remainder
            new_longest_interval = (longest_interval[0], new_end)



        longest_intervals[frame] = {
            'interval': new_longest_interval,
            'sequence': current_seq,
            'sequence_length': len(current_seq)
        }


    return longest_intervals



for id, seq in sequences.items():


    if '@Chromosome-77234/1' in id:
        print()
    # else:
    #     continue

    frames_covered = collections.defaultdict(int)
    #Pos Strand
    stops = get_stops(seq)
    predicted_genes = find_longest_interval(seq,stops,'+')


    rev_seq = revCompIterative(seq)
    rev_stops = get_stops(rev_seq)



    key_mapping = {0: 3, 1: 4, 2: 5}
    tmp_predicted_genes = find_longest_interval(rev_seq,rev_stops,'-')
    tmp_predicted_genes = {key_mapping[old_key]: value for old_key, value in tmp_predicted_genes.items()}
    predicted_genes.update(tmp_predicted_genes)

    longest_prediction = max(predicted_genes, key=lambda k: predicted_genes[k].get('sequence_length', 0))

    longest_interval = predicted_genes[longest_prediction].get('interval')
    longest_sequence = predicted_genes[longest_prediction].get('sequence')

    start_position = longest_interval[0]
    stop_position = longest_interval[1]

    if longest_prediction in [0,1,2]:
        direction = '+'
    elif longest_prediction in [3,4,5]:
        direction = '-'

    if longest_interval[0] == 1:
        if longest_prediction in [0,3]:
            sequence_for_frame = longest_sequence
        elif longest_prediction in [1]:
            sequence_for_frame = longest_sequence[1:]
            sequence_for_frame = sequence_for_frame[:-2]
            start_position += 1
            stop_position -= 2
        elif longest_prediction in [4]:
            if longest_interval[1] in [148,149,150]:
                sequence_for_frame = longest_sequence[1:]
                sequence_for_frame = sequence_for_frame[:-2]
            else:
                sequence_for_frame = longest_sequence
            # sequence_for_frame = longest_sequence[1:] # why different?
            # sequence_for_frame = sequence_for_frame[:-2]
            #sequence_for_frame = longest_sequence
            start_position += 2
            stop_position -= 1
            #sequence_for_frame, remainder = trim_sequence(sequence_for_frame,direction)
            #stop_position = longest_interval[1] - remainder
        elif longest_prediction in [2]:
            sequence_for_frame = longest_sequence[2:]
            sequence_for_frame = sequence_for_frame[:-1]
            start_position += 2
            stop_position -= 1
        elif longest_prediction in [5]:
            if longest_interval[1] in [148,149,150]:
                sequence_for_frame = longest_sequence[2:]
                sequence_for_frame = sequence_for_frame[:-1]
            else:
                sequence_for_frame = longest_sequence
            start_position += 1
            stop_position -= 2

            #sequence_for_frame, remainder = trim_sequence(sequence_for_frame,direction)
            #stop_position = longest_interval[1] - remainder

       # if longest_prediction >= 3: # might need to do +1
       #     corrected_start_position = max(len(seq) - int(stop_position), 1)
       #     corrected_stop_position = max(len(seq) - int(start_position - 1), 1)
        #else:
        corrected_start_position = start_position
        corrected_stop_position = stop_position

    else:
        if longest_prediction in [4]:
            if longest_interval[1] in [148, 149, 150]:
                sequence_for_frame = longest_sequence[1:]
                sequence_for_frame = sequence_for_frame[:-2]
        elif longest_prediction in [5]:
            if longest_interval[1] in [148,149,150]:
                sequence_for_frame = longest_sequence[2:]
                sequence_for_frame = sequence_for_frame[:-1]
        else:
            sequence_for_frame = longest_sequence

        if longest_prediction >= 3:
            corrected_start_position = start_position #+ 1
            corrected_stop_position = stop_position #- 1
        else:
            corrected_start_position = start_position + 1
            corrected_stop_position = stop_position #+ 1



    aa_seq = translate_frame(sequence_for_frame)
    if aa_seq[0] == '*':
        aa_seq = aa_seq[1:]

    if len(aa_seq) >= 20:
        id = id.replace('@','')
        predictions_fasta.write('>'+id+'|'+str(corrected_start_position)+'_'+str(corrected_stop_position)
                          +'|Frame:'+str(longest_prediction+1)
                          +'\n'+aa_seq+'\n')
        if longest_prediction in [0, 1, 2]:
            strand = "+"
        elif longest_prediction in [3, 4, 5]:
            strand = "-"
        predictions_gff.write(id+'\tNS\tCDS\t'+str(corrected_start_position)+'\t'+str(corrected_stop_position)+
                              '\t.\t'+strand+'\t.\t'+id+'|'+str(corrected_start_position)+'_'+str(corrected_stop_position)
                          +'|Frame:'+str(longest_prediction+1)+'\n')


    else:
        print("Sequence under 20 aa")




