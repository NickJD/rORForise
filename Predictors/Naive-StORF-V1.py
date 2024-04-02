import re,sys,random
import collections


reads_in = open('../Genome_Processing/Mycoplasma_genitalium_G37//Processing/ART_Simulated_Reads/Myco_ART_errFree_Combined.fasta', 'r')

predictions_fasta = open('../Genome_Processing/Mycoplasma_genitalium_G37/Naive-StORF-V1/Naive-StORF-V1_ART_errFree_Combined.faa','w')
predictions_gff = open('../Genome_Processing/Mycoplasma_genitalium_G37/Naive-StORF-V1/Naive-StORF-V1_ART_errFree_Combined.gff','w')


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

def trim_sequence(seq):
    remainder = len(seq) % 3
    if remainder == 1:
        return seq[:-1],remainder
    elif remainder == 2:
        return seq[:-2],remainder
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


def find_longest_interval(sequence, stop_codon_positions):
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
        #
        if longest_interval[0] == 1:
            current_seq = sequence[longest_interval[0] - 1:longest_interval[1]]
        else:
            current_seq = sequence[longest_interval[0]:longest_interval[1]]

        current_seq,remainder = trim_sequence(current_seq)
        new_end = longest_interval[1] - remainder
        new_longest_interval = (longest_interval[0], new_end)

        #aa_seq = translate_frame(current_seq)
        longest_intervals[frame] = {
            'interval': new_longest_interval,
            'sequence': current_seq,
            'sequence_length': len(current_seq)
        }


    return longest_intervals



for id, seq in sequences.items():
    frames_covered = collections.defaultdict(int)
    #Pos Strand

    stops = get_stops(seq)

    if id == '@Chromosome-77292/1':
        print("answer")
    if id == '@Chromosome-77336/1':
        print("ss")

    predicted_genes = find_longest_interval(seq,stops)


    rev_seq = revCompIterative(seq)
    rev_stops = get_stops(rev_seq)

    key_mapping = {0: 3, 1: 4, 2: 5}
    tmp_predicted_genes = find_longest_interval(rev_seq,rev_stops)
    tmp_predicted_genes = {key_mapping[old_key]: value for old_key, value in tmp_predicted_genes.items()}
    predicted_genes.update(tmp_predicted_genes)

    longest_prediction = max(predicted_genes, key=lambda k: predicted_genes[k].get('sequence_length', 0))
    #longest_prediction = random.randint(0, 5)

    # # Get the unique 'sequence_length' values and sort them to find the second longest
    # unique_lengths = sorted(set(gene_data['sequence_length'] for gene_data in predicted_genes.values()), reverse=True)
    # second_longest_length = unique_lengths[1] if len(unique_lengths) > 1 else None
    #
    # # Find the key corresponding to the second longest 'sequence_length'
    # second_longest_prediction = next(
    #     (key for key, val in predicted_genes.items() if val['sequence_length'] == second_longest_length), None)

    longest_interval = predicted_genes[longest_prediction].get('interval')
    longest_sequence = predicted_genes[longest_prediction].get('sequence')

    start_position = longest_interval[0]
    stop_position = longest_interval[1]

    if id == '@Chromosome-77292/1':
        print("answer")

    if longest_interval[0] == 1:
        if longest_prediction in [0,3]:
            sequence_for_frame = longest_sequence
            #start_position += 1
        elif longest_prediction in [1,4]:
            sequence_for_frame = longest_sequence[1:]
            start_position += 1
            sequence_for_frame, remainder = trim_sequence(sequence_for_frame)
            stop_position = longest_interval[1] - remainder
            #start_position +=2
        elif longest_prediction in [2,5]:
            #start_position +=3
            sequence_for_frame = longest_sequence[2:]
            start_position += 2
            sequence_for_frame, remainder = trim_sequence(sequence_for_frame)
            stop_position = longest_interval[1] - remainder
            #new_longest_interval = (longest_interval[0], new_end)

        if longest_prediction >= 3:
            corrected_start_position = max(len(seq) - int(stop_position - 1), 1)
            corrected_stop_position = max(len(seq) - int(start_position - 1), 1)
        else:
            corrected_start_position = start_position
            corrected_stop_position = stop_position


    else:
        sequence_for_frame = longest_sequence

        if longest_prediction >= 3:
            corrected_start_position = start_position #+ 1
            corrected_stop_position = stop_position #- 1
        else:
            corrected_start_position = start_position #+ 1
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




