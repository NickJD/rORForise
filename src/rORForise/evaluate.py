import collections

import check_pred as cp
import csv
import os, glob
import gzip
import math



nucleotides = ['A','C','G','T']

def make_kmers(k):
    if k == 1:
        return nucleotides
    return [n+c for n in nucleotides for c in make_kmers(k-1)]

codons = make_kmers(3)

# What is the expected number of occurrences of this codon given the GC probability
# in the genome and the total number of codons?
def expect(codon, total, GC_prob):
    prob = math.prod([GC_prob/2 if (c == 'G' or c=='C') else (1-GC_prob)/2 for c in codon])
    return prob*total

def print_with_expect(answer_type, codon_counts, GC_prob):
    print(answer_type+":")
    chosen = codon_counts[cp.answers[answer_type]]
    total = sum(chosen.values())
    for codon in codons: 
        print(f'{codon}\t{chosen[codon]:4d}\t{expect(codon, total, GC_prob):4.0f}')

def print_without_expect(answer_type, codon_counts):
    print(answer_type+":")
    chosen = codon_counts[cp.answers[answer_type]]
    for codon in codons: 
        print(f'{codon}\t{chosen[codon]:4d}')


        
def evaluate(genome_name, GC_prob, preds, intersect_filename, method):

    number_of_CDS_mappings_with_predictions = 0
    number_of_on_target_preds = 0

    reads_without_predictions = set()
    seen_read_names = set()
    good_overlap_read_names = set()

    answer_counts = collections.defaultdict(int)
    codon_counts = {}

    with gzip.open(intersect_filename, 'rt') as f:
        csvr = csv.reader(f, delimiter="\t") 
        header = next(csvr) # ignore single header line

        for bed_row in csvr:
            if bed_row[6] == "CDS":
                read_start = int(bed_row[2]) 
                read_end   = int(bed_row[3])
                read_name  = bed_row[1]
                read_dir   = bed_row[4]
                seen_read_names.add(read_name)
                cds_start = int(bed_row[7])
                cds_end = int(bed_row[8])

                if min(read_end,cds_end) - max(read_start,cds_start) + 1 < 60:
                    # the read does not overlap the CDS by enough for a pred
                    answer_counts[cp.answers["not enough read-CDS overlap"]] += 1
                    if read_name in preds:
                        pass # no more to add if read doesn't overlap CDS
                    else:
                        reads_without_predictions.add(read_name)
                else:
                    # the read does overlap the CDS by enough
                    good_overlap_read_names.add(read_name)
                    
                    if read_name in preds:
                        cds_dir = bed_row[9]
                        read_seq = bed_row[10]
                        number_of_CDS_mappings_with_predictions += 1
                        for (pred_start, pred_end, pred_dir) in preds[read_name]:

                            number_of_on_target_preds += 1

                            # answer_details is a set of pairs: (num, codon)
                            answer_details = cp.check_pred(cds_start, cds_end, cds_dir, read_start, read_end, read_dir, pred_start, pred_end, pred_dir,read_seq)

                            # Count the types of answers
                            for (answer, codon) in answer_details:
                                answer_counts[answer] += 1
                                if codon is not None: 
                                    if answer not in codon_counts:
                                        codon_counts[answer] = collections.defaultdict(int)
                                    codon_counts[answer][codon] += 1

                    else:
                        reads_without_predictions.add(read_name)

    print("Number of CDS-aligned reads", len(seen_read_names), sep = "\t")
    print("Number of CDS-aligned reads without predictions", len(reads_without_predictions), sep = "\t")
    print("Number of reads that overlap CDS enough", len(good_overlap_read_names), sep = "\t")
    print("Number of reads that overlap CDS enough with at least one pred", number_of_CDS_mappings_with_predictions, sep = "\t")
    print("Number of predictions", len([value for value in preds.values()]), sep = "\t")
    print("Number of preds on reads that overlap CDS enough", number_of_on_target_preds, sep = "\t")
    
    # Print the counted answers
    for (answer,count) in sorted(answer_counts.items()):
        print(f'{cp.inverse_answers[answer]}\t{count}')

    print_without_expect("correct stop", codon_counts)
    print_without_expect("alternative stop", codon_counts)
    print_with_expect("incorrect stop", codon_counts, GC_prob)
    print_without_expect("correct start", codon_counts)
    print_without_expect("alternative start", codon_counts)
    print_with_expect("incorrect start", codon_counts, GC_prob)




# This returns a dictionary of lists. Key is read name. List contains all preds that made for this read. Each pred has (start, stop, dir).
def read_preds(directory, genome_name, method, fragmentation_type, group):
    total_preds = 0
    preds = collections.defaultdict(list)
    pred_length_distn = collections.defaultdict(int) # length->frequency

    directory_path = os.path.join(directory, genome_name, method)
    file_pattern = f"*_{fragmentation_type}_{group}.gff.gz"
    gff_name = glob.glob(os.path.join(directory_path, file_pattern))
    print(directory_path, file_pattern, gff_name)
    print(gff_name[0])
    
    with gzip.open(gff_name[0], 'rt', encoding='utf-8')  as f:
        csvr = csv.reader(f, delimiter="\t")
        for row in csvr:
            if row[0].startswith("#"): # ignore the rows that are comments
                continue
            else:
                # do we need to check if row[2] == "CDS"?
                read_name = row[0].replace('@','') #.replace('ID=', '').split(';')[0] # This is not great
                #prediction_name = row[0].replace('ID=', '').split(';')[0].replace('@','')
                pred_start = int(row[3])
                pred_end = int(row[4])
                pred_dir = row[6]
                pred_length_distn[pred_end - pred_start+1] += 1
                preds[read_name].append((pred_start, pred_end, pred_dir))
                total_preds += 1

    # Return the count and the dictionary

    print("Prediction length distribution")
    for l in sorted(pred_length_distn.keys()):
        print(l, pred_length_distn[l])
    print()
    
    return (total_preds, preds)
                
    
