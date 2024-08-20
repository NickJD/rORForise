import collections

import check_pred as cp
import csv
import os, glob
import gzip
import math

dir = "Genome_Processing"
#genomes = ["Mycoplasma_genitalium_G37"] #,"Staphylococcus_aureus_502A"] 
genomes = ["Staphylococcus_aureus_502A"]
fragmentation_types = ["ART_errFree"]
subgroups = ['Combined']
methods = ["FragGeneScan"] #,"Pyrodigal","Naive-StORF-V1"]#,"Naive-StORF-V2","Naive-StORF-V3"] 

#methods = ["Pyrodigal"]
#methods = ["Naive-StORF-V1"]
# Hardcoding Myco/Mycoplasma for now


#cor = {"Pyrodigal":[],"FragGeneScan":[],"Naive-StORF-V1":[]}
#incor = {"Pyrodigal":[],"FragGeneScan":[],"Naive-StORF-V1":[]}
#alt = {"Pyrodigal":[],"FragGeneScan":[],"Naive-StORF-V1":[]}

#GC_prob = 0.3169 # Myco
GC_prob = 0.3292 # Staph


def expect(codon, total):
    prob = math.prod([GC_prob/2 if (c == 'G' or c=='C') else (1-GC_prob)/2 for c in codon])
    return prob*total



def evaluate(genome_name, preds, intersect_filename, method):

    number_of_CDS_mappings_with_predictions = 0
    number_of_on_target_preds = 0

    reads_without_predictions = []
    seen_read_names = []

    answer_counts = collections.defaultdict(int)
    incorrect_stop_codons = collections.defaultdict(int)
    incorrect_stop_code = cp.answers["incorrect stop"]

    with gzip.open(intersect_filename, 'rt') as f:
        csvr = csv.reader(f, delimiter="\t") 
        header = next(csvr) # ignore single header line

        for bed_row in csvr:
            if bed_row[6] == "CDS":
                read_start = int(bed_row[2]) 
                read_end   = int(bed_row[3])
                read_name  = bed_row[1]
                read_dir   = bed_row[4]

                if read_name in preds:
                    cds_start = int(bed_row[7])
                    cds_end = int(bed_row[8])
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
                            #incor[method].append([read_name,pred_start,pred_end])
                            if answer == incorrect_stop_code and codon is not None:
                                incorrect_stop_codons[codon] += 1
                else:
                    reads_without_predictions.append(read_name)

    set_a = set(seen_read_names)
    set_b = set(reads_without_predictions)
    unique_items_in_b =  set_b - set_a
    count_reads_without_predictions = len(unique_items_in_b)

    #number_of_predictions = len([value for value in preds.values()])

    print("Number of CDS mappings with at least one pred:", number_of_CDS_mappings_with_predictions)
    print("Number of on-target preds:", number_of_on_target_preds)
    # Print the counted answers
    for (answer,count) in sorted(answer_counts.items()):
        print(f'{cp.inverse_answers[answer]}: {count}')

    print("incorrect stop codons:")
    print(" codon\tcount\texpected")
    for (codon,count) in sorted(incorrect_stop_codons.items()):
        print(f' {codon}\t{count:4d}\t{expect(codon, sum(incorrect_stop_codons.values())):4.0f}')


                        

# This returns a dictionary of lists. Key is read name. List contains all preds that made for this read. Each pred has (start, stop, dir).
def read_preds(genome_name, method, fragmentation_type, group):
    total_preds = 0
    preds = collections.defaultdict(list)
    directory_path = os.path.join(dir, genome_name, method)
    file_pattern = f"*_{fragmentation_type}_{group}.gff.gz"
    gff_name = glob.glob(os.path.join(directory_path, file_pattern))
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
                preds[read_name].append((pred_start, pred_end, pred_dir))
                total_preds += 1

    # Return the count and the dictionary
    return (total_preds, preds)
                
    
def main():
    for genome_name in genomes:
        print(genome_name)
        datadir = dir + "/" + genome_name + "/Processing"
        intersect_bed_filename = datadir + "/"+genome_name+"_Reads_Intersect.tsv.gz"

        for method in methods:
            print(method)

            for fragmentation_type in fragmentation_types:
                print(fragmentation_type)

                for group in subgroups:
                    (total_preds, preds) = read_preds(genome_name, method, fragmentation_type, group)
                    print("Number of predictions made for reads:", total_preds)
                    evaluate(genome_name, preds, intersect_bed_filename, method)

        import sys
        sys.exit()
        for correct in cor["FragGeneScan"]:
            for incorrect in incor["Naive-StORF-V1"]:
                if correct[0] in incorrect[0]:
                    print()
                    #print(incorrect)
            for incorrect in alt["Naive-StORF-V1"]:
                if correct[0] in incorrect[0]:
                    print()
                    #print(incorrect)
                    

print("")


                    
main()
# use as follows: python evaluate.py
