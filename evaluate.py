import collections

import check_pred as cp
import csv
import os, glob
import gzip


dir = "Genome_Processing"
genomes = ["Mycoplasma_genitalium_G37"] #,"Staph"]
fragmentation_types = ["ART_errFree"]
subgroups = ['Combined']
methods = ["FragGeneScan","Pyrodigal","Naive-StORF-V1"]#,"Naive-StORF-V2","Naive-StORF-V3"] #,"Pyrodigal","FragGeneScan","FrameRate"]

# Hardcoding Myco/Mycoplasma for now

""" bedfile:
1 Chromosome
2 848
3 1032
4 Chromosome-41280/1
5 42
6 +
7 848
8 1032
9 0,0,0
10 1
11 184,
12 0,
13 Chromosome
14 ena
15 CDS
16 686
17 1828
18 .
19 +
20 0
21 ID=CDS:AAC71217;Parent=transcript:AAC71217;protein_id=AAC71217  184
"""

def evaluate(genome_name, preds, bed_intersect_filename, method):

    correct_starts, incorrect_starts, alternative_starts, middle_alternative_starts, correct_stops, \
        incorrect_stops, alternative_stops, middle_alternative_stops,  reads_with_predictions, \
        reads_without_predictions, correct_frames, incorrect_frames, correct_directions, incorrect_directions, \
        prediction_ends_before_cds_starts, prediction_starts_before_cds_ends = 0,0, 0, 0, 0, 0 , 0, 0 , 0, 0, 0,0,0,0,0,0

    with gzip.open(bed_intersect_filename, 'rt') as f:
        csvr = csv.reader(f, delimiter="\t")
        for bed_row in csvr:
            if bed_row[14] == "CDS":

                read_start = int(bed_row[1])
                read_end   = int(bed_row[2])
                read_name  = bed_row[3]
                read_dir   = bed_row[5]



                if read_name in preds: # Not counting the number of preds but the number of reads with at least one pred
                    reads_with_predictions += 1 #len(preds[read_name])
                    cds_start = int(bed_row[15])
                    cds_end   = int(bed_row[16])
                    cds_dir   = bed_row[18]

                    if read_name == 'Chromosome-77340/1':
                        print("answer")

                    for (pred_start, pred_end, pred_dir) in preds[read_name]:
                        #print(cds_start, cds_end, cds_dir, read_start, read_end, read_dir, pred_start, pred_end, pred_dir)
                        answer = cp.check_pred(cds_start, cds_end, cds_dir, read_start, read_end, read_dir, pred_start, pred_end, pred_dir)

                        if read_name == 'Chromosome-77340/1':
                            print(answer)

                        if 0 in answer:
                            correct_starts +=1
                            # if 'Naive' in method:
                            #     print(read_name)
                        elif 1 in answer:
                            alternative_starts +=1
                        elif 2 in answer:
                            middle_alternative_starts +=1
                        elif 3 in answer:
                            incorrect_starts +=1
                        if 4 in answer:
                            correct_stops +=1
                        elif 5 in answer:
                            alternative_stops +=1
                        elif 6 in answer:
                            middle_alternative_stops +=1
                        elif 7 in answer:
                            incorrect_stops +=1
                        elif 8 in answer:
                            correct_frames +=1
                        elif 9 in answer:
                            incorrect_frames +=1
                        elif 10 in answer:
                            correct_directions +=1
                        elif 11 in answer:
                            incorrect_directions +=1
                        elif 12 in answer:
                            prediction_ends_before_cds_starts +=1
                        elif 13 in answer:
                            prediction_starts_before_cds_ends +=1





#                        print(read_name,"CDS Start/End: " + str(cds_start) + "/" + str(cds_end) + " Read Start/End: " + str(read_start) + "/" +
#                              str(read_end),preds[read_name],answer)
                else:
                    reads_without_predictions +=1
#                    print(read_name)
    print("reads_with_predictions: " + str(reads_with_predictions))
    print("reads_without_predictions: " + str(reads_without_predictions))
    print("correct_starts: " + str(correct_starts))
    print("alternative_starts: " + str(alternative_starts))
    print("middle_alternative_starts: " + str(middle_alternative_starts))
    print("incorrect_starts: " + str(incorrect_starts))
    print("correct_stops: " + str(correct_stops))
    print("alternative_stops: " + str(alternative_stops))
    print("middle_alternative_stops: " + str(middle_alternative_stops))
    print("incorrect_stops: " + str(incorrect_stops))
    print("correct_frames: " + str(correct_frames))
    print("incorrect_frames: " + str(incorrect_frames))
    print("correct_directions: " + str(correct_directions))
    print("incorrect_directions: " + str(incorrect_directions))
    print("prediction_ends_before_cds_starts: " + str(prediction_ends_before_cds_starts))
    print("prediction_starts_before_cds_ends: " + str(prediction_starts_before_cds_ends))
                        

# This returns a dictionary of lists. Key is read name. List contains all preds that made for this read. Each pred has (start, stop, dir).
def read_preds(genome_name, method, fragmentation_type, group):
    preds = collections.defaultdict(list)
    directory_path = os.path.join(dir, genome_name, method)
    file_pattern = f"*_{fragmentation_type}_{group}.gff"
    gff_name = glob.glob(os.path.join(directory_path, file_pattern))
    print(gff_name[0])
    with open(gff_name[0]) as f:
        csvr = csv.reader(f, delimiter="\t")
        for row in csvr:
            if row[0].startswith("#"):
                continue
            else:
                # do we need to check if row[2] == "CDS"?
                read_name = row[0].replace('@','') #.replace('ID=', '').split(';')[0] # This is not great
                #prediction_name = row[0].replace('ID=', '').split(';')[0].replace('@','')
                pred_start = int(row[3])
                pred_end = int(row[4])
                pred_dir = row[6]
                preds[read_name].append((pred_start, pred_end, pred_dir))

    #print(preds)
    return preds
                
    
def main():
    for genome_name in genomes:
        print(genome_name)
        datadir = dir + "/" + genome_name + "/Processing"
        intersect_bed_filename = datadir + "/Myco_Reads_Intersect.bed.gz"

        for method in methods:
            print(method)

            for fragmentation_type in fragmentation_types:
                print(fragmentation_type)

                for group in subgroups:
                    preds = read_preds(genome_name, method, fragmentation_type, group)
                    evaluate(genome_name, preds, intersect_bed_filename, method)

        
main()

# use as follows: python evaluate.py
