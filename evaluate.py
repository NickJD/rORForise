import check_pred as cp
import csv
import os, glob


dir = "Genome_Processing"
genomes = ["Mycoplasma_genitalium_G37"] #,"Staph"]
fragmentation_types = ["ART_errFree_paired"]
methods = ["FragGeneScan","Pyrodigal","Naive-StORF-V1","Naive-StORF-V2","Naive-StORF-V3"] #,"Pyrodigal","FragGeneScan","FrameRate"]

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

def evaluate(genome_name, preds, bed_intersect_filename):
    with open(bed_intersect_filename) as f:
        csvr = csv.reader(f, delimiter="\t")
        for bed_row in csvr:
            if bed_row[14] == "CDS":

                read_start = int(bed_row[1])
                read_end   = int(bed_row[2])
                read_name  = bed_row[3]
                read_dir   = bed_row[5]

                if read_name in preds:
                    cds_start = int(bed_row[15])
                    cds_end   = int(bed_row[16])
                    cds_dir   = bed_row[18]

                    for (pred_start, pred_end, pred_dir) in preds[read_name]:
                        #print(cds_start, cds_end, cds_dir, read_start, read_end, read_dir, pred_start, pred_end, pred_dir)
                        answer = cp.check_pred(cds_start, cds_end, cds_dir, read_start, read_end, read_dir, pred_start, pred_end, pred_dir)
                        print(read_name,"CDS Start/End: " + str(cds_start) + "/" + str(cds_end) + " Read Start/End: " + str(read_start) + "/" +
                              str(read_end),preds[read_name],answer)
                else:
                    print(read_name)
                        

# This returns a dictionary of lists. Key is read name. List contains all preds that made for this read. Each pred has (start, stop, dir).
def read_preds(genome_name, method, fragmentation_type):
    preds = {}
    directory_path = os.path.join(dir, genome_name, method)
    file_pattern = f"*_{fragmentation_type}.gff"
    gff_name = glob.glob(os.path.join(directory_path, file_pattern))
    print(gff_name[0])
    with open(gff_name[0]) as f:
        csvr = csv.reader(f, delimiter="\t")
        for row in csvr:
            if row[0].startswith("#"):
                continue
            else:
                # do we need to check if row[2] == "CDS"?
                read_name = row[0]
                if read_name.startswith('@'):
                    read_name  = row[0][1:] # remove the leading '@' char

                pred_start = int(row[3])
                pred_end   = int(row[4])
                pred_dir   = row[6]
                preds[read_name] = preds.get(read_name, []) + [(pred_start,pred_end,pred_dir)]
    #print(preds)
    return preds
                
    
def main():
    for genome_name in genomes:
        print(genome_name)
        datadir = dir + "/" + genome_name + "/Processing"
        intersect_bed_filename = datadir + "/Myco_Reads_Intersect.bed"

        for method in methods:
            print(method)

            for fragmentation_type in fragmentation_types:
                print(fragmentation_type)
                preds = read_preds(genome_name, method, fragmentation_type)
                evaluate(genome_name, preds, intersect_bed_filename)

        
main()

# use as follows: python evaluate.py
