import evaluate

directory = "../../Genome_Processing"

genomes = ["Mycoplasma_genitalium_G37"] 
#genomes = ["Staphylococcus_aureus_502A"]
#genomes = ["Escherichia_coli_k_12"]

GC_prob = 0.3169 # Myco
#GC_prob = 0.3292 # Staph
#GC_prob = 0.5080 # Ecoli

fragmentation_types = ["ART_errFree"]
subgroups = ['Combined']
#subgroups = ['Combined_FGS_sanger_5']

methods = ["FragGeneScan"] #,"Pyrodigal","Naive-StORF-V1"]#,"Naive-StORF-V2","Naive-StORF-V3"] 
#methods = ["Pyrodigal"]
#methods = ["Naive-StORF-V1"]

def main():
    for genome_name in genomes:
        print(genome_name)
        datadir = directory + "/" + genome_name + "/Processing"
        intersect_bed_filename = datadir + "/"+genome_name+"_Reads_Intersect_CDS.tsv.gz"

        for method in methods:
            print(method)

            for fragmentation_type in fragmentation_types:
                print(fragmentation_type)

                for group in subgroups:
                    (total_preds, preds) = evaluate.read_preds(directory, genome_name, method, fragmentation_type, group)
                    print("Number of predictions made for reads", total_preds, sep = "\t")
                    evaluate.evaluate(genome_name, GC_prob, preds, intersect_bed_filename, method)


                    
main()
# use as follows: python rORForise_Testing.py
