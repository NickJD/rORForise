import argparse
import evaluate


def parse_args():
    parser = argparse.ArgumentParser(description="Process genome evaluation parameters.")

    parser.add_argument('-d', '--directory', type=str, required=True,
                        help="Path to the main directory containing genome processing subdirectories.")

    parser.add_argument('-g', '--genomes', type=str, nargs='+', required=True,
                        help="List of genome names to process. Example: 'Mycoplasma_genitalium_G37 Staphylococcus_aureus_502A'")

    parser.add_argument('--gc_prob', type=float, required=True,
                        help="GC probability of the genome being processed (e.g., 0.3169 for Mycoplasma genitalium).")

    parser.add_argument('-f', '--fragmentation_types', type=str, nargs='+', required=True,
                        help="List of fragmentation types to use (e.g., 'ART_errFree').")

    parser.add_argument('-s', '--subgroups', type=str, nargs='+', required=True,
                        help="List of subgroups for evaluation (e.g., 'Combined').")

    parser.add_argument('-m', '--methods', type=str, nargs='+', required=True,
                        help="List of methods for evaluation (e.g., 'FragGeneScan').")

    return parser.parse_args()


def main():
    args = parse_args()

    for genome_name in args.genomes:
        print(genome_name)
        datadir = f"{args.directory}/{genome_name}/Processing"
        intersect_bed_filename = f"{datadir}/{genome_name}_Reads_Intersect_CDS.tsv.gz"

        for method in args.methods:
            print(method)

            for fragmentation_type in args.fragmentation_types:
                print(fragmentation_type)

                for group in args.subgroups:
                    total_preds, preds = evaluate.read_preds(args.directory, genome_name, method, fragmentation_type,
                                                             group)
                    print("Number of predictions made for reads", total_preds, sep="\t")
                    evaluate.evaluate(genome_name, args.gc_prob, preds, intersect_bed_filename, method)


if __name__ == "__main__":
    print("Running rORForise")
    main()
