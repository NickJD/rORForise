import argparse
#import evaluate

try: # Try to import from the package if available
    from .constants import *
    #from .evaluate import evaluate as evaluate
    from .evaluate import *
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from constants import *
    #from evaluate import evaluate as evaluate
    import evaluate

def parse_args():
    parser = argparse.ArgumentParser(description="rORForise version: " + rORForise_VERSION + " - Process genome evaluation parameters.")

    parser.add_argument('-int_bed', '--intersect_bed', type=str, nargs='+', required=True,
                        help="File path to intersect bedfile with read to CDS mappings.")

    parser.add_argument('-p_gff', '--predictions_gff', type=str, required=True,
                        help="Path to the predictions in GFF format.")

    parser.add_argument('-gc_prob', type=float, required=True,
                        help="GC probability of the genome being processed (e.g., 0.3169 for Mycoplasma genitalium).")


    return parser.parse_args()


def main():
    options = parse_args()

    total_preds, preds = evaluate.read_preds(options.predictions_gff)
    print("Number of predictions made for reads", total_preds, sep="\t")
    #evaluate.evaluate(genome_name, args.gc_prob, preds, intersect_bed_filename, method)

    intersect_bed_filename = options.intersect_bed[0]

    evaluate.evaluate(intersect_bed_filename, preds, options.gc_prob)


if __name__ == "__main__":
    print("Running rORForise")
    main()
