import argparse
import shutil
import os

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

    parser.add_argument('-o', '--output_dir', type=str, required=True,
                        help="Directory to store output files.")

    parser.add_argument('-gc_prob', type=float, required=True,
                        help="GC probability of the genome being processed (e.g., 0.3169 for Mycoplasma genitalium).")

    parser.add_argument('-l', '--overlap_threshold', type=int, default=60, help="Minimum number of bases of overlap that are required for the read to overlap the CDS by before a prediction is inspected.")
                        
    return parser.parse_args()


def main():
    options = parse_args()

    # Ensure output directory exists and is empty
    if os.path.exists(options.output_dir):
        shutil.rmtree(options.output_dir)
    os.makedirs(options.output_dir)

    total_preds, preds = evaluate.read_preds(options.predictions_gff, options.output_dir)
    print("Number of predictions made for reads", total_preds, sep="\t")

    intersect_bed_filename = options.intersect_bed[0]

    evaluate.evaluate(intersect_bed_filename, preds, options.gc_prob, options.output_dir, options.overlap_threshold)


if __name__ == "__main__":
    print("Running rORForise")
    main()
