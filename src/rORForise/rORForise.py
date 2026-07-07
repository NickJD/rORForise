import argparse
import shutil
from pathlib import Path

try:
    from .utils import *
    from .constants import *
    from .evaluate import read_preds, evaluate
    from .report import build_single_tool_report
except (ImportError, ModuleNotFoundError):
    from utils import *
    from constants import *
    from evaluate import read_preds, evaluate
    from report import build_single_tool_report



def parse_args():
    parser = argparse.ArgumentParser(description="rORForise version: " + rORForise_VERSION + " - Process genome evaluation parameters.")

    parser.add_argument('-int_bed', '--intersect_bed', type=str, required=True,
                        help="File path to intersect bedfile with read to CDS mappings.")

    parser.add_argument('-p_gff', '--predictions_gff', type=str, required=True,
                        help="Path to the predictions in GFF format.")

    parser.add_argument('-o', '--output_dir', type=str, required=True,
                        help="Directory to store output files.")

    parser.add_argument('-gc_prob', type=float, required=True,
                        help="GC probability of the genome being processed (e.g., 0.3169 for Mycoplasma genitalium).")

    parser.add_argument('-prefix', '--output_prefix', type=str, default="orf_evaluation", help="Prefix to add to output files. Default is orf_evaluation.") 

    parser.add_argument('-l', '--overlap_threshold', type=int, default=60, help="Minimum number of bases of overlap that are required for the read to overlap the CDS by before a prediction is inspected.")

    parser.add_argument('--verbose', action='store_true', help='Verbose output')
    parser.add_argument('--clear-output', action='store_true',
                        help='Delete the output directory before writing results. By default existing files are preserved unless overwritten by this run.')
    parser.add_argument('--report', choices=['none', 'html'], default='none',
                        help='Generate an optional visual report. Use html to write report.html and figure files.')
    parser.add_argument('--report-tool-name', default='rORForise',
                        help='Tool name to display in the optional report.')
    parser.add_argument('--max-tracks', type=int, default=12,
                        help='Maximum read-level track examples in the optional HTML report.')

    return parser.parse_args()


def main():
    options = parse_args()

    if options.clear_output and os.path.exists(options.output_dir):
        shutil.rmtree(options.output_dir)
    os.makedirs(options.output_dir, exist_ok=True)

    total_preds, preds = read_preds(options.predictions_gff, options.output_dir, options.output_prefix, verbose=options.verbose)
    evaluate(options.intersect_bed, preds, options.gc_prob, options.output_dir, options.output_prefix, options.overlap_threshold, verbose=options.verbose)
    if options.report == 'html':
        report_path = build_single_tool_report(
            Path(options.output_dir),
            options.output_prefix,
            options.report_tool_name,
            options.intersect_bed,
            preds,
            overlap_threshold=options.overlap_threshold,
            max_tracks=options.max_tracks,
        )
        print(f"HTML report written to: {report_path}")


if __name__ == "__main__":
    main()
