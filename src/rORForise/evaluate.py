import collections
import check_pred as cp
import csv
import os
import gzip
import math
from datetime import datetime

nucleotides = ['A', 'C', 'G', 'T']


def make_kmers(k):
    if k == 1:
        return nucleotides
    return [n + c for n in nucleotides for c in make_kmers(k - 1)]


codons = make_kmers(3)


# What is the expected number of occurrences of this codon given the GC probability
# in the genome and the total number of codons?
def expect(codon, total, GC_prob):
    prob = math.prod([GC_prob / 2 if (c == 'G' or c == 'C') else (1 - GC_prob) / 2 for c in codon])
    return prob * total


def write_codon_analysis_csv(codon_counts, gc_prob, output_dir, output_prefix):

    # Get the codon-related answer types that we have data for
    codon_answer_types = {
        cp.answers['correct stop']: 'correct stop',
        cp.answers['alternative stop']: 'alternative stop',
        cp.answers['incorrect stop']: 'incorrect stop',
        cp.answers['correct start']: 'correct start',
        cp.answers['alternative start']: 'alternative start',
        cp.answers['incorrect start']: 'incorrect start',
        cp.answers['middle or alternative start']: 'middle or alternative start',
        cp.answers['middle or alternative stop']: 'middle or alternative stop',
        cp.answers['middle incorrect start']: 'middle incorrect start',
        cp.answers['middle incorrect stop']: 'middle incorrect stop'
    }

    # Only write codon analysis if we have codon data
    if not codon_counts:
        print("No codon data available for analysis")
        return

    # Summary CSV with all answer types
    summary_file = os.path.join(output_dir, f"{output_prefix}_codon_summary.csv")
    with open(summary_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['answer_type', 'codon', 'observed', 'expected', 'obs_exp_ratio'])

        for answer_code, answer_name in codon_answer_types.items():
            if answer_code in codon_counts:
                chosen = codon_counts[answer_code]
                total = sum(chosen.values())

                for codon in codons:
                    observed = chosen[codon]
                    if answer_name in ['incorrect stop', 'incorrect start', 'middle incorrect start',
                                       'middle incorrect stop']:
                        expected = expect(codon, total, gc_prob)
                        ratio = observed / expected if expected > 0 else float('inf')
                    else:
                        expected = 'N/A'
                        ratio = 'N/A'

                    writer.writerow([answer_name, codon, observed, expected, ratio])

    # Individual detailed CSV files for each answer type
    for answer_code, answer_name in codon_answer_types.items():
        if answer_code in codon_counts:
            filename = os.path.join(output_dir, f"{output_prefix}_codons_{answer_name.replace(' ', '_')}.csv")
            chosen = codon_counts[answer_code]
            total = sum(chosen.values())

            with open(filename, 'w', newline='') as f:
                writer = csv.writer(f)
                if answer_name in ['incorrect stop', 'incorrect start', 'middle incorrect start',
                                   'middle incorrect stop']:
                    writer.writerow(['codon', 'observed', 'expected', 'obs_exp_ratio', 'enrichment'])
                else:
                    writer.writerow(['codon', 'observed', 'frequency'])

                for codon in codons:
                    observed = chosen[codon]
                    frequency = observed / total if total > 0 else 0

                    if answer_name in ['incorrect stop', 'incorrect start', 'middle incorrect start',
                                       'middle incorrect stop']:
                        expected = expect(codon, total, gc_prob)
                        ratio = observed / expected if expected > 0 else float('inf')
                        enrichment = 'enriched' if ratio > 1.5 else 'depleted' if ratio < 0.5 else 'normal'
                        writer.writerow([codon, observed, f"{expected:.1f}", f"{ratio:.2f}", enrichment])
                    else:
                        writer.writerow([codon, observed, f"{frequency:.3f}"])


def write_summary_csv(stats, output_dir, output_prefix):

    summary_file = os.path.join(output_dir, f"{output_prefix}_run_summary.csv")


    with open(summary_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['metric', 'value'])

        # Basic statistics
        writer.writerow(['timestamp', datetime.now().isoformat()])
        writer.writerow(['total_number_of_reads', stats['total_number_of_reads']])
        writer.writerow(['total_cds_aligned_reads', stats['total_cds_aligned_reads']])
        writer.writerow(['reads_without_predictions', stats['reads_without_predictions']])
        writer.writerow(['reads_good_overlap', stats['reads_good_overlap']])
        writer.writerow(['reads_good_overlap_with_preds', stats['reads_good_overlap_with_preds']])
        writer.writerow(['total_predictions', stats['total_predictions']])
        writer.writerow(['on_target_predictions', stats['on_target_predictions']])

        # Coverage metrics
        if stats['reads_good_overlap'] > 0:
            writer.writerow(['prediction_coverage',
                             f"{stats['reads_good_overlap_with_preds'] / stats['reads_good_overlap']:.3f}"])
        writer.writerow(['avg_preds_per_read',
                         f"{stats['total_predictions'] / stats['total_number_of_reads']:.2f}"]) # Currently total reads which might not be correct




def write_prediction_length_csv(pred_length_distn, output_dir, output_prefix):
    length_file = os.path.join(output_dir, f"{output_prefix}_prediction_lengths.csv")

    with open(length_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['length', 'count', 'frequency'])

        total = sum(pred_length_distn.values())

        for length in sorted(pred_length_distn.keys()):
            count = pred_length_distn[length]
            frequency = count / total if total > 0 else 0
            writer.writerow([length, count, f"{frequency:.4f}"])

def write_track_predictions(track_preds, output_dir, output_prefix):
    filename = os.path.join(output_dir, f"{output_prefix}_track_predictions.csv")
    with open(filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['read_name', 'pred_start', 'pred_end', 'answer_code', 'answer_name', 'codon'])
        #pred_counter = 0
        for read_name in sorted(track_preds.keys()):
            for pred_start, pred_end, answer_details in track_preds[read_name]:
                #pred_counter += 1
                for answer_code, codon in answer_details:
                    answer_name = cp.inverse_answers.get(answer_code, str(answer_code))

                    writer.writerow([read_name, pred_start, pred_end, answer_code, answer_name, codon if codon is not None else ''])


def write_detailed_results_csv(answer_counts, output_dir, output_prefix, stats, context_counts):
    """
    Write detailed results with percentages calculated relative to the appropriate context.

    Args:
        answer_counts: Dictionary of answer_code -> count
        output_dir: Directory for output files
        output_prefix: Prefix for output filenames
        stats: Dictionary of general statistics
        context_counts: Dictionary tracking how many predictions overlapped start/stop/middle
    """
    results_file = os.path.join(output_dir, f"{output_prefix}_detailed_results.csv")

    with open(results_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['answer_type', 'count', 'percentage', 'denominator', 'context'])

        # Define which answer types belong to which context
        # Start-related answers (only for predictions overlapping CDS start)
        start_only_answer_types = [
            cp.answers['correct start'],
            cp.answers['incorrect start'],
            cp.answers['alternative start']
        ]

        # Stop-related answers (only for predictions overlapping CDS stop)
        stop_only_answer_types = [
            cp.answers['correct stop'],
            cp.answers['incorrect stop'],
            cp.answers['alternative stop']
        ]

        # Middle start-related (only for middle predictions)
        middle_start_answer_types = [
            cp.answers['middle or alternative start'],
            cp.answers['middle incorrect start']
        ]

        # Middle stop-related (only for middle predictions)
        middle_stop_answer_types = [
            cp.answers['middle or alternative stop'],
            cp.answers['middle incorrect stop']
        ]

        # Frame answers (all on-target predictions can have these)
        frame_answer_types = [
            cp.answers['correct frame'],
            cp.answers['incorrect frame']
        ]

        # Direction answers (all predictions including off-target)
        direction_answer_types = [
            cp.answers['correct direction'],
            cp.answers['incorrect direction']
        ]

        # Middle-only answers (only for predictions in middle of CDS)
        middle_only_answer_types = [
            cp.answers['middle']
        ]

        # Boundary issues (predictions that extend beyond CDS boundaries)
        boundary_answer_types = [
            cp.answers['prediction ends before cds starts'],
            cp.answers['prediction starts after cds ends']
        ]

        # Overlap issues (insufficient read-CDS overlap)
        overlap_answer_types = [
            cp.answers['not enough read-CDS overlap']
        ]

        # Create list of all possible answer codes  - but could be copied from check_pred?
        all_answer_codes = [
            cp.answers['correct start'],
            cp.answers['correct stop'],
            cp.answers['correct frame'],
            cp.answers['correct direction'],
            cp.answers['incorrect stop'],
            cp.answers['incorrect start'],
            cp.answers['incorrect frame'],
            cp.answers['incorrect direction'],
            cp.answers['prediction ends before cds starts'],
            cp.answers['prediction starts after cds ends'],
            cp.answers['alternative start'],
            cp.answers['middle or alternative start'],
            cp.answers['alternative stop'],
            cp.answers['middle or alternative stop'],
            cp.answers['middle'],
            cp.answers['middle incorrect start'],
            cp.answers['middle incorrect stop'],
            cp.answers['not enough read-CDS overlap']
        ]

        # Write each answer with its contextual percentage (including zeros)
        for answer_code in sorted(all_answer_codes):
            answer_name = cp.inverse_answers[answer_code]
            count = answer_counts.get(answer_code, 0)  # Use .get() to default to 0

            # Determine the appropriate denominator based on answer type
            if answer_code in start_only_answer_types:
                # Only predictions that overlap the CDS start can give these answers
                denominator = context_counts['overlaps_start']
                context = 'predictions overlapping CDS start'

            elif answer_code in stop_only_answer_types:
                # Only predictions that overlap the CDS stop can give these answers
                denominator = context_counts['overlaps_stop']
                context = 'predictions overlapping CDS stop'

            elif answer_code in middle_start_answer_types:
                # Only middle predictions (not overlapping start or stop) can give these
                denominator = context_counts['middle_only']
                context = 'predictions in CDS middle only'

            elif answer_code in middle_stop_answer_types:
                # Only middle predictions (not overlapping start or stop) can give these
                denominator = context_counts['middle_only']
                context = 'predictions in CDS middle only'

            elif answer_code in middle_only_answer_types:
                # Middle answer is only for predictions entirely within CDS
                denominator = context_counts['middle_only']
                context = 'predictions in CDS middle only'

            elif answer_code in frame_answer_types:
                # Frame can be assessed for all on-target predictions
                denominator = stats['on_target_predictions']
                context = 'all on-target predictions'

            elif answer_code in direction_answer_types:
                # Direction can be assessed for all predictions
                denominator = stats['on_target_predictions']
                context = 'all on-target predictions'

            elif answer_code in boundary_answer_types:
                # Boundary issues can occur with any on-target prediction
                denominator = stats['on_target_predictions']
                context = 'all on-target predictions'

            elif answer_code in overlap_answer_types:
                # Overlap issues are calculated differently - out of total reads
                denominator = stats['total_cds_aligned_reads']
                context = 'all CDS-aligned reads'

            else:
                # Fallback - shouldn't happen if all 18 types are covered
                denominator = stats['on_target_predictions']
                context = 'all on-target predictions (uncategorized)'
                print(f"Warning: Uncategorized answer type: {answer_name}")

            percentage = (count / denominator * 100) if denominator > 0 else 0
            writer.writerow([answer_name, count, f"{percentage:.2f}", denominator, context])



def evaluate(intersect_bed_filename, preds, gc_prob, output_dir, output_prefix="orf_evaluation"):
    number_of_CDS_mappings_with_predictions = 0
    number_of_on_target_preds = 0

    track_preds = collections.defaultdict(list)

    reads_without_predictions = set()
    seen_read_names = set()
    good_overlap_read_names = set()

    answer_counts = collections.defaultdict(int)
    codon_counts = {}

    # Track context: how many predictions overlapped start, stop, or were in middle
    context_counts = {
        'overlaps_start': 0,
        'overlaps_stop': 0,
        'middle_only': 0
    }

    #with gzip.open(intersect_bed_filename, 'rt') as f:
    with (gzip.open(intersect_bed_filename, 'rt', encoding='utf-8') if intersect_bed_filename.endswith('.gz') else open(intersect_bed_filename,
                                                                                                        'r', encoding='utf-8')) as f:

        csvr = csv.reader(f, delimiter="\t")
        header = next(csvr)  # ignore single header line
        read_num = 0
        for bed_row in csvr:
            read_num +=1
            if bed_row[6] == "CDS":
                read_start = int(bed_row[2])
                read_end = int(bed_row[3])
                read_name = bed_row[1]
                read_dir = bed_row[4]
                seen_read_names.add(read_name)
                cds_start = int(bed_row[7])
                cds_end = int(bed_row[8])

                if min(read_end, cds_end) - max(read_start, cds_start) + 1 < 60: # 60 could be a user-threshold
                    # the read does not overlap the CDS by enough for a pred
                    answer_counts[cp.answers["not enough read-CDS overlap"]] += 1
                    if read_name not in preds:
                        reads_without_predictions.add(read_name)
                else:
                    # the read does overlap the CDS by enough
                    good_overlap_read_names.add(read_name)

                    if read_name in preds:
                        cds_dir = bed_row[9]
                        read_seq = bed_row[10]
                        number_of_CDS_mappings_with_predictions += 1

                        # Process each prediction for this read
                        for prediction_name, (pred_start, pred_end, pred_dir) in preds[read_name].items():
                            number_of_on_target_preds += 1

                            # Determine if this prediction overlaps start, stop, or is in middle
                            # This logic needs to account for strand direction
                            if cds_dir == '+' and read_dir == '+':
                                read_captures_cds_start = read_start <= cds_start
                                read_captures_cds_end = read_end >= cds_end
                            elif cds_dir == '+' and read_dir == '-':
                                read_captures_cds_start = read_start <= cds_start
                                read_captures_cds_end = read_end >= cds_end
                            elif cds_dir == '-' and read_dir == '+':
                                read_captures_cds_start = read_end >= cds_end
                                read_captures_cds_end = read_start <= cds_start
                            elif cds_dir == '-' and read_dir == '-':
                                read_captures_cds_start = read_end >= cds_end
                                read_captures_cds_end = read_start <= cds_start
                            else:
                                read_captures_cds_start = False
                                read_captures_cds_end = False

                            # Count the context for this prediction
                            if read_captures_cds_start:
                                context_counts['overlaps_start'] += 1
                            if read_captures_cds_end:
                                context_counts['overlaps_stop'] += 1
                            if not read_captures_cds_start and not read_captures_cds_end:
                                context_counts['middle_only'] += 1

                            # Get answer_details for this prediction
                            answer_details = cp.check_pred(cds_start, cds_end, cds_dir, read_start, read_end, read_dir,
                                                           pred_start, pred_end, pred_dir, read_seq)
                            track_preds[prediction_name].append((pred_start, pred_end, answer_details))
                            # Count the types of answers
                            for (answer, codon) in answer_details:
                                answer_counts[answer] += 1
                                if codon is not None:
                                    if answer not in codon_counts:
                                        codon_counts[answer] = collections.defaultdict(int)
                                    codon_counts[answer][codon] += 1

                    else:
                        reads_without_predictions.add(prediction_name)

    # Collect statistics
    stats = {
        'total_number_of_reads': read_num,  # excluding header
        'total_cds_aligned_reads': len(seen_read_names),
        'reads_without_predictions': len(reads_without_predictions),
        'reads_good_overlap': len(good_overlap_read_names),
        'reads_good_overlap_with_preds': number_of_CDS_mappings_with_predictions,
        'total_predictions': sum(len(pred_list) for pred_list in preds.values()),
        'on_target_predictions': number_of_on_target_preds
    }

    # Write CSV outputs
    write_summary_csv(stats, output_dir, output_prefix)
    write_detailed_results_csv(answer_counts, output_dir, output_prefix, stats, context_counts)
    write_codon_analysis_csv(codon_counts, gc_prob, output_dir, output_prefix)
    write_track_predictions(track_preds, output_dir, output_prefix)



    # Print summary to stdout for immediate feedback
    print(f"Evaluation complete. CSV files written with prefix: {output_prefix}")
    print(
        f"Summary: {stats['reads_good_overlap_with_preds']}/{stats['reads_good_overlap']} reads with good overlap have predictions")
    print(
        f"Context: {context_counts['overlaps_start']} predictions overlap start, {context_counts['overlaps_stop']} overlap stop, {context_counts['middle_only']} are middle-only")
    print(
        f"Prediction accuracy: {answer_counts.get(cp.answers.get('correct stop', 0), 0) + answer_counts.get(cp.answers.get('correct start', 0), 0)} correct predictions")

    print(f"\nAnalysis complete! Generated CSV files:")
    print("Output Directory: " + str(os.path.join(output_dir)))
    print(f"  - {output_prefix}_run_summary.csv")
    print(f"  - {output_prefix}_detailed_results.csv")
    print(f"  - {output_prefix}_codon_summary.csv")
    print(f"  - {output_prefix}_prediction_lengths.csv")
    print(f"Individual codon analysis files: {output_prefix}_codons_*.csv")

    return stats, answer_counts, codon_counts, context_counts




def read_preds(predictions_gff, output_dir, output_prefix="orf_evaluation"):
    total_preds = 0
    preds = collections.defaultdict(dict)
    pred_length_distn = collections.defaultdict(int)  # length->frequency

    print("Predictions GFF: " + predictions_gff)

    with (gzip.open(predictions_gff, 'rt', encoding='utf-8') if predictions_gff.endswith('.gz') else open(predictions_gff,
                                                                                                    'r', encoding='utf-8')) as f:
        csvr = csv.reader(f, delimiter="\t")
        for row in csvr:
            if row[0].startswith("#"):  # ignore the rows that are comments
                continue
            else:
                # do we need to check if row[2] == "CDS"?
                read_name = row[0].replace('@', '')  # .replace('ID=', '').split(';')[0] # This is not great
                prediction_name = row[8].replace('ID=', '').split(';')[0].replace('@', '')
                pred_start = int(row[3])
                pred_end = int(row[4])
                pred_dir = row[6]
                pred_length_distn[pred_end - pred_start + 1] += 1
                preds[read_name][prediction_name] = (pred_start, pred_end, pred_dir)
                total_preds += 1

    # Write prediction length distribution to CSV
    write_prediction_length_csv(pred_length_distn, output_dir, output_prefix)

    # Print brief summary to stdout
    print(f"Loaded {total_preds} predictions from {len(preds)} reads")
    print(f"Prediction length distribution written to {output_prefix}_prediction_lengths.csv")
    print()

    return (total_preds, preds)


# def run_analysis(intersect_bed_filename, predictions_gff, gc_prob, output_dir, output_prefix=None):
#
#     if output_prefix is None:
#         # Generate output prefix from input filenames
#         bed_base = os.path.splitext(os.path.basename(intersect_bed_filename))[0]
#         gff_base = os.path.splitext(os.path.basename(predictions_gff))[0]
#         output_prefix = f"{bed_base}_{gff_base}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
#
#     print(f"Starting analysis with output prefix: {output_prefix}")
#
#     # Read predictions
#     total_preds, preds = read_preds(predictions_gff, output_dir, output_prefix)
#
#     # Run evaluation
#     stats, answer_counts, codon_counts, context_counts = evaluate(intersect_bed_filename, preds, gc_prob, output_dir, output_prefix)
#
#     print(f"\nAnalysis complete! Generated CSV files:")
#     print(os.path.join(output_dir, f"  - {output_prefix}_run_summary.csv"))
#     print(os.path.join(output_dir, f"  - {output_prefix}_detailed_results.csv"))
#     print(os.path.join(output_dir, f"  - {output_prefix}_codon_summary.csv"))
#     print(os.path.join(output_dir, f"  - {output_prefix}_prediction_lengths.csv"))
#     print(f"  - Individual codon analysis files: {output_prefix}_codons_*.csv")
#
#     return stats, answer_counts, codon_counts, context_counts