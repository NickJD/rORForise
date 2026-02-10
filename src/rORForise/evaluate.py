import collections
import csv
import os
import gzip
import math
from datetime import datetime


try:
    import rORForise.check_pred as cp
except Exception:
    try:
        # Fallback to local relative import when executed inside package context
        from . import check_pred as cp
    except Exception:
        # Final fallback: load from file path (works when script is executed from arbitrary cwd)
        import importlib.util
        from pathlib import Path
        candidate = Path(__file__).resolve().parent / 'check_pred.py'
        if candidate.exists():
            spec = importlib.util.spec_from_file_location('rORForise_check_pred', str(candidate))
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            cp = mod
        else:
            raise


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

    if not codon_counts:
        print("No codon data available for analysis")
        return

    summary_file = os.path.join(output_dir, f"{output_prefix}_codon_summary.csv")
    with open(summary_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['answer_type', 'codon', 'observed', 'expected', 'obs_exp_ratio'])

        for answer_code, answer_name in codon_answer_types.items():
            if answer_code in codon_counts:
                chosen = codon_counts[answer_code]
                total = sum(chosen.values())

                for codon in codons:
                    observed = chosen.get(codon, 0)
                    if answer_name in ['incorrect stop', 'incorrect start', 'middle incorrect start',
                                       'middle incorrect stop']:
                        expected_val = expect(codon, total, gc_prob)
                        ratio = observed / expected_val if expected_val > 0 else 0
                    else:
                        expected_val = 'N/A'
                        ratio = 'N/A'

                    writer.writerow([answer_name, codon, observed, expected_val, ratio])

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
                    observed = chosen.get(codon, 0)
                    frequency = observed / total if total > 0 else 0

                    if answer_name in ['incorrect stop', 'incorrect start', 'middle incorrect start',
                                       'middle incorrect stop']:
                        expected_val = expect(codon, total, gc_prob)
                        ratio = observed / expected_val if expected_val > 0 else 0
                        enrichment = 'enriched' if ratio > 1.5 else 'depleted' if ratio < 0.5 else 'normal'
                        writer.writerow([codon, observed, f"{expected_val:.1f}", f"{ratio:.2f}", enrichment])
                    else:
                        writer.writerow([codon, observed, f"{frequency:.3f}"])


def write_summary_csv(stats, output_dir, output_prefix):

    summary_file = os.path.join(output_dir, f"{output_prefix}_run_summary.csv")

    with open(summary_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['metric', 'value'])

        writer.writerow(['timestamp', datetime.now().isoformat()])
        writer.writerow(['total_number_of_reads', stats['total_number_of_reads']])
        writer.writerow(['total_cds_aligned_reads', stats['total_cds_aligned_reads']])
        writer.writerow(['reads_without_predictions', stats['reads_without_predictions']])
        writer.writerow(['reads_good_overlap', stats['reads_good_overlap']])
        writer.writerow(['reads_good_overlap_with_preds', stats['reads_good_overlap_with_preds']])
        writer.writerow(['total_predictions', stats['total_predictions']])
        writer.writerow(['on_target_predictions', stats['on_target_predictions']])
        writer.writerow(['shared_reads', stats['shared_reads']])
        writer.writerow(['shared_predictions', stats['shared_predictions']])

        if stats['reads_good_overlap'] > 0:
            writer.writerow(['prediction_coverage',
                             f"{stats['reads_good_overlap_with_preds'] / stats['reads_good_overlap']:.3f}"])
        writer.writerow(['avg_preds_per_read',
                         f"{stats['total_predictions'] / stats['total_number_of_reads']:.2f}"])


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
        for read_name in sorted(track_preds.keys()):
            for pred_start, pred_end, answer_details in track_preds[read_name]:
                for answer_code, codon in answer_details:
                    answer_name = cp.inverse_answers.get(answer_code, str(answer_code))
                    writer.writerow([read_name, pred_start, pred_end, answer_code, answer_name, codon if codon is not None else ''])


def write_detailed_results_csv(answer_counts, output_dir, output_prefix, stats, context_counts):
    # Implementation ported from previous corrupted version, but fixed
    results_file = os.path.join(output_dir, f"{output_prefix}_detailed_results.csv")

    with open(results_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['answer_type', 'count', 'percentage', 'denominator', 'context'])

        start_only_answer_types = [
            cp.answers['correct start'],
            cp.answers['incorrect start'],
            cp.answers['alternative start']
        ]

        stop_only_answer_types = [
            cp.answers['correct stop'],
            cp.answers['incorrect stop'],
            cp.answers['alternative stop']
        ]

        middle_start_answer_types = [
            cp.answers['middle or alternative start'],
            cp.answers['middle incorrect start']
        ]

        middle_stop_answer_types = [
            cp.answers['middle or alternative stop'],
            cp.answers['middle incorrect stop']
        ]

        frame_answer_types = [
            cp.answers['correct frame'],
            cp.answers['incorrect frame']
        ]

        direction_answer_types = [
            cp.answers['correct direction'],
            cp.answers['incorrect direction']
        ]

        middle_only_answer_types = [
            cp.answers['middle']
        ]

        boundary_answer_types = [
            cp.answers['prediction ends before cds starts'],
            cp.answers['prediction starts after cds ends']
        ]

        overlap_answer_types = [
            cp.answers['not enough read-CDS overlap']
        ]

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

        for answer_code in sorted(all_answer_codes):
            answer_name = cp.inverse_answers[answer_code]
            count = answer_counts.get(answer_code, 0)

            if answer_code in start_only_answer_types:
                # start-related answers may also appear when a prediction is middle-only
                denominator = max(context_counts.get('overlaps_start', 0), context_counts.get('middle_only', 0))
                context = 'predictions overlapping CDS start or in CDS middle'

            elif answer_code in stop_only_answer_types:
                # stop-related answers may also appear when a prediction is middle-only
                denominator = max(context_counts.get('overlaps_stop', 0), context_counts.get('middle_only', 0))
                context = 'predictions overlapping CDS stop or in CDS middle'

            elif answer_code in middle_start_answer_types:
                denominator = context_counts.get('middle_only', 0)
                context = 'predictions in CDS middle only'

            elif answer_code in middle_stop_answer_types:
                denominator = context_counts.get('middle_only', 0)
                context = 'predictions in CDS middle only'

            elif answer_code in middle_only_answer_types:
                denominator = context_counts.get('middle_only', 0)
                context = 'predictions in CDS middle only'

            elif answer_code in frame_answer_types:
                denominator = stats.get('on_target_predictions', 0)
                context = 'all on-target predictions'

            elif answer_code in direction_answer_types:
                denominator = stats.get('on_target_predictions', 0)
                context = 'all on-target predictions'

            elif answer_code in boundary_answer_types:
                denominator = stats.get('on_target_predictions', 0)
                context = 'all on-target predictions'

            elif answer_code in overlap_answer_types:
                denominator = stats.get('total_cds_aligned_reads', 0)
                context = 'all CDS-aligned reads'

            else:
                denominator = stats.get('on_target_predictions', 0)
                context = 'all on-target predictions (uncategorised)'

            percentage = (count / denominator * 100) if denominator > 0 else 0
            writer.writerow([answer_name, count, f"{percentage:.2f}", denominator, context])


def evaluate(intersect_bed_filename, preds, gc_prob, output_dir, output_prefix="orf_evaluation", overlap_threshold=60, verbose=False, output_format='csv'):
    number_of_CDS_mappings_with_predictions = 0
    number_of_on_target_preds = 0

    track_preds = collections.defaultdict(list)

    reads_without_predictions = set()
    seen_read_names = set()
    good_overlap_read_names = set()

    answer_counts = collections.defaultdict(int)
    codon_counts = {}

    context_counts = {
        'overlaps_start': 0,
        'overlaps_stop': 0,
        'middle_only': 0
    }

    # Read intersect BED (supports gz)
    opener = (gzip.open if str(intersect_bed_filename).endswith('.gz') else open)
    with opener(intersect_bed_filename, 'rt', encoding='utf-8') as f:
        csvr = csv.reader(f, delimiter='\t')
        try:
            header = next(csvr)
        except StopIteration:
            header = []
        read_num = 0
        unique_read_names = set()
        for bed_row in csvr:
            read_num += 1
            if len(bed_row) < 11:
                continue
            # count each read once regardless of multiple BED rows
            read_name_all = bed_row[1]
            unique_read_names.add(read_name_all)
            if bed_row[6] == 'CDS':
                read_start = int(bed_row[2])
                read_end = int(bed_row[3])
                read_name = bed_row[1]
                read_dir = bed_row[4]
                seen_read_names.add(read_name)
                cds_start = int(bed_row[7])
                cds_end = int(bed_row[8])

                if min(read_end, cds_end) - max(read_start, cds_start) + 1 < overlap_threshold:
                    answer_counts[cp.answers["not enough read-CDS overlap"]] += 1
                    if read_name not in preds:
                        reads_without_predictions.add(read_name)
                else:
                    good_overlap_read_names.add(read_name)

                    if read_name in preds:
                        cds_dir = bed_row[9]
                        read_seq = bed_row[10]
                        number_of_CDS_mappings_with_predictions += 1

                        for prediction_name, (pred_start, pred_end, pred_dir) in preds[read_name].items():
                            number_of_on_target_preds += 1

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

                            if read_captures_cds_start:
                                context_counts['overlaps_start'] += 1
                            if read_captures_cds_end:
                                context_counts['overlaps_stop'] += 1
                            if not read_captures_cds_start and not read_captures_cds_end:
                                context_counts['middle_only'] += 1

                            answer_details = cp.check_pred(cds_start, cds_end, cds_dir, read_start, read_end, read_dir,
                                                           pred_start, pred_end, pred_dir, read_seq)
                            track_preds[read_name].append((pred_start, pred_end, answer_details))

                            for (answer, codon) in answer_details:
                                answer_counts[answer] += 1
                                if codon is not None:
                                    if answer not in codon_counts:
                                        codon_counts[answer] = collections.defaultdict(int)
                                    codon_counts[answer][codon] += 1

                    else:
                        reads_without_predictions.add(read_name)

    stats = {
        'total_number_of_reads': len(unique_read_names),
        'total_cds_aligned_reads': len(seen_read_names),
        'reads_without_predictions': len(reads_without_predictions),
        'reads_good_overlap': len(good_overlap_read_names),
        'reads_good_overlap_with_preds': number_of_CDS_mappings_with_predictions,
        'total_predictions': sum(len(pred_list) for pred_list in preds.values()),
        'on_target_predictions': number_of_on_target_preds
    }

    shared_reads = set(preds.keys()).intersection(seen_read_names)
    shared_predictions = sum(len(preds[read]) for read in shared_reads)
    stats['shared_reads'] = len(shared_reads)
    stats['shared_predictions'] = shared_predictions

    os.makedirs(output_dir, exist_ok=True)

    write_summary_csv(stats, output_dir, output_prefix)
    write_detailed_results_csv(answer_counts, output_dir, output_prefix, stats, context_counts)
    write_codon_analysis_csv(codon_counts, gc_prob, output_dir, output_prefix)
    write_track_predictions(track_preds, output_dir, output_prefix)

    # Prediction length distribution
    pred_length_distn = collections.defaultdict(int)
    for pred_list in preds.values():
        for _, (ps, pe, _) in pred_list.items() if isinstance(pred_list, dict) else []:
            pred_length_distn[pe - ps + 1] += 1

    write_prediction_length_csv(pred_length_distn, output_dir, output_prefix)

    # Compute simple context-specific accuracies (start/stop/middle)
    metrics = {}
    start_den = context_counts.get('overlaps_start', 0)
    stop_den = context_counts.get('overlaps_stop', 0)
    middle_den = context_counts.get('middle_only', 0)

    correct_start = answer_counts.get(cp.answers.get('correct start', -1), 0)
    correct_stop = answer_counts.get(cp.answers.get('correct stop', -1), 0)
    correct_middle = answer_counts.get(cp.answers.get('middle', -1), 0)

    def pct(num, den):
        if den <= 0:
            return 'N/A'
        return f"{(num / den * 100):.1f}%"

    metrics['start'] = (correct_start, start_den, pct(correct_start, start_den))
    metrics['stop'] = (correct_stop, stop_den, pct(correct_stop, stop_den))
    metrics['middle'] = (correct_middle, middle_den, pct(correct_middle, middle_den))

    # write metrics CSV
    metrics_file = os.path.join(output_dir, f"{output_prefix}_context_accuracy.csv")
    with open(metrics_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['context', 'correct_count', 'denominator', 'accuracy_pct'])
        for ctx in ['start', 'stop', 'middle']:
            c, d, p = metrics[ctx]
            writer.writerow([ctx, c, d, p])

    # print friendly normalised metrics
    if verbose:
        print('\nContext accuracies:')
        for ctx in ['start', 'stop', 'middle']:
            c, d, p = metrics[ctx]
            print(f"  {ctx.title():6}: {p} ({c}/{d})")
    else:
        print('\nContext accuracies written to: ' + metrics_file)

    print(f"\nAnalysis complete! Generated CSV files:")
    print("Output Directory: " + str(os.path.join(output_dir)))
    print(f"  - {output_prefix}_run_summary.csv")
    print(f"  - {output_prefix}_detailed_results.csv")
    print(f"  - {output_prefix}_codon_summary.csv")
    print(f"  - {output_prefix}_prediction_lengths.csv")
    print(f"Individual codon analysis files: {output_prefix}_codons_*.csv")


def read_preds(predictions_gff, output_dir, output_prefix="orf_evaluation", verbose=False):
    preds = {}
    total_preds = 0
    pred_length_distn = collections.defaultdict(int)

    opener = (gzip.open if str(predictions_gff).endswith('.gz') else open)
    with opener(predictions_gff, 'rt', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            row = line.split('\t')
            if len(row) < 9:
                continue
            read_name = row[0].lstrip('@')
            pred_start = int(row[3])
            pred_end = int(row[4])
            pred_dir = row[6]
            attr = row[8]
            prediction_name = attr.replace('ID=', '').split(';')[0].replace('@', '')

            if read_name not in preds:
                preds[read_name] = {}
            preds[read_name][prediction_name] = (pred_start, pred_end, pred_dir)
            total_preds += 1
            pred_length_distn[pred_end - pred_start + 1] += 1

    write_prediction_length_csv(pred_length_distn, output_dir, output_prefix)

    if verbose:
        print(f"Loaded {total_preds} predictions from {len(preds)} reads")
        print(f"Prediction length distribution written to {output_prefix}_prediction_lengths.csv")
        print()

    return (total_preds, preds)


if __name__ == '__main__':
    print('This module provides evaluate(), read_preds(), and CSV writers for analysis.')
