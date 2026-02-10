"""
Generate synthetic test data for rORForise evaluation.

Creates:
1. An intersect BED file with read-to-CDS mappings
2. A predictions GFF file with ORF predictions
3. A validation CSV file documenting expected outcomes

This allows testing of the evaluation pipeline with controlled scenarios.
"""

import random
import gzip
import argparse
import csv
from pathlib import Path

# robust import of check_pred so generator can compute canonical expected answers
try:
    from rORForise import check_pred as cp
except Exception:
    try:
        import check_pred as cp
    except Exception:
        import importlib.util, os
        repo_root = Path(__file__).resolve().parents[1]
        candidate = repo_root / 'src' / 'rORForise' / 'check_pred.py'
        if candidate.exists():
            spec = importlib.util.spec_from_file_location('check_pred', str(candidate))
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            cp = mod
        else:
            cp = None


def generate_sequence(length, gc_content=0.32):
    bases = []
    for _ in range(length):
        if random.random() < gc_content:
            bases.append(random.choice(["G", "C"]))
        else:
            bases.append(random.choice(["A", "T"]))
    return ''.join(bases)


def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))


def insert_codon_at_position(seq, position, codon):
    # Replace three bases starting at position with codon (0-based)
    if position < 0 or position + 3 > len(seq):
        return seq
    return seq[:position] + codon + seq[position + 3:]


def generate_cds_features(num_cds=50, min_length=300, max_length=1500, genome_size=580076):
    cds_features = []
    for _ in range(num_cds):
        length = random.randint(min_length, max_length)
        length = (length // 3) * 3
        start = random.randint(1000, genome_size - length - 1000)
        end = start + length - 1
        # allow CDS on either strand
        strand = random.choice(['+', '-'])
        cds_features.append((start, end, strand))
    cds_features.sort(key=lambda x: x[0])
    return cds_features


def generate_read_mapping(cds_start, cds_end, cds_strand, read_length=150, scenario='correct_start'):
    # Simple deterministic mappings that approximate scenarios
    if scenario in ('correct_start', 'left_overhang'):
        overhang = random.randint(10, 50)
        read_start = cds_start - overhang
        read_end = read_start + read_length - 1
        read_strand = cds_strand
    elif scenario in ('correct_stop', 'right_overhang'):
        overhang = random.randint(10, 50)
        read_end = cds_end + overhang
        read_start = read_end - read_length + 1
        read_strand = cds_strand
    elif scenario == 'correct_both':
        read_length = (cds_end - cds_start + 1) + random.randint(20, 100)
        overhang = random.randint(10, 50)
        read_start = cds_start - overhang
        read_end = read_start + read_length - 1
        read_strand = cds_strand
    elif scenario == 'middle':
        # Put read entirely inside CDS
        if cds_end - cds_start + 1 <= read_length + 20:
            read_start = cds_start + 5
        else:
            read_start = random.randint(cds_start + 10, cds_end - read_length - 10)
        read_end = read_start + read_length - 1
        read_strand = cds_strand
    elif scenario == 'poor_overlap':
        overlap = random.randint(20, 55)
        read_start = cds_start - (read_length - overlap)
        read_end = read_start + read_length - 1
        read_strand = cds_strand
    else:
        read_start = random.randint(cds_start + 10, max(cds_start + 10, cds_end - read_length - 10))
        read_end = read_start + read_length - 1
        read_strand = cds_strand

    # occasionally flip read strand to simulate reads mapping to opposite strand
    if random.random() < 0.3:
        read_strand = '-' if read_strand == '+' else '+'

    return read_start, read_end, read_strand


def generate_prediction(read_start, read_end, read_strand, cds_start, cds_end, cds_strand,
                        read_seq, scenario='correct_start', frame_offset=0):
    read_length = read_end - read_start + 1
    category = (cds_strand, read_strand)

    expected = {
        'expected_answers': [],
        'pred_scenario': scenario,
        'should_overlap_start': False,
        'should_overlap_stop': False,
        'is_middle_only': False,
        'pred_start_codon': None,
        'pred_stop_codon': None
    }

    # Map (cds_strand, read_strand) -> preferred pred_strand so that check_pred
    # handles the prediction as being in the "correct direction" category.
    pred_strand_map = {
        ('+', '+'): '+',
        ('+', '-'): '-',
        ('-', '+'): '-',
        ('-', '-'): '+'
    }
    default_pred_strand = pred_strand_map.get(category, read_strand)

    def make_pred(start_on_read, length=120, strand=read_strand):
        end_on_read = min(read_length, start_on_read + length - 1)
        return start_on_read, end_on_read, strand

    def map_cds_coord_to_pred(cds_coord, which='start', length=120, pred_strand=default_pred_strand):
        """Map a genome coordinate (cds_coord) to a prediction (start,end) on the read.
        This computes positions relative to the read sequence passed into this
        function (the read_seq should already be oriented according to read_strand).
        """
        if read_strand == '+':
            pos = cds_coord - read_start + 1
        else:
            # read is reverse-oriented relative to genome; first base of read_seq
            # corresponds to genome coordinate read_end
            pos = read_end - cds_coord + 1

        if which == 'start':
            pstart = pos
            pend = pos + length - 1
        else:
            pend = pos
            pstart = pos - length + 1

        # clamp into read coordinates
        pstart = max(1, pstart)
        pend = max(1, min(read_length, pend))
        return pstart, pend

    # Simple heuristics for predictions depending on scenario
    if scenario == 'correct_start':
        if category in (('+', '+'), ('+', '-')) and read_start <= cds_start:
            pred_start_on_read = max(1, cds_start - read_start + 1 + frame_offset)
            pred_start_on_read, pred_end_on_read = map_cds_coord_to_pred(cds_start, 'start', 120, default_pred_strand)
            pred = (pred_start_on_read, pred_end_on_read, default_pred_strand)
            # insert ATG
            pos = pred[0] - 1
            if 0 <= pos <= len(read_seq) - 3:
                read_seq = insert_codon_at_position(read_seq, pos, 'ATG')
                expected['pred_start_codon'] = 'ATG'
            expected['should_overlap_start'] = True
            expected['expected_answers'] = ['correct direction']
            if frame_offset == 0:
                expected['expected_answers'].extend(['correct start', 'correct frame'])
            elif frame_offset % 3 == 0:
                expected['expected_answers'].extend(['alternative start', 'correct frame'])
            else:
                expected['expected_answers'].extend(['incorrect start', 'incorrect frame'])
            # return pred_start, pred_end, pred_strand, read_seq, expected
            return (pred[0], pred[1], pred[2], read_seq, expected)

    if scenario == 'correct_stop':
        if category in (('+', '+'), ('+', '-')) and read_end >= cds_end:
            pred_start_on_read, pred_end_on_read = map_cds_coord_to_pred(cds_end, 'end', 150, default_pred_strand)
            pos = pred_end_on_read - 3
            if 0 <= pos <= len(read_seq) - 3:
                read_seq = insert_codon_at_position(read_seq, pos, 'TAA')
                expected['pred_stop_codon'] = 'TAA'
            expected['should_overlap_stop'] = True
            expected['expected_answers'] = ['correct direction']
            if frame_offset == 0:
                expected['expected_answers'].extend(['correct stop', 'correct frame'])
            elif frame_offset % 3 == 0:
                expected['expected_answers'].extend(['alternative stop', 'correct frame'])
            else:
                expected['expected_answers'].extend(['incorrect stop', 'incorrect frame'])
            return (pred_start_on_read, pred_end_on_read, default_pred_strand, read_seq, expected)

    if scenario == 'alternative_start':
        if category in (('+', '+'), ('+', '-')) and read_start <= cds_start:
            pred_start_on_read, pred_end_on_read = map_cds_coord_to_pred(cds_start, 'start', 120, default_pred_strand)
            pos = pred_start_on_read - 1
            alt_codon = random.choice(['GTG', 'TTG'])
            if 0 <= pos <= len(read_seq) - 3:
                read_seq = insert_codon_at_position(read_seq, pos, alt_codon)
                expected['pred_start_codon'] = alt_codon
            expected['should_overlap_start'] = True
            expected['expected_answers'] = ['correct direction', 'alternative start', 'correct frame']
            return (pred_start_on_read, pred_end_on_read, default_pred_strand, read_seq, expected)

    if scenario == 'incorrect_start':
        if category in (('+', '+'), ('+', '-')) and read_start <= cds_start:
            # place an out-of-frame start by offsetting a base or two from true start
            pstart, pend = map_cds_coord_to_pred(cds_start, 'start', 120, default_pred_strand)
            offset = random.choice([1, 2])
            pred_start_on_read = pstart + offset
            expected['should_overlap_start'] = True
            expected['expected_answers'] = ['correct direction', 'incorrect start', 'incorrect frame']
            return (pred_start_on_read, min(read_length, pred_start_on_read + 150 - 1), default_pred_strand, read_seq, expected)

    if scenario == 'middle':
        if not (read_start <= cds_start or read_end >= cds_end):
            pred_start_on_read = random.randint(1, max(1, read_length - 100))
            pred_end_on_read = min(read_length, pred_start_on_read + 100)
            expected['is_middle_only'] = True
            expected['expected_answers'] = ['correct direction', 'middle']
            return (pred_start_on_read, pred_end_on_read, default_pred_strand, read_seq, expected)

    if scenario == 'wrong_direction':
        pred_start_on_read = random.randint(1, max(1, read_length - 100))
        pred_end_on_read = min(read_length, pred_start_on_read + 100)
        pred_strand = '-' if read_strand == '+' else '+'
        expected['expected_answers'] = ['incorrect direction']
        return (pred_start_on_read, pred_end_on_read, pred_strand, read_seq, expected)

    # fallback: random prediction on read
    pred_start_on_read = random.randint(1, max(1, read_length - 60))
    pred_end_on_read = min(read_length, pred_start_on_read + random.randint(60, 150))
    expected['expected_answers'] = ['correct direction']
    return (pred_start_on_read, pred_end_on_read, default_pred_strand, read_seq, expected)


def write_intersect_bed(filename, entries, gzip_out=False):
    header = ['chrom', 'read_name', 'read_start', 'read_end', 'read_strand',
              'score', 'feature_type', 'cds_start', 'cds_end', 'cds_strand', 'read_sequence']
    if gzip_out or str(filename).endswith('.gz'):
        with gzip.open(filename, 'wt', encoding='utf-8') as f:
            f.write('\t'.join(header) + '\n')
            for entry in entries:
                f.write('\t'.join(map(str, entry)) + '\n')
    else:
        with open(filename, 'w', encoding='utf-8') as f:
            f.write('\t'.join(header) + '\n')
            for entry in entries:
                f.write('\t'.join(map(str, entry)) + '\n')


def write_predictions_gff(filename, predictions, gzip_out=False):
    if gzip_out or str(filename).endswith('.gz'):
        with gzip.open(filename, 'wt', encoding='utf-8') as f:
            f.write('##gff-version 3\n')
            for pred in predictions:
                f.write('\t'.join(map(str, pred)) + '\n')
    else:
        with open(filename, 'w', encoding='utf-8') as f:
            f.write('##gff-version 3\n')
            for pred in predictions:
                f.write('\t'.join(map(str, pred)) + '\n')


def write_validation_csv(filename, validation_data, gzip_out=False):
    header = ['read_name',  'read_start', 'read_end', 'read_strand',
              'pred_id', 'pred_start', 'pred_end', 'pred_strand',
              'cds_id', 'cds_start', 'cds_end', 'cds_strand',
              'pred_scenario', 'should_overlap_start', 'should_overlap_stop',
              'is_middle_only', 'expected_answers', 'pred_start_codon', 'pred_stop_codon']

    if gzip_out or str(filename).endswith('.gz'):
        with gzip.open(filename, 'wt', encoding='utf-8', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(header)
            for entry in validation_data:
                writer.writerow(entry)
    else:
        with open(filename, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow(header)
            for entry in validation_data:
                writer.writerow(entry)


def main():
    parser = argparse.ArgumentParser(description="Generate synthetic test data for rORForise")
    parser.add_argument('-n', '--num_reads', type=int, default=200, help='Number of reads to generate (default: 200)')
    parser.add_argument('-c', '--num_cds', type=int, default=50, help='Number of CDS features to generate (default: 50)')
    parser.add_argument('-o', '--output_dir', type=str, default='test_output', help='Directory to write output files')
    parser.add_argument('--output_prefix', type=str, default='test', help='Prefix for output files')
    parser.add_argument('-gc', '--gc_content', type=float, default=0.32, help='GC content for sequence generation (default: 0.32)')
    parser.add_argument('--seed', type=int, default=42, help='Random seed for reproducibility (default: 42)')
    parser.add_argument('--gzip', action='store_true', help='Write output files compressed with gzip (.gz)')
    parser.add_argument('--min_pred_len', type=int, default=30, help='Minimum length of predictions for GFF/validation (default: 30)')
    parser.add_argument('--strand_mode', choices=['mixed','plus','minus'], default='mixed', help='Strand generation mode: mixed (default), plus, or minus')
    # Fragment and read sizing for paired-end simulation
    parser.add_argument('--fragment_length_mean', type=float, default=500.0, help='Mean fragment (insert) length in bases (default: 500).')
    parser.add_argument('--fragment_length_sd', type=float, default=30.0, help='SD for fragment length sampling (default: 30). If 0, fragments are fixed at the mean.')
    parser.add_argument('--read_length', type=int, default=150, help='Length of one read in a pair (default: 150).')
    # Paired-end options and orientation
    parser.add_argument('--paired', dest='paired', action='store_true', help='Simulate paired-end reads (default).')
    parser.add_argument('--no-paired', dest='paired', action='store_false', help='Do not simulate paired reads; treat reads as single contiguous fragments.')
    parser.add_argument('--pair_orientation', choices=['FR','RF','FF'], default='FR', help='Pair orientation: FR (typical Illumina inward-facing), RF, or FF (both same strand).')
    parser.set_defaults(paired=True)

    args = parser.parse_args()
    random.seed(args.seed)

    # strand mode controls whether CDS/read strands are generated on both strands
    global STRAND_MODE
    STRAND_MODE = args.strand_mode

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    cds_features = generate_cds_features(num_cds=args.num_cds)

    bed_entries = []
    gff_predictions = []
    validation_data = []

    scenarios = [
        'correct_start', 'correct_stop', 'alternative_start',
        'incorrect_start', 'middle', 'wrong_direction'
    ]
    scenario_weights = [0.25, 0.25, 0.15, 0.15, 0.15, 0.05]

    scenario_counts = {s: 0 for s in scenarios}

    for read_idx in range(args.num_reads):
        cds_start, cds_end, cds_strand = random.choice(cds_features)
        pred_scenario = random.choices(scenarios, weights=scenario_weights)[0]
        scenario_counts[pred_scenario] += 1

        if pred_scenario in ['correct_start', 'alternative_start', 'incorrect_start']:
            mapping_scenario = 'left_overhang'
        elif pred_scenario == 'correct_stop':
            mapping_scenario = 'right_overhang'
        elif pred_scenario == 'middle':
            mapping_scenario = 'middle'
        else:
            mapping_scenario = 'middle'

        if args.paired:
            # sample fragment (insert) length per fragment from normal distribution
            single_len = args.read_length
            if args.fragment_length_sd <= 0:
                frag_len = int(round(args.fragment_length_mean))
            else:
                frag_len = int(round(random.normalvariate(args.fragment_length_mean, args.fragment_length_sd)))
            # enforce minimum fragment length to contain two reads
            if frag_len < single_len * 2 + 1:
                frag_len = single_len * 2 + 1

            # map a fragment (using existing helper) then derive mate positions
            frag_start, frag_end, frag_strand = generate_read_mapping(cds_start, cds_end, cds_strand, frag_len, mapping_scenario)
            fragment_seq = generate_sequence(frag_end - frag_start + 1, args.gc_content)

            # read 1 at fragment start, read 2 at fragment end
            r1_start = frag_start
            r1_end = frag_start + single_len - 1
            r2_end = frag_end
            r2_start = frag_end - single_len + 1

            # strands based on requested pair_orientation
            if args.pair_orientation == 'FR':
                r1_strand = frag_strand
                r2_strand = '-' if frag_strand == '+' else '+'
            elif args.pair_orientation == 'RF':
                r1_strand = '-' if frag_strand == '+' else '+'
                r2_strand = frag_strand
            else:  # 'FF'
                r1_strand = frag_strand
                r2_strand = frag_strand

            # extract sequences oriented to their reported strand
            r1_seq = fragment_seq[0:single_len]
            r2_seq = fragment_seq[-single_len:]
            if r1_strand == '-':
                r1_seq = reverse_complement(r1_seq)
            if r2_strand == '-':
                r2_seq = reverse_complement(r2_seq)

            read_name_base = f"read_{read_idx:06d}"
            read1_name = read_name_base + "/1"
            read2_name = read_name_base + "/2"

            # create BED entries for both mates
            bed_entries.append(['Chromosome-0', read1_name, r1_start, r1_end, r1_strand, '.', 'CDS', cds_start, cds_end, cds_strand, r1_seq])
            bed_entries.append(['Chromosome-0', read2_name, r2_start, r2_end, r2_strand, '.', 'CDS', cds_start, cds_end, cds_strand, r2_seq])

            # generate predictions for each mate independently
            for mate_name, rs, re, rstrand, rseq, mate_idx in [
                (read1_name, r1_start, r1_end, r1_strand, r1_seq, 1),
                (read2_name, r2_start, r2_end, r2_strand, r2_seq, 2)
            ]:
                frame_offset = 0 if pred_scenario != 'incorrect_start' else random.choice([1, 2])
                pred_start, pred_end, pred_strand, modified_seq, expected = generate_prediction(
                    rs, re, rstrand, cds_start, cds_end, cds_strand, rseq, pred_scenario, frame_offset
                )
                # update seq if modified
                rseq = modified_seq
                if pred_end - pred_start + 1 >= args.min_pred_len:
                    pred_id = f"{mate_name}_pred_{0}"
                    gff_predictions.append([f"@{mate_name}", 'TestGenerator', 'CDS', pred_start, pred_end, '.', pred_strand, '0', f"ID={pred_id}"])
                    if cp is not None:
                        answers = cp.check_pred(cds_start, cds_end, cds_strand, rs, re, rstrand, pred_start, pred_end, pred_strand, rseq)
                        answer_names = '; '.join(sorted({cp.inverse_answers[a] for a, _ in answers}))
                        pred_start_codon = ''
                        pred_stop_codon = ''
                        for a, cod in answers:
                            name = cp.inverse_answers.get(a)
                            if name and 'start' in name and cod:
                                pred_start_codon = cod
                            if name and 'stop' in name and cod:
                                pred_stop_codon = cod
                    else:
                        answer_names = '; '.join(expected['expected_answers'])
                        pred_start_codon = expected.get('pred_start_codon', '')
                        pred_stop_codon = expected.get('pred_stop_codon', '')
                    validation_data.append([mate_name, rs, re, rstrand, pred_id, pred_start, pred_end, pred_strand, f"cds_{cds_start}_{cds_end}", cds_start, cds_end, cds_strand, expected['pred_scenario'], expected['should_overlap_start'], expected['should_overlap_stop'], expected['is_middle_only'], answer_names, pred_start_codon, pred_stop_codon])

        else:
            # single-read behavior (backwards-compatible): args.read_length is read length
            read_length = args.read_length
            read_start, read_end, read_strand = generate_read_mapping(cds_start, cds_end, cds_strand, read_length, mapping_scenario)
            read_name = f"read_{read_idx:06d}"
            read_seq_length = read_end - read_start + 1
            read_seq = generate_sequence(read_seq_length, args.gc_content)

            frame_offset = 0 if pred_scenario != 'incorrect_start' else random.choice([1, 2])

            if read_strand == '-':
                read_seq = reverse_complement(read_seq)

            pred_start, pred_end, pred_strand, modified_seq, expected = generate_prediction(read_start, read_end, read_strand, cds_start, cds_end, cds_strand, read_seq, pred_scenario, frame_offset)
            read_seq = modified_seq

            bed_entry = ['Chromosome-0', read_name, read_start, read_end, read_strand, '.', 'CDS', cds_start, cds_end, cds_strand, read_seq]
            bed_entries.append(bed_entry)

            if pred_end - pred_start + 1 >= args.min_pred_len:
                pred_id = f"{read_name}_pred_0"
                gff_entry = [f"@{read_name}", 'TestGenerator', 'CDS', pred_start, pred_end, '.', pred_strand, '0', f"ID={pred_id}"]
                gff_predictions.append(gff_entry)

                if cp is not None:
                    answers = cp.check_pred(cds_start, cds_end, cds_strand, read_start, read_end, read_strand, pred_start, pred_end, pred_strand, read_seq)
                    answer_names = '; '.join(sorted({cp.inverse_answers[a] for a, _ in answers}))
                    pred_start_codon = ''
                    pred_stop_codon = ''
                    for a, cod in answers:
                        name = cp.inverse_answers.get(a)
                        if name and 'start' in name and cod:
                            pred_start_codon = cod
                        if name and 'stop' in name and cod:
                            pred_stop_codon = cod
                else:
                    answer_names = '; '.join(expected['expected_answers'])
                    pred_start_codon = expected.get('pred_start_codon', '')
                    pred_stop_codon = expected.get('pred_stop_codon', '')

                validation_entry = [read_name, read_start, read_end, read_strand, pred_id, pred_start, pred_end, pred_strand, f"cds_{cds_start}_{cds_end}", cds_start, cds_end, cds_strand, expected['pred_scenario'], expected['should_overlap_start'], expected['should_overlap_stop'], expected['is_middle_only'], answer_names, pred_start_codon, pred_stop_codon]
                validation_data.append(validation_entry)

    bed_file = outdir / f"{args.output_prefix}_intersect.bed"
    gff_file = outdir / f"{args.output_prefix}_predictions.gff"
    validation_file = outdir / f"{args.output_prefix}_validation.csv"

    write_intersect_bed(str(bed_file) + ('.gz' if args.gzip else ''), bed_entries, gzip_out=args.gzip)
    write_predictions_gff(str(gff_file) + ('.gz' if args.gzip else ''), gff_predictions, gzip_out=args.gzip)
    write_validation_csv(str(validation_file) + ('.gz' if args.gzip else ''), validation_data, gzip_out=args.gzip)

    print(f"Wrote {len(bed_entries)} BED entries to {bed_file}")
    print(f"Wrote {len(gff_predictions)} predictions to {gff_file}")
    print(f"Wrote {len(validation_data)} validation rows to {validation_file}")
    print("Scenario distribution:")
    for s, c in scenario_counts.items():
        print(f"  {s}: {c}")


if __name__ == '__main__':
    main()
