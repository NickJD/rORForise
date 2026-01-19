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


def generate_sequence(length, gc_content=0.32):
    # Generate a random DNA sequence with specified GC content.
    bases = []
    for _ in range(length):
        if random.random() < gc_content:
            bases.append(random.choice(['G', 'C']))
        else:
            bases.append(random.choice(['A', 'T']))
    return ''.join(bases)


def reverse_complement(seq):
    # Return the reverse complement of a DNA sequence.
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))


def insert_codon_at_position(seq, position, codon):
    # Insert a specific codon at a position in the sequence (0-based).
    if position < 0 or position + 3 > len(seq):
        return seq
    return seq[:position] + codon + seq[position + 3:]


def generate_cds_features(genome_name, num_cds=50, min_length=300, max_length=1500, genome_size=580076):
    # Generate synthetic CDS features distributed across a genome.
    # Returns list of tuples: (cds_start, cds_end, cds_strand)

    cds_features = []

    for i in range(num_cds):
        length = random.randint(min_length, max_length)
        # Ensure length is divisible by 3
        length = (length // 3) * 3

        # Random start position, ensuring we don't go past genome end
        start = random.randint(1000, genome_size - length - 1000)
        end = start + length - 1
        strand = random.choice(['+', '-'])

        cds_features.append((start, end, strand))

    # Sort by start position
    cds_features.sort(key=lambda x: x[0])

    return cds_features


def generate_read_mapping(cds_start, cds_end, cds_strand, read_length=150, scenario='correct_start'):
    """
    Generate a read mapping for a specific scenario.

    Scenarios:
    - 'correct_start': Read captures CDS start with correct prediction
    - 'correct_stop': Read captures CDS stop with correct prediction
    - 'correct_both': Read captures both start and stop
    - 'middle': Read is in the middle of CDS
    - 'left_overhang': Read overhangs left of CDS start
    - 'right_overhang': Read overhangs right of CDS end
    - 'poor_overlap': Read has insufficient overlap with CDS
    """

    if scenario == 'correct_start' or scenario == 'left_overhang':
        # Read starts before CDS and extends into it
        overhang = random.randint(10, 50)
        read_start = cds_start - overhang
        read_end = read_start + read_length - 1
        read_strand = '+'

    elif scenario == 'correct_stop' or scenario == 'right_overhang':
        # Read ends after CDS and starts within it
        overhang = random.randint(10, 50)
        read_end = cds_end + overhang
        read_start = read_end - read_length + 1
        read_strand = '+'

    elif scenario == 'correct_both':
        # Read is longer than CDS (rare but possible)
        read_length = (cds_end - cds_start + 1) + random.randint(20, 100)
        overhang = random.randint(10, 50)
        read_start = cds_start - overhang
        read_end = read_start + read_length - 1
        read_strand = '+'

    elif scenario == 'middle':
        # Read is entirely within CDS
        max_start = cds_end - read_length - 10
        min_start = cds_start + 10
        if max_start <= min_start:
            min_start = cds_start + 5
            max_start = cds_end - read_length - 5
        read_start = random.randint(min_start, max_start)
        read_end = read_start + read_length - 1
        read_strand = '+'

    elif scenario == 'poor_overlap':
        # Read has less than 60bp overlap with CDS
        overlap = random.randint(20, 55)
        read_start = cds_start - (read_length - overlap)
        read_end = read_start + read_length - 1
        read_strand = '+'

    else:
        # Default to middle
        read_start = random.randint(cds_start + 10, cds_end - read_length - 10)
        read_end = read_start + read_length - 1
        read_strand = '+'

    # Randomly flip to reverse strand
    if random.random() < 0.3:
        read_strand = '-'

    return read_start, read_end, read_strand


def generate_prediction(read_start, read_end, read_strand, cds_start, cds_end, cds_strand,
                        read_seq, scenario='correct_start', frame_offset=0):
    """
    Generate an ORF prediction on a read with explicit expected outcomes.

    Returns:
        tuple: (pred_start_on_read, pred_end_on_read, pred_strand, read_seq, expected_outcomes_dict)

    expected_outcomes_dict contains:
        - 'expected_answers': list of expected answer strings from check_pred
        - 'pred_scenario': description of what was generated
        - 'should_overlap_start': bool
        - 'should_overlap_stop': bool
        - 'is_middle_only': bool
    """

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

    # Determine if read captures CDS start/stop
    if category == ('+', '+'):
        read_captures_start = read_start <= cds_start
        read_captures_stop = read_end >= cds_end
    elif category == ('+', '-'):
        read_captures_start = read_start <= cds_start
        read_captures_stop = read_end >= cds_end
    elif category == ('-', '+'):
        read_captures_start = read_end >= cds_end
        read_captures_stop = read_start <= cds_start
    elif category == ('-', '-'):
        read_captures_start = read_end >= cds_end
        read_captures_stop = read_start <= cds_start
    else:
        read_captures_start = False
        read_captures_stop = False

    # CORRECT START scenarios
    if scenario == 'correct_start':
        if category == ('+', '+') and read_captures_start:
            pred_start_on_read = (cds_start - read_start + 1) + frame_offset
            pred_end_on_read = min(read_length, pred_start_on_read + 150)
            pred_strand = '+'

            # Insert ATG at prediction start
            pos = pred_start_on_read - 1
            if pos >= 0 and pos + 3 <= len(read_seq):
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

            return pred_start_on_read, pred_end_on_read, pred_strand, read_seq, expected

        elif category == ('-', '+') and read_captures_start:
            pred_end_on_read = cds_end - read_start + 1
            pred_start_on_read = max(1, pred_end_on_read - 150)
            pred_strand = '-'

            # Insert TAC (reverse comp of ATG) at prediction end
            pos = pred_end_on_read - 3
            if pos >= 0 and pos + 3 <= len(read_seq):
                read_seq = insert_codon_at_position(read_seq, pos, 'TAC')
                expected['pred_start_codon'] = 'ATG'

            expected['should_overlap_start'] = True
            expected['expected_answers'] = ['correct direction', 'correct start', 'correct frame']
            return pred_start_on_read, pred_end_on_read, pred_strand, read_seq, expected

    # CORRECT STOP scenarios
    elif scenario == 'correct_stop':
        if category == ('+', '+') and read_captures_stop:
            pred_end_on_read = (cds_end - read_start + 1) + frame_offset
            pred_start_on_read = max(1, pred_end_on_read - 150)
            pred_strand = '+'

            # Insert TAA at prediction end
            pos = pred_end_on_read - 3
            if pos >= 0 and pos + 3 <= len(read_seq):
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

            return pred_start_on_read, pred_end_on_read, pred_strand, read_seq, expected

        elif category == ('-', '+') and read_captures_stop:
            pred_start_on_read = cds_start - read_start + 1
            pred_end_on_read = min(read_length, pred_start_on_read + 150)
            pred_strand = '-'

            # Insert TTA (reverse comp of TAA) at prediction start
            pos = pred_start_on_read - 1
            if pos >= 0 and pos + 3 <= len(read_seq):
                read_seq = insert_codon_at_position(read_seq, pos, 'TTA')
                expected['pred_stop_codon'] = 'TAA'

            expected['should_overlap_stop'] = True
            expected['expected_answers'] = ['correct direction', 'correct stop', 'correct frame']
            return pred_start_on_read, pred_end_on_read, pred_strand, read_seq, expected

    # ALTERNATIVE START (in-frame but downstream)
    elif scenario == 'alternative_start':
        if category == ('+', '+') and read_captures_start:
            offset = random.choice([3, 6, 9])  # In-frame alternatives
            pred_start_on_read = (cds_start - read_start + 1) + offset
            pred_end_on_read = min(read_length, pred_start_on_read + 150)
            pred_strand = '+'

            # Insert GTG or TTG (alternative start codons)
            alt_codon = random.choice(['GTG', 'TTG'])
            pos = pred_start_on_read - 1
            if pos >= 0 and pos + 3 <= len(read_seq):
                read_seq = insert_codon_at_position(read_seq, pos, alt_codon)
                expected['pred_start_codon'] = alt_codon

            expected['should_overlap_start'] = True
            expected['expected_answers'] = ['correct direction', 'alternative start', 'correct frame']
            return pred_start_on_read, pred_end_on_read, pred_strand, read_seq, expected

    # INCORRECT START (wrong frame)
    elif scenario == 'incorrect_start':
        if category == ('+', '+') and read_captures_start:
            offset = random.choice([1, 2])  # Out of frame
            pred_start_on_read = (cds_start - read_start + 1) + offset
            pred_end_on_read = min(read_length, pred_start_on_read + 150)
            pred_strand = '+'

            expected['should_overlap_start'] = True
            expected['expected_answers'] = ['correct direction', 'incorrect start', 'incorrect frame']
            return pred_start_on_read, pred_end_on_read, pred_strand, read_seq, expected

    # MIDDLE PREDICTION
    elif scenario == 'middle':
        if not read_captures_start and not read_captures_stop:
            pred_start_on_read = random.randint(1, max(1, read_length - 100))
            pred_end_on_read = min(read_length, pred_start_on_read + 100)
            pred_strand = read_strand

            expected['is_middle_only'] = True
            expected['expected_answers'] = ['correct direction', 'middle']
            return pred_start_on_read, pred_end_on_read, pred_strand, read_seq, expected

    # WRONG DIRECTION
    elif scenario == 'wrong_direction':
        pred_start_on_read = random.randint(1, max(1, read_length - 100))
        pred_end_on_read = min(read_length, pred_start_on_read + 100)
        pred_strand = '-' if read_strand == '+' else '+'

        expected['expected_answers'] = ['incorrect direction']
        return pred_start_on_read, pred_end_on_read, pred_strand, read_seq, expected

    # Default fallback
    pred_start_on_read = random.randint(1, max(1, read_length - 60))
    pred_end_on_read = min(read_length, pred_start_on_read + random.randint(60, 150))
    pred_strand = read_strand
    expected['expected_answers'] = ['correct direction']
    return pred_start_on_read, pred_end_on_read, pred_strand, read_seq, expected


def write_intersect_bed(filename, entries):
    # Write the intersect BED file.
    with gzip.open(filename, 'wt') as f:
        # Write header
        f.write('\t'.join(['chrom', 'read_name', 'read_start', 'read_end', 'read_strand',
                           'score', 'feature_type', 'cds_start', 'cds_end', 'cds_strand',
                           'read_sequence']) + '\n')

        for entry in entries:
            f.write('\t'.join(map(str, entry)) + '\n')


def write_predictions_gff(filename, predictions):
    # Write the predictions GFF file.
    with gzip.open(filename, 'wt') as f:
        # Write GFF3 header
        f.write('##gff-version 3\n')

        for pred in predictions:
            f.write('\t'.join(map(str, pred)) + '\n')


def write_validation_csv(filename, validation_data):
    # Write validation data showing expected outcomes for each prediction.
    with open(filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['read_name',  'read_start', 'read_end', 'read_strand',
                         'pred_id', 'pred_start', 'pred_end', 'pred_strand',
                         'cds_id', 'cds_start', 'cds_end', 'cds_strand',
                         'pred_scenario', 'should_overlap_start', 'should_overlap_stop',
                         'is_middle_only', 'expected_answers', 'pred_start_codon', 'pred_stop_codon'])

        for entry in validation_data:
            writer.writerow(entry)


def main():
    parser = argparse.ArgumentParser(description='Generate synthetic test data for rORForise')
    parser.add_argument('-n', '--num_reads', type=int, default=1000,
                        help='Number of reads to generate (default: 1000)')
    parser.add_argument('-c', '--num_cds', type=int, default=50,
                        help='Number of CDS features to generate (default: 50)')
    parser.add_argument('-o', '--output_prefix', type=str, default='test_data',
                        help='Prefix for output files (default: test_data)')
    parser.add_argument('-g', '--genome_name', type=str, default='Chromosome-0',
                        help='Genome/chromosome name (default: Chromosome-0')
    parser.add_argument('-gc', '--gc_content', type=float, default=0.32,
                        help='GC content for sequence generation (default: 0.32)')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed for reproducibility (default: 42)')

    args = parser.parse_args()

    random.seed(args.seed)

    print(f"Generating test data with {args.num_reads} reads and {args.num_cds} CDS features...")

    # Generate CDS features
    cds_features = generate_cds_features(args.genome_name, num_cds=args.num_cds)
    print(f"Generated {len(cds_features)} CDS features")

    bed_entries = []
    gff_predictions = []
    validation_data = []

    # Define scenario distribution
    scenarios = [
        'correct_start', 'correct_stop', 'alternative_start',
        'incorrect_start', 'middle', 'wrong_direction'
    ]
    scenario_weights = [0.25, 0.25, 0.15, 0.15, 0.15, 0.05]

    scenario_counts = {s: 0 for s in scenarios}

    for read_idx in range(args.num_reads):
        # Randomly select a CDS
        cds_start, cds_end, cds_strand = random.choice(cds_features)

        # Select scenario
        pred_scenario = random.choices(scenarios, weights=scenario_weights)[0]
        scenario_counts[pred_scenario] += 1

        # Determine read mapping type based on prediction scenario
        if pred_scenario in ['correct_start', 'alternative_start', 'incorrect_start']:
            mapping_scenario = 'left_overhang'
        elif pred_scenario == 'correct_stop':
            mapping_scenario = 'right_overhang'
        elif pred_scenario == 'middle':
            mapping_scenario = 'middle'
        else:
            mapping_scenario = random.choice(['left_overhang', 'right_overhang', 'middle'])

        # Generate read mapping
        read_length = random.randint(100, 250)
        read_start, read_end, read_strand = generate_read_mapping(
            cds_start, cds_end, cds_strand, read_length, mapping_scenario
        )

        read_name = f"read_{read_idx:06d}"

        # Generate read sequence
        read_seq_length = read_end - read_start + 1
        read_seq = generate_sequence(read_seq_length, args.gc_content)

        # Generate prediction with expected outcomes
        frame_offset = 0 if pred_scenario != 'incorrect_start' else random.choice([1, 2])

        pred_start, pred_end, pred_strand, modified_seq, expected = generate_prediction(
            read_start, read_end, read_strand, cds_start, cds_end, cds_strand,
            read_seq, pred_scenario, frame_offset
        )

        # Update read sequence with modifications
        read_seq = modified_seq

        # If read is reverse strand, the sequence in the file should be reverse complemented
        if read_strand == '-':
            read_seq = reverse_complement(read_seq)

        # Add to BED entries
        bed_entry = [
            args.genome_name,  # chrom
            read_name,  # read_name
            read_start,  # read_start
            read_end,  # read_end
            read_strand,  # read_strand
            '.',  # score
            'CDS',  # feature_type
            cds_start,  # cds_start
            cds_end,  # cds_end
            cds_strand,  # cds_strand
            read_seq  # read_sequence
        ]
        bed_entries.append(bed_entry)

        # Ensure prediction is at least 60bp
        if pred_end - pred_start + 1 >= 60:
            pred_id = f"{read_name}_pred_0"
            gff_entry = [
                f"@{read_name}",  # sequence name
                'TestGenerator',  # source
                'CDS',  # feature
                pred_start,  # start (1-based)
                pred_end,  # end (1-based)
                '.',  # score
                pred_strand,  # strand
                '0',  # frame
                f"ID={pred_id}"  # attributes
            ]
            gff_predictions.append(gff_entry)

            # Add validation entry
            validation_entry = [
                read_name,  # read_name
                read_start,  # read_start
                read_end,  # read_end
                read_strand,  # read_strand
                pred_id,  # pred_id
                pred_start,  # pred_start
                pred_end,  # pred_end
                pred_strand,  # pred_strand
                f"cds_{cds_start}_{cds_end}",  # cds_id
                cds_start,  # cds_start
                cds_end,  # cds_end
                cds_strand,  # cds_strand
                expected['pred_scenario'],  # pred_scenario
                expected['should_overlap_start'],  # should_overlap_start
                expected['should_overlap_stop'],  # should_overlap_stop
                expected['is_middle_only'],  # is_middle_only
                '; '.join(expected['expected_answers']),  # expected_answers
                expected.get('pred_start_codon', ''),  # pred_start_codon
                expected.get('pred_stop_codon', '')  # pred_stop_codon
            ]
            validation_data.append(validation_entry)

    # Write output files
    bed_file = f"{args.output_prefix}_intersect.bed.gz"
    gff_file = f"{args.output_prefix}_predictions.gff.gz"
    validation_file = f"{args.output_prefix}_validation.csv"

    write_intersect_bed(bed_file, bed_entries)
    write_predictions_gff(gff_file, gff_predictions)
    write_validation_csv(validation_file, validation_data)

    print(f"\nGenerated files:")
    print(f"  BED file: {bed_file} ({len(bed_entries)} entries)")
    print(f"  GFF file: {gff_file} ({len(gff_predictions)} predictions)")
    print(f"  Validation file: {validation_file} ({len(validation_data)} entries)")

    print(f"\nScenario distribution:")
    for scenario, count in scenario_counts.items():
        print(f"  {scenario}: {count} ({count / args.num_reads * 100:.1f}%)")

    print(f"\nYou can now run rORForise with:")
    print(f"  rORForise -int_bed {bed_file} -p_gff {gff_file} \\")
    print(f"    -o output_dir -gc_prob {args.gc_content}")
    print(f"\nThen compare the detailed_results.csv against {validation_file} to validate correctness")


if __name__ == '__main__':
    main()