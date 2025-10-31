#!/usr/bin/env python3
"""
Fragment a genome into read-like chunks and generate realistic alignment files.

Creates:
1. FASTA file with fragmented reads
2. BED file simulating genome alignments (similar to bedtools intersect output)
3. Optional SAM file for more realistic simulation

This is useful for testing ORF prediction tools on simulated sequencing data.
"""

import argparse
import random
import gzip
from pathlib import Path


def parse_fasta(fasta_file):
    """Parse FASTA file and return sequence with header."""
    sequences = {}
    current_header = None
    current_seq = []

    with open(fasta_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if current_header:
                    sequences[current_header] = ''.join(current_seq)
                current_header = line[1:].split()[0]  # Take first part of header
                current_seq = []
            else:
                current_seq.append(line)

        if current_header:
            sequences[current_header] = ''.join(current_seq)

    return sequences


def parse_gff_for_cds(gff_file):
    """
    Parse GFF file and extract CDS annotations.

    Returns list of tuples: (chrom, start, end, strand, attributes)
    """
    cds_annotations = []
    with open(gff_file, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                parts = line.strip().split('\t')
                if len(parts) >= 9 and parts[2] == 'CDS':
                    chrom = parts[0]
                    start = int(parts[3])
                    end = int(parts[4])
                    strand = parts[6]
                    attributes = parts[8] if len(parts) > 8 else ''
                    cds_annotations.append((chrom, start, end, strand, attributes))
    return cds_annotations


def reverse_complement(seq):
    # Return the reverse complement of a DNA sequence.
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g', 'n': 'n'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))


def introduce_sequencing_errors(sequence, error_rate=0.01):
    # Introduce random sequencing errors into a sequence.
    if error_rate <= 0:
        return sequence

    seq_list = list(sequence)
    bases = ['A', 'T', 'G', 'C']

    for i in range(len(seq_list)):
        if random.random() < error_rate:
            # Choose error type: substitution (0.7), deletion (0.15), insertion (0.15)
            error_type = random.choices(['sub', 'del', 'ins'], weights=[0.7, 0.15, 0.15])[0]

            if error_type == 'sub':
                # Substitution
                original = seq_list[i].upper()
                if original in bases:
                    seq_list[i] = random.choice([b for b in bases if b != original])
            elif error_type == 'del':
                # Deletion
                seq_list[i] = ''
            else:
                # Insertion
                seq_list[i] = seq_list[i] + random.choice(bases)

    return ''.join(seq_list)


def split_genome(sequences, cds_annotations, chunk_size=150, overlap=20,
                 reverse_prob=0.3, error_rate=0.0):
    # Split genome into overlapping chunks simulating sequencing reads.
    chunks = []
    read_counter = 0

    for chrom_name, sequence in sequences.items():
        # Get CDS annotations for this chromosome
        chrom_cds = [(start, end, strand, attrs) for c, start, end, strand, attrs
                     in cds_annotations if c == chrom_name]

        # Generate chunks with sliding window
        step_size = chunk_size - overlap
        for i in range(0, len(sequence), step_size):
            chunk_end = min(i + chunk_size, len(sequence))

            # Skip if chunk is too short
            if chunk_end - i < chunk_size // 2:
                continue

            chunk_start = i
            chunk_seq = sequence[chunk_start:chunk_end]

            # Determine strand (simulate forward/reverse reads)
            read_strand = '-' if random.random() < reverse_prob else '+'

            # If reverse strand, reverse complement the sequence
            if read_strand == '-':
                chunk_seq = reverse_complement(chunk_seq)

            # Introduce sequencing errors if requested
            if error_rate > 0:
                chunk_seq = introduce_sequencing_errors(chunk_seq, error_rate)

            # Find overlapping CDS annotations
            overlapping_cds = []
            for cds_start, cds_end, cds_strand, cds_attrs in chrom_cds:
                # Check if chunk overlaps with CDS
                if chunk_start < cds_end and chunk_end > cds_start:
                    overlapping_cds.append((cds_start, cds_end, cds_strand, cds_attrs))

            read_name = f"read_{read_counter:06d}"

            # Genomic coordinates are 1-based in GFF, but we use 0-based internally
            # and convert to 1-based when writing
            chunks.append((
                chrom_name,
                read_name,
                chunk_start + 1,  # Convert to 1-based
                chunk_end,  # 1-based end (inclusive in GFF format)
                read_strand,
                chunk_seq,
                overlapping_cds
            ))

            read_counter += 1

    return chunks


def write_output(options, chunks, compress_output=False):
    # Write output FASTA and BED files.
    # Determine if we should compress
    fasta_mode = 'wt' if not compress_output else 'wt'
    bed_mode = 'wt' if not compress_output else 'wt'

    fasta_file = options.output_fasta if not compress_output else f"{options.output_fasta}.gz"
    bed_file = options.output_bed if not compress_output else f"{options.output_bed}.gz"

    # Open files (with or without compression)
    if compress_output:
        fasta_out = gzip.open(fasta_file, fasta_mode)
        bed_out = gzip.open(bed_file, bed_mode)
    else:
        fasta_out = open(fasta_file, fasta_mode)
        bed_out = open(bed_file, bed_mode)

    try:
        # Write BED header (matching bedtools intersect output format)
        bed_out.write('\t'.join([
            'chrom', 'read_name', 'read_start', 'read_end', 'read_strand',
            'score', 'feature_type', 'cds_start', 'cds_end', 'cds_strand',
            'read_sequence'
        ]) + '\n')

        for chrom, read_name, read_start, read_end, read_strand, read_seq, cds_list in chunks:
            # Write to FASTA
            fasta_out.write(f">{read_name}\n{read_seq}\n")

            # Write to BED (one line per CDS overlap)
            if cds_list:
                for cds_start, cds_end, cds_strand, cds_attrs in cds_list:
                    # BED format matching your intersect file structure
                    bed_out.write('\t'.join([
                        str(chrom),  # chrom
                        str(read_name),  # read_name
                        str(read_start),  # read_start (1-based)
                        str(read_end),  # read_end (1-based)
                        str(read_strand),  # read_strand
                        '.',  # score
                        'CDS',  # feature_type
                        str(cds_start),  # cds_start (1-based)
                        str(cds_end),  # cds_end (1-based)
                        str(cds_strand),  # cds_strand
                        str(read_seq)  # read_sequence
                    ]) + '\n')
            elif options.non_coding:
                # Read doesn't overlap any CDS - optionally include these
                bed_out.write('\t'.join([
                    str(chrom), str(read_name), str(read_start), str(read_end),
                    str(read_strand), '.', 'non-coding', '.', '.', '.', str(read_seq)
                ]) + '\n')
                pass

    finally:
        fasta_out.close()
        bed_out.close()

    print(f"Written {len(chunks)} reads to {fasta_file}")
    print(f"Written BED file to {bed_file}")


def write_statistics(chunks, output_stats):
    # Write statistics about the fragmentation.

    total_reads = len(chunks)
    reads_with_cds = sum(1 for c in chunks if c[6])  # c[6] is overlapping_cds list
    total_cds_overlaps = sum(len(c[6]) for c in chunks)

    strand_counts = {'+': 0, '-': 0}
    for chunk in chunks:
        strand_counts[chunk[4]] += 1

    read_lengths = [len(c[5]) for c in chunks]  # c[5] is read_seq

    with open(output_stats, 'w') as f:
        f.write("Fragmentation Statistics\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Total reads generated: {total_reads}\n")
        f.write(f"Reads overlapping CDS: {reads_with_cds} ({reads_with_cds / total_reads * 100:.1f}%)\n")
        f.write(f"Total CDS overlaps: {total_cds_overlaps}\n")
        f.write(f"Average CDS overlaps per read: {total_cds_overlaps / total_reads:.2f}\n")
        f.write(f"\nStrand distribution:\n")
        f.write(f"  Forward (+): {strand_counts['+']} ({strand_counts['+'] / total_reads * 100:.1f}%)\n")
        f.write(f"  Reverse (-): {strand_counts['-']} ({strand_counts['-'] / total_reads * 100:.1f}%)\n")
        f.write(f"\nRead length statistics:\n")
        f.write(f"  Min: {min(read_lengths)}\n")
        f.write(f"  Max: {max(read_lengths)}\n")
        f.write(f"  Mean: {sum(read_lengths) / len(read_lengths):.1f}\n")

    print(f"Written statistics to {output_stats}")


def main():
    parser = argparse.ArgumentParser(
        description="Fragment a genome FASTA file into chunks and generate realistic alignment files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python fragment_genome.py -fasta genome.fa -gff annotations.gff -output_fasta reads.fa -output_bed alignments.bed

  # With custom parameters and statistics
  python fragment_genome.py -fasta genome.fa -gff annotations.gff -chunk_size 200 -overlap 30 \\
      -reverse_prob 0.5 -error_rate 0.01 -output_fasta reads.fa -output_bed alignments.bed \\
      -output_stats fragmentation_stats.txt

  # With compression
  python fragment_genome.py -fasta genome.fa -gff annotations.gff -output_fasta reads.fa \\
      -output_bed alignments.bed --compress
        """
    )

    parser.add_argument(
        "-fasta", dest="fasta", required=True,
        help="Path to the input FASTA file."
    )
    parser.add_argument(
        "-gff", dest="gff", required=True,
        help="Path to the input GFF file with CDS annotations."
    )
    parser.add_argument(
            "-non_coding", dest="non_coding", default=False, action='store_true',
            help="Include non-coding regions in the fragmentation process (default: False)."
    )
    parser.add_argument(
        "-chunk_size", dest="chunk_size", type=int, default=150,
        help="Size of each genome chunk/read length (default: 150)."
    )
    parser.add_argument(
        "-overlap", dest="overlap", type=int, default=20,
        help="Number of overlapping bases between chunks (default: 20)."
    )
    parser.add_argument(
        "-reverse_prob", dest="reverse_prob", type=float, default=0.3,
        help="Probability of generating reverse-strand reads (default: 0.3)."
    )
    parser.add_argument(
        "-error_rate", dest="error_rate", type=float, default=0.0,
        help="Sequencing error rate to introduce (default: 0.0, range: 0.0-0.1)."
    )
    parser.add_argument(
        "-output_fasta", dest="output_fasta", required=True,
        help="Output FASTA file for genome chunks/reads."
    )
    parser.add_argument(
        "-output_bed", dest="output_bed", required=True,
        help="Output BED file (bedtools intersect format) for read-CDS mappings."
    )
    parser.add_argument(
        "-output_stats", dest="output_stats",
        help="Optional output file for fragmentation statistics."
    )
    parser.add_argument(
        "--compress", action="store_true",
        help="Compress output files with gzip."
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed for reproducibility (default: 42)."
    )

    options = parser.parse_args()

    # Set random seed
    random.seed(options.seed)

    # Validate error rate
    if not 0.0 <= options.error_rate <= 0.1:
        parser.error("Error rate must be between 0.0 and 0.1")

    if not 0.0 <= options.reverse_prob <= 1.0:
        parser.error("Reverse probability must be between 0.0 and 1.0")

    print("Parsing genome FASTA file...")
    sequences = parse_fasta(options.fasta)
    print(f"Loaded {len(sequences)} sequence(s)")

    print("Parsing GFF file...")
    annotations = parse_gff_for_cds(options.gff)
    print(f"Found {len(annotations)} CDS annotations")

    print(f"Fragmenting genome with chunk_size={options.chunk_size}, overlap={options.overlap}...")
    chunks = split_genome(
        sequences,
        annotations,
        options.chunk_size,
        options.overlap,
        options.reverse_prob,
        options.error_rate
    )

    # Use first chromosome name as genome name
    #genome_name = list(sequences.keys())[0] if sequences else "genome"

    print("Writing output files...")
    write_output(options, chunks)

    if options.output_stats:
        write_statistics(chunks, options.output_stats)

    print("\nFragmentation complete!")


if __name__ == "__main__":
    main()