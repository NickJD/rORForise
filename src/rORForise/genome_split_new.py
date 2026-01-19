
import argparse
import random
import logging
from pathlib import Path
from collections import defaultdict, namedtuple
import sys

try:
    from .constants import *
except (ModuleNotFoundError, ImportError, NameError, TypeError):
    try:
        from constants import *
    except:
        rORForise_VERSION = "1.0.0"  # Fallback if constants not available

# Data structures for better organization
CDSAnnotation = namedtuple('CDSAnnotation', ['chrom', 'start', 'end', 'strand', 'gene_id', 'attributes'])
SyntheticRead = namedtuple('SyntheticRead', ['name', 'sequence', 'start', 'end', 'strand', 'source_info'])


def setup_logging(verbose=False):
    """Set up logging configuration"""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s',
        level=level
    )
    return logging.getLogger(__name__)


def parse_fasta(fasta_file):
    """
    Parse FASTA file handling multiple chromosomes/contigs
    Returns dict: {chrom_name: sequence}
    """
    logger = logging.getLogger(__name__)
    sequences = {}
    current_chrom = None
    current_seq = []

    try:
        with open(fasta_file, 'r') as file:
            for line_num, line in enumerate(file, 1):
                line = line.strip()
                if line.startswith('>'):
                    # Save previous sequence
                    if current_chrom and current_seq:
                        sequences[current_chrom] = ''.join(current_seq)
                        logger.debug(f"Loaded chromosome {current_chrom}: {len(sequences[current_chrom])} bp")

                    # Start new sequence
                    current_chrom = line[1:].split()[0]  # Take first word after >
                    current_seq = []
                elif line and not line.startswith('#'):
                    # Remove any non-DNA characters and uppercase
                    clean_seq = ''.join(c.upper() for c in line if c.upper() in 'ATCGN')
                    current_seq.append(clean_seq)

        # Don't forget the last sequence
        if current_chrom and current_seq:
            sequences[current_chrom] = ''.join(current_seq)
            logger.debug(f"Loaded chromosome {current_chrom}: {len(sequences[current_chrom])} bp")

    except Exception as e:
        logger.error(f"Error parsing FASTA file {fasta_file}: {e}")
        raise

    logger.info(f"Loaded {len(sequences)} chromosome(s), total: {sum(len(s) for s in sequences.values())} bp")
    return sequences


def parse_gff_for_cds(gff_file):
    """
    Enhanced GFF parser that extracts comprehensive CDS information
    Returns list of CDSAnnotation objects
    """
    logger = logging.getLogger(__name__)
    cds_annotations = []

    try:
        with open(gff_file, 'r') as file:
            for line_num, line in enumerate(file, 1):
                line = line.strip()
                if line.startswith('#') or not line:
                    continue

                try:
                    parts = line.split('\t')
                    if len(parts) < 9:
                        continue

                    chrom, source, feature_type, start, end, score, strand, phase, attributes = parts

                    if feature_type.upper() == 'CDS':
                        # Parse attributes for gene information
                        attr_dict = {}
                        for attr in attributes.split(';'):
                            if '=' in attr:
                                key, value = attr.split('=', 1)
                                attr_dict[key.strip()] = value.strip()

                        gene_id = attr_dict.get('ID', attr_dict.get('gene_id', f'gene_{len(cds_annotations)}'))

                        cds_annotations.append(CDSAnnotation(
                            chrom=chrom,
                            start=int(start),
                            end=int(end),
                            strand=strand,
                            gene_id=gene_id,
                            attributes=attr_dict
                        ))

                except (ValueError, IndexError) as e:
                    logger.warning(f"Skipping malformed GFF line {line_num}: {e}")
                    continue

    except Exception as e:
        logger.error(f"Error parsing GFF file {gff_file}: {e}")
        raise

    logger.info(f"Loaded {len(cds_annotations)} CDS annotations")
    return cds_annotations


def reverse_complement(seq):
    """Generate reverse complement of DNA sequence"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))


def introduce_sequencing_errors(sequence, error_rate=0.01):
    """
    Introduce realistic sequencing errors (substitutions, indels)
    """
    if error_rate <= 0:
        return sequence

    seq_list = list(sequence)
    bases = ['A', 'T', 'C', 'G']

    for i in range(len(seq_list)):
        if random.random() < error_rate:
            error_type = random.choice(['substitution', 'deletion', 'insertion'])

            if error_type == 'substitution':
                # Replace with different base
                current_base = seq_list[i]
                new_base = random.choice([b for b in bases if b != current_base])
                seq_list[i] = new_base
            elif error_type == 'deletion' and len(seq_list) > 1:
                # Delete this base
                seq_list[i] = ''
            elif error_type == 'insertion':
                # Insert random base
                seq_list[i] = seq_list[i] + random.choice(bases)

    return ''.join(seq_list)


def generate_synthetic_reads(sequences, cds_annotations, read_length=150, coverage=10,
                             overlap_min=20, error_rate=0.0, include_non_coding=False,
                             fragment_size_mean=300, fragment_size_std=50):
    """
    Generate synthetic reads with various sampling strategies
    """
    logger = logging.getLogger(__name__)
    synthetic_reads = []
    read_counter = 0

    # Group CDS by chromosome
    cds_by_chrom = defaultdict(list)
    for cds in cds_annotations:
        cds_by_chrom[cds.chrom].append(cds)

    # Generate reads for each chromosome
    for chrom, seq in sequences.items():
        chrom_cds = cds_by_chrom.get(chrom, [])
        logger.info(f"Generating reads for {chrom}: {len(seq)} bp, {len(chrom_cds)} CDS regions")

        # Strategy 1: CDS-focused reads (ensuring good CDS coverage)
        for cds in chrom_cds:
            cds_length = cds.end - cds.start + 1
            reads_needed = max(1, int(cds_length * coverage / read_length))

            for _ in range(reads_needed):
                # Random position within or around CDS
                margin = read_length // 2
                start_pos = random.randint(
                    max(0, cds.start - margin),
                    min(len(seq) - read_length, cds.end + margin - read_length)
                )
                end_pos = start_pos + read_length

                if end_pos > len(seq):
                    continue

                read_seq = seq[start_pos:end_pos]
                read_strand = random.choice(['+', '-']) if cds.strand == '.' else cds.strand

                if read_strand == '-':
                    read_seq = reverse_complement(read_seq)

                # Add sequencing errors
                if error_rate > 0:
                    read_seq = introduce_sequencing_errors(read_seq, error_rate)

                read_name = f"synthetic_read_{read_counter:08d}_cds_{cds.gene_id}"
                source_info = {
                    'source': 'cds_focused',
                    'target_cds': cds.gene_id,
                    'overlap_length': min(end_pos, cds.end) - max(start_pos, cds.start - 1)
                }

                synthetic_reads.append(SyntheticRead(
                    name=read_name,
                    sequence=read_seq,
                    start=start_pos + 1,  # 1-based
                    end=end_pos,  # 1-based
                    strand=read_strand,
                    source_info=source_info
                ))
                read_counter += 1

        # Strategy 2: Random genome-wide reads (if requested)
        if include_non_coding:
            genome_reads_needed = int(len(seq) * coverage / read_length / 10)  # 10x less than CDS

            for _ in range(genome_reads_needed):
                start_pos = random.randint(0, len(seq) - read_length)
                end_pos = start_pos + read_length

                read_seq = seq[start_pos:end_pos]
                read_strand = random.choice(['+', '-'])

                if read_strand == '-':
                    read_seq = reverse_complement(read_seq)

                if error_rate > 0:
                    read_seq = introduce_sequencing_errors(read_seq, error_rate)

                read_name = f"synthetic_read_{read_counter:08d}_random"
                source_info = {'source': 'random_genome'}

                synthetic_reads.append(SyntheticRead(
                    name=read_name,
                    sequence=read_seq,
                    start=start_pos + 1,
                    end=end_pos,
                    strand=read_strand,
                    source_info=source_info
                ))
                read_counter += 1

        # Strategy 3: Overlapping reads (sliding window)
        step_size = read_length - overlap_min
        for start_pos in range(0, len(seq) - read_length + 1, step_size):
            end_pos = start_pos + read_length

            # Only create if overlaps with CDS
            overlaps_cds = any(
                cds.start <= end_pos and cds.end >= start_pos + 1
                for cds in chrom_cds
            )

            if overlaps_cds:
                read_seq = seq[start_pos:end_pos]
                read_strand = '+'  # Keep original orientation for sliding window

                if error_rate > 0:
                    read_seq = introduce_sequencing_errors(read_seq, error_rate)

                read_name = f"synthetic_read_{read_counter:08d}_sliding"
                source_info = {'source': 'sliding_window'}

                synthetic_reads.append(SyntheticRead(
                    name=read_name,
                    sequence=read_seq,
                    start=start_pos + 1,
                    end=end_pos,
                    strand=read_strand,
                    source_info=source_info
                ))
                read_counter += 1

    logger.info(f"Generated {len(synthetic_reads)} synthetic reads")
    return synthetic_reads


def write_output(synthetic_reads, cds_annotations, output_fasta, output_bed,
                 output_mapping=None, detailed_bed=False):
    """
    Write outputs with comprehensive annotation tracking
    """
    logger = logging.getLogger(__name__)

    # Write FASTA
    with open(output_fasta, "w") as fasta_out:
        for read in synthetic_reads:
            fasta_out.write(f">{read.name}\n{read.sequence}\n")

    # Write BED with read-CDS overlaps
    with open(output_bed, "w") as bed_out:
        if detailed_bed:
            bed_out.write(
                "# chrom\tstart\tend\tread_name\tscore\tstrand\tcds_start\tcds_end\tcds_strand\tgene_id\toverlap_bp\tsource\n")

        for read in synthetic_reads:
            # Find overlapping CDS regions
            read_start_0based = read.start - 1  # Convert to 0-based for calculations
            read_end_0based = read.end

            overlapping_cds = []
            for cds in cds_annotations:
                if (cds.start <= read.end and cds.end >= read.start):
                    overlap_start = max(read.start, cds.start)
                    overlap_end = min(read.end, cds.end)
                    overlap_length = overlap_end - overlap_start + 1

                    if overlap_length > 0:
                        overlapping_cds.append((cds, overlap_length))

            # Write BED entry for each overlapping CDS
            if overlapping_cds:
                for cds, overlap_length in overlapping_cds:
                    if detailed_bed:
                        bed_out.write(
                            f"{cds.chrom}\t{read_start_0based}\t{read_end_0based}\t{read.name}\t{overlap_length}\t{read.strand}\t{cds.start}\t{cds.end}\t{cds.strand}\t{cds.gene_id}\t{overlap_length}\t{read.source_info.get('source', 'unknown')}\n")
                    else:
                        bed_out.write(
                            f"{cds.chrom}\t{read_start_0based}\t{read_end_0based}\t{read.name}\t{overlap_length}\t{read.strand}\n")
            else:
                # Non-coding read
                if detailed_bed:
                    bed_out.write(
                        f"unknown\t{read_start_0based}\t{read_end_0based}\t{read.name}\t0\t{read.strand}\t.\t.\t.\t.\t0\t{read.source_info.get('source', 'unknown')}\n")

    # Write mapping file (read -> CDS relationships)
    if output_mapping:
        with open(output_mapping, "w") as map_out:
            map_out.write(
                "read_name\tread_start\tread_end\tread_strand\tcds_gene_id\tcds_start\tcds_end\tcds_strand\toverlap_length\tsource_type\n")

            for read in synthetic_reads:
                overlapping_cds = []
                for cds in cds_annotations:
                    if (cds.start <= read.end and cds.end >= read.start):
                        overlap_start = max(read.start, cds.start)
                        overlap_end = min(read.end, cds.end)
                        overlap_length = overlap_end - overlap_start + 1

                        if overlap_length > 0:
                            overlapping_cds.append((cds, overlap_length))

                if overlapping_cds:
                    for cds, overlap_length in overlapping_cds:
                        map_out.write(
                            f"{read.name}\t{read.start}\t{read.end}\t{read.strand}\t{cds.gene_id}\t{cds.start}\t{cds.end}\t{cds.strand}\t{overlap_length}\t{read.source_info.get('source')}\n")
                else:
                    map_out.write(
                        f"{read.name}\t{read.start}\t{read.end}\t{read.strand}\tno_overlap\t.\t.\t.\t0\t{read.source_info.get('source')}\n")

    logger.info(f"Output written to: {output_fasta}, {output_bed}")
    if output_mapping:
        logger.info(f"Mapping file written to: {output_mapping}")


def generate_statistics_report(synthetic_reads, cds_annotations, output_stats):
    """Generate comprehensive statistics about the synthetic dataset"""
    logger = logging.getLogger(__name__)

    stats = {
        'total_reads': len(synthetic_reads),
        'read_length_mean': sum(len(r.sequence) for r in synthetic_reads) / len(synthetic_reads),
        'source_distribution': defaultdict(int),
        'strand_distribution': defaultdict(int),
        'cds_coverage': defaultdict(int),
        'total_cds': len(cds_annotations)
    }

    for read in synthetic_reads:
        stats['source_distribution'][read.source_info.get('source', 'unknown')] += 1
        stats['strand_distribution'][read.strand] += 1

    # Calculate CDS coverage
    for read in synthetic_reads:
        for cds in cds_annotations:
            if (cds.start <= read.end and cds.end >= read.start):
                stats['cds_coverage'][cds.gene_id] += 1

    with open(output_stats, 'w') as f:
        f.write("# Synthetic Read Generation Statistics\n")
        f.write(f"Total reads generated: {stats['total_reads']}\n")
        f.write(f"Average read length: {stats['read_length_mean']:.1f}\n")
        f.write(f"Total CDS regions: {stats['total_cds']}\n")
        f.write(f"CDS regions with coverage: {len(stats['cds_coverage'])}\n")
        f.write(f"Coverage percentage: {len(stats['cds_coverage']) / stats['total_cds'] * 100:.1f}%\n\n")

        f.write("Source distribution:\n")
        for source, count in stats['source_distribution'].items():
            f.write(f"  {source}: {count} ({count / stats['total_reads'] * 100:.1f}%)\n")

        f.write("\nStrand distribution:\n")
        for strand, count in stats['strand_distribution'].items():
            f.write(f"  {strand}: {count} ({count / stats['total_reads'] * 100:.1f}%)\n")

    logger.info(f"Statistics report written to: {output_stats}")


def main():
    parser = argparse.ArgumentParser(
        description=f"Enhanced synthetic read generator v{getattr(sys.modules[__name__], 'rORForise_VERSION', '1.0.0')} - "
                    "Generate realistic synthetic reads for validating gene prediction algorithms.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Input files
    parser.add_argument('-f', '--fasta_file', type=Path, required=True,
                        help="Path to input FASTA file(s) containing genome sequence(s)")
    parser.add_argument('-g', '--gff_file', type=Path, required=True,
                        help="Path to input GFF file containing CDS annotations")

    # Read generation parameters
    parser.add_argument('-l', '--read_length', type=int, default=150,
                        help="Length of synthetic reads to generate")
    parser.add_argument('-c', '--coverage', type=float, default=10.0,
                        help="Target coverage depth for CDS regions")
    parser.add_argument('-o', '--overlap_min', type=int, default=20,
                        help="Minimum overlap between consecutive reads in sliding window mode")

    # Advanced options
    parser.add_argument('--error_rate', type=float, default=0.0,
                        help="Sequencing error rate (0.0-1.0)")
    parser.add_argument('--include_non_coding', action='store_true',
                        help="Include random non-coding reads")
    parser.add_argument('--fragment_size_mean', type=int, default=300,
                        help="Mean fragment size for paired-end simulation")
    parser.add_argument('--fragment_size_std', type=int, default=50,
                        help="Standard deviation of fragment size")

    # Output files
    parser.add_argument('--output_fasta', type=Path, default="synthetic_reads.fasta",
                        help="Output FASTA file for synthetic reads")
    parser.add_argument('--output_bed', type=Path, default="synthetic_annotations.bed",
                        help="Output BED file with read-CDS overlaps")
    parser.add_argument('--output_mapping', type=Path, default=None,
                        help="Output detailed read-to-CDS mapping file")
    parser.add_argument('--output_stats', type=Path, default="generation_stats.txt",
                        help="Output statistics report")

    # Options
    parser.add_argument('--detailed_bed', action='store_true',
                        help="Generate detailed BED file with additional annotation columns")
    parser.add_argument('--seed', type=int, default=None,
                        help="Random seed for reproducible results")
    parser.add_argument('-v', '--verbose', action='store_true',
                        help="Enable verbose logging")

    options = parser.parse_args()

    # Setup
    logger = setup_logging(options.verbose)
    if options.seed is not None:
        random.seed(options.seed)
        logger.info(f"Random seed set to: {options.seed}")

    try:
        # Parse inputs
        logger.info("Parsing input files...")
        genome_sequences = parse_fasta(options.fasta_file)
        cds_annotations = parse_gff_for_cds(options.gff_file)

        if not genome_sequences:
            raise ValueError("No sequences found in FASTA file")
        if not cds_annotations:
            raise ValueError("No CDS annotations found in GFF file")

        # Generate synthetic reads
        logger.info("Generating synthetic reads...")
        synthetic_reads = generate_synthetic_reads(
            genome_sequences, cds_annotations,
            read_length=options.read_length,
            coverage=options.coverage,
            overlap_min=options.overlap_min,
            error_rate=options.error_rate,
            include_non_coding=options.include_non_coding,
            fragment_size_mean=options.fragment_size_mean,
            fragment_size_std=options.fragment_size_std
        )

        # Write outputs
        logger.info("Writing output files...")
        write_output(
            synthetic_reads, cds_annotations,
            options.output_fasta, options.output_bed,
            options.output_mapping, options.detailed_bed
        )

        # Generate statistics
        generate_statistics_report(synthetic_reads, cds_annotations, options.output_stats)

        logger.info("Synthetic read generation completed successfully!")

    except Exception as e:
        logger.error(f"Error during execution: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())