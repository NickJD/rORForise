import argparse
import gzip

try: # Try to import from the package if available
    from .adapters import parse_gff_attributes
    from .utils import reverse_complement
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from adapters import parse_gff_attributes
    from utils import reverse_complement


def _open_text(path):
    return gzip.open(path, 'rt', encoding='utf-8') if str(path).endswith('.gz') else open(path, 'r', encoding='utf-8')


def _load_pysam():
    try:
        import pysam
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "bam-gff-intersect requires pysam. Install it with `pip install rORForise[bam]` "
            "or install pysam in the active environment."
        ) from exc
    return pysam

def parse_gff(gff_file):
    # Parse GFF file and return a list of tuples.
    gene_list = []
    with _open_text(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            cols = line.strip().split('\t')
            if len(cols) < 9:
                continue
            seqid = cols[0]
            source = cols[1]
            type = cols[2]
            start = int(cols[3])
            end = int(cols[4])
            score = cols[5]
            strand = cols[6]
            phase = cols[7]
            attributes = parse_gff_attributes(cols[8])
            gene_id = attributes.get('ID') or attributes.get('Name') or f"{seqid}:{start}-{end}"
            gene_list.append((seqid, source, type, start, end, score, strand, phase, gene_id))
    return gene_list

def process_bam_and_gff(bam_file, gff_file, features, output_file):
    # Process BAM and GFF files and output information in a tab-separated file.
    # Parse GFF file
    gff_data = parse_gff(gff_file)
    feature_filter = set(features.split(',')) if features else None

    # Open BAM file
    pysam = _load_pysam()
    bam = pysam.AlignmentFile(bam_file, "rb")

    # Open output file
    with open(output_file, 'w') as out:
        # Write header
        out.write("Chromosome\tReadID\tStart\tEnd\tDirection\tMappingQuality\tFeatureType\tGeneStart\tGeneEnd\tGeneStrand\tReadSequence\n")

        # Process each read in BAM file
        for read in bam.fetch():
            if read.is_unmapped:
                continue
            read_id = read.query_name
            if read.is_paired:
                if read.is_read1:
                    read_id = f"{read_id}/1"
                elif read.is_read2:
                    read_id = f"{read_id}/2"
            else:
                read_id = read.query_name
            chrom = read.reference_name
            start = read.reference_start +1 # bam file starts are base-0
            end = read.reference_end
            read_strand = '-' if read.is_reverse else '+'
            mapq = read.mapping_quality
            sequence = read.query_sequence or ''
            if read_strand == '-':
                sequence = reverse_complement(sequence)

            # Find corresponding gene in GFF data

            for gene in gff_data:
                (seqid, source, type, gene_start, gene_end, score, strand, phase, gene_id) = gene
                if feature_filter is not None and type not in feature_filter:
                    continue
                if seqid == chrom and min(end, gene_end) - max(start, gene_start) + 1 > 0:
                    gene_strand = strand
                    out.write(f"{chrom}\t{read_id}\t{start}\t{end}\t{read_strand}\t{mapq}\t{type}\t{gene_start}\t{gene_end}\t{gene_strand}\t{sequence}\n")

    # Close BAM file
    bam.close()

def main():
    parser = argparse.ArgumentParser(description='Process BAM and GFF files to output information in a tab-separated file.')
    parser.add_argument('-b', '--bam', required=True, help='Path to the BAM file.')
    parser.add_argument('-g', '--gff', required=True, help='Path to the GFF file.')
    parser.add_argument('-f', '--features', required=False, help='Features to use (CDS).', default='CDS')
    parser.add_argument('-o', '--output', required=True, help='Path to the output file.')

    args = parser.parse_args()

    bam_file = args.bam
    gff_file = args.gff
    features = args.features
    output_file = args.output

    process_bam_and_gff(bam_file, gff_file, features, output_file)

if __name__ == "__main__":
    main()
