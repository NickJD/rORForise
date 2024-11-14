import pysam
import argparse
from utils import *

def parse_gff(gff_file):
    """Parse GFF file and return a list of tuples."""
    gene_list = []
    with open(gff_file, 'r') as f:
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
            attributes = cols[8]
            gene_id = attributes.split(";")[0].split("=")[1]
            gene_list.append((seqid, source, type, start, end, score, strand, phase, gene_id))
    return gene_list

def process_bam_and_gff(bam_file, gff_file, features, output_file):
    """Process BAM and GFF files and output information in a tab-separated file."""
    # Parse GFF file
    gff_data = parse_gff(gff_file)

    # Open BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")

    # Open output file
    with open(output_file, 'w') as out:
        # Write header
        out.write("Chromosome\tReadID\tStart\tEnd\tDirection\tMappingQuality\tGeneID\tGeneStart\tGeneEnd\tGeneStrand\tReadSequence\n")

        # Process each read in BAM file
        for read in bam.fetch():
            read_id = read.query_name
            if '502A-748/2' in read_id:
                print("s")
            if read.is_paired:
                if read.is_read1:
                    read_id = f"{read_id}/1"
                elif read.is_read2:
                    read_id = f"{read_id}/2"
            else:
                read_id = f"{read_id}/1"  # Or handle single-end reads differently if needed
            chrom = read.reference_name
            start = read.reference_start +1 # bam file starts are base-0
            end = read.reference_end
            read_strand = '-' if read.is_reverse else '+'
            mapq = read.mapping_quality
            sequence = read.query_sequence
            if read_strand == '-':
                sequence = reverse_complement(sequence)

            # Find corresponding gene in GFF data

            for gene in gff_data:
                (seqid, source, type, gene_start, gene_end, score, strand, phase, gene_id) = gene
                if features != None:
                    if any(type in s for s in features.split(',')):
                        if seqid == chrom and gene_start <= start <= gene_end:
                            gene_strand = strand
                            # Write information to output file
                            out.write(f"{chrom}\t{read_id}\t{start}\t{end}\t{read_strand}\t{mapq}\t{type}\t{gene_start}\t{gene_end}\t{gene_strand}\t{sequence}\n")
                else:
                    if seqid == chrom and gene_start <= start <= gene_end:
                        gene_strand = strand
                        # Write information to output file
                        out.write(
                            f"{chrom}\t{read_id}\t{start}\t{end}\t{read_strand}\t{mapq}\t{type}\t{gene_start}\t{gene_end}\t{gene_strand}\t{sequence}\n")

    # Close BAM file
    bam.close()

def main():
    parser = argparse.ArgumentParser(description='Process BAM and GFF files to output information in a tab-separated file.')
    parser.add_argument('-b', '--bam', required=True, help='Path to the BAM file.')
    parser.add_argument('-g', '--gff', required=True, help='Path to the GFF file.')
    parser.add_argument('-f', '--features', required=False, help='Features to use (CDS).', default=None)
    parser.add_argument('-o', '--output', required=True, help='Path to the output file.')

    args = parser.parse_args()

    bam_file = args.bam
    gff_file = args.gff
    features = args.features
    output_file = args.output

    process_bam_and_gff(bam_file, gff_file, features, output_file)

if __name__ == "__main__":
    main()
