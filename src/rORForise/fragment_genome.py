import argparse


def parse_fasta(fasta_file):
    with open(fasta_file, 'r') as file:
        seq = ''
        for line in file:
            if not line.startswith('>'):
                seq += line.strip()
    return seq

def parse_gff_for_cds(gff_file):
    cds_annotations = []
    with open(gff_file, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                parts = line.strip().split('\t')
                if len(parts) > 8 and parts[2] == 'CDS':
                    start, end = int(parts[3]), int(parts[4])
                    cds_annotations.append((start, end))
    return cds_annotations

def split_genome(sequence, cds_annotations, chunk_size=150, overlap=20):
    chunks = []
    for i in range(0, len(sequence), chunk_size - overlap):
        end = min(i + chunk_size, len(sequence))
        chunk_seq = sequence[i:end]
        chunk_cds = [(start, end) for start, end in cds_annotations if start < i + chunk_size and end > i]
        chunks.append((i, end, chunk_seq, chunk_cds))
    return chunks

def write_output(chunks, output_fasta, output_bed):
    with open(output_fasta, "w") as fasta_out, open(output_bed, "w") as bed_out:
        for idx, (chunk_start, chunk_end, chunk_seq, cds_list) in enumerate(chunks):
            chunk_name = f"chunk_{idx}_{chunk_start}_{chunk_end}"
            fasta_out.write(f">{chunk_name}\n{chunk_seq}\n")
            for cds_start, cds_end in cds_list:
                if cds_start < chunk_end and cds_end > chunk_start:
                    adjusted_start = max(cds_start, chunk_start) - chunk_start
                    adjusted_end = min(cds_end, chunk_end) - chunk_start
                    bed_out.write(f"{chunk_name}\t{adjusted_start}\t{adjusted_end}\tCDS\t{cds_start}\t{cds_end}\n")




def main():
    parser = argparse.ArgumentParser(
        description="Fragment a genome FASTA file into chunks and annotate CDS regions from a GFF file."
    )
    parser.add_argument(
        "-fasta", dest="fasta", required=True, help="Path to the input FASTA file."
    )
    parser.add_argument(
        "-gff", dest="gff", required=True, help="Path to the input GFF file."
    )
    parser.add_argument(
        "-chunk-size", dest="chunk-size", type=int, default=150, help="Size of each genome chunk (default: 150)."
    )
    parser.add_argument(
        "-overlap", dest="overlap", type=int, default=20, help="Number of overlapping bases between chunks (default: 20)."
    )
    parser.add_argument(
        "-output-fasta", dest="output-fasta", help="Output FASTA file for genome chunks."
    )
    parser.add_argument(
        "-output-bed", dest="output-bed", help="Output BED file for CDS annotations."
    )

    args = parser.parse_args()

    genome_seq = parse_fasta(args.fasta)
    annotations = parse_gff_for_cds(args.gff)
    chunks = split_genome(genome_seq, annotations, args.chunk_size, args.overlap)
    write_output(chunks, args.output_fasta, args.output_bed)

if __name__ == "__main__":
    main()

# # Example usage
# fasta_file = "Mycoplasma_genitalium_G37/Myco.fa"
# gff_file = "Mycoplasma_genitalium_G37/Myco.gff"
# genome_seq = parse_fasta(fasta_file)
# annotations = parse_gff_for_cds(gff_file)
# chunks = split_genome(genome_seq, annotations, 150, 20)
# write_output(chunks, "output_chunks.fasta", "output_annotations.bed")