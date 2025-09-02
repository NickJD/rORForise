

try: # Try to import from the package if available
    from .constants import *
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from constants import *

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


import argparse


def main():
    parser = argparse.ArgumentParser(description="rORForise version: " + rORForise_VERSION + " - genome-split.")

    parser.add_argument('-f', '--fasta_file', type=str, required=True,
                        help="Path to the input FASTA file containing the genome sequence.")

    parser.add_argument('-g', '--gff_file', type=str, required=True,
                        help="Path to the input GFF file containing genome annotations.")

    parser.add_argument('-c', '--chunk_size', type=int, default=150,
                        help="Size of genome chunks to split into (default: 150).")

    parser.add_argument('-o', '--overlap', type=int, default=20,
                        help="Overlap size between genome chunks (default: 20).")

    parser.add_argument('-output_fasta', type=str, default="output_chunks.fasta",
                        help="Path to the output FASTA file for genome chunks (default: output_chunks.fasta).")

    parser.add_argument('-output_bed', type=str, default="output_annotations.bed",
                        help="Path to the output BED file for annotations (default: output_annotations.bed).")


    options = parser.parse_args()


    genome_seq = parse_fasta(options.fasta_file)
    annotations = parse_gff_for_cds(options.gff_file)
    chunks = split_genome(genome_seq, annotations, options.chunk_size, options.overlap)
    write_output(chunks, options.output_fasta, options.output_bed)





if __name__ == "__main__":
    main()
    print("Complete")