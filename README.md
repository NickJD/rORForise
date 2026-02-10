# rORForise
Read-based gene coverage evaluation

## Test data is available in the Genome_Processing directory.

### Menu - (rORForise -h):  
```bash
usage: rORForise.py [-h] -int_bed INTERSECT_BED -p_gff PREDICTIONS_GFF -o OUTPUT_DIR -gc_prob GC_PROB [-prefix OUTPUT_PREFIX] [-l OVERLAP_THRESHOLD] [--verbose]

rORForise version: v0.0.4 - Process genome evaluation parameters.

options:
  -h, --help            show this help message and exit
  -int_bed, --intersect_bed INTERSECT_BED
                        File path to intersect bedfile with read to CDS mappings.
  -p_gff, --predictions_gff PREDICTIONS_GFF
                        Path to the predictions in GFF format.
  -o, --output_dir OUTPUT_DIR
                        Directory to store output files.
  -gc_prob GC_PROB      GC probability of the genome being processed (e.g., 0.3169 for Mycoplasma genitalium).
  -prefix, --output_prefix OUTPUT_PREFIX
                        Prefix to add to output files. Default is orf_evaluation.
  -l, --overlap_threshold OVERLAP_THRESHOLD
                        Minimum number of bases of overlap that are required for the read to overlap the CDS by before a prediction is inspected.
  --verbose             Verbose output

```
## Example Run:
```bash
rORForise-int_bed ~/Genome_Processing/Escherichia_coli_k_12/Processing/Escherichia_coli_k_12_Reads_Intersect.tsv.gz -p_gff ~/Genome_Processing/Escherichia_coli_k_12/FragGeneScan/FragGeneScan_ART_errFree_Combined.gff.gz -o ~/rORForise_FGS_Ecoli_testing  -gc 0.39```
```