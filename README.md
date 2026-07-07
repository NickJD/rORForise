# rORForise
Read-based gene coverage evaluation

## Test data is available in the Genome_Processing directory.

## Install
```bash
pip install rORForise
```

### Menu - (rORForise -h):  
```bash
usage: rORForise.py [-h] -int_bed INTERSECT_BED -p_gff PREDICTIONS_GFF -o OUTPUT_DIR -gc_prob GC_PROB [-prefix OUTPUT_PREFIX] [-l OVERLAP_THRESHOLD] [--verbose] [--clear-output] [--report {none,html}]

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
  --clear-output        Delete the output directory before writing results.
  --report {none,html}  Generate report.html and static figures.

```
## Example Run:
```bash
rORForise -int_bed ~/Genome_Processing/Escherichia_coli_k_12/Processing/Escherichia_coli_k_12_Reads_Intersect.tsv.gz -p_gff ~/Genome_Processing/Escherichia_coli_k_12/FragGeneScan/FragGeneScan_ART_errFree_Combined.gff.gz -o ~/rORForise_FGS_Ecoli_testing -gc_prob 0.39 --report html --report-tool-name FragGeneScan
```

## Multiple-Tool Benchmarking
```bash
rORForise-benchmark \
  -b reads_intersect.tsv.gz \
  --tool FragGeneScan=fgs_predictions.gff.gz:read:start_aware \
  --tool Pyrodigal=pyrodigal_predictions.gff.gz:genome:start_aware \
  --tool Naive=naive_predictions.gff.gz:read:stop_to_stop \
  -o benchmark_out \
  --gc_prob 0.39 \
  --report html
```

For BED12/GFF intersects without embedded read sequence, pass the reads separately:
```bash
rORForise-benchmark \
  -b reads_intersect.bed.gz \
  --reads simulated_reads.fasta.gz \
  --tool FragGeneScan=fgs_predictions.gff.gz:read:start_aware \
  -o benchmark_out \
  --gc_prob 0.39
```

Benchmark manifests can be JSON or simple YAML:
```yaml
intersect: Processing/Reads_Intersect_CDS.tsv.gz
reads: Processing/ART_Simulated_Reads/reads.fasta.gz
tools:
  - name: FragGeneScan
    predictions: FragGeneScan/predictions.gff.gz
    coords: read
    type: start_aware
  - name: NaiveStORF
    predictions: Naive/predictions.gff.gz
    coords: read
    type: stop_to_stop
```

Run a manifest with:
```bash
rORForise-benchmark --manifest benchmark.yaml -o benchmark_out --gc_prob 0.39
```

Benchmark outputs include `benchmark_summary.csv`, `per_tool_metrics.csv`, `per_prediction_results.csv`, `tool_scorecards.csv`, `boundary_diagnostics.csv`, `stratified_metrics.csv`, `pairwise_read_overlap.csv`, `pairwise_prediction_overlap.csv`, `tool_rankings.csv`, `benchmark_run_manifest.json`, legacy per-tool evaluation outputs, `figures/*.png`, and an interactive `report.html`.

## Development Checks
```bash
PYTHONPATH=src python -m unittest discover -s tests -v
PYTHONPATH=src python -m compileall -q src tests
```
