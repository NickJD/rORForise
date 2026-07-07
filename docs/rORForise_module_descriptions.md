rORForise — module and tool descriptions

This document lists the modules and CLI tools present under `src/rORForise` and gives a brief description, main functions/classes, and CLI flags for each.

Index
- adapters.py — Shared BED/GFF/FASTA parsing, coordinate conversion, and prediction/CDS overlap utilities.
- bam_gff_intersect.py — Helper to build rORForise-compatible read/CDS mapping files from BAM and GFF.
- check_pred.py — Core prediction checker: given CDS/read/pred positions, categorises correctness.
- compare_tools.py — Wrapper to run evaluations for multiple prediction tools and produce comparisons.
- compute_start_stop_accuracy.py — Script to compute start/stop accuracy per-tool using adapters and check_pred.
- constants.py — Version string.
- count_Read_Per_Gene.py — (brief) helper to compute read counts per gene.
- evaluate.py — Main evaluation logic; computes per-prediction answers, summaries and CSV writers.
- metrics.py — Reusable long-form prediction metrics, edge diagnostics, overlap decomposition, stratified summaries, confidence intervals, pairwise overlap tables, and rankings.
- plots.py — Static Matplotlib figure generation for reports.
- report.py — Interactive HTML report generation from metrics and figures.
- benchmark.py — Multi-tool benchmark orchestration and report generation.
- examples_checking.py — (brief) examples / small checks.
- filter_read_type.py — (brief) filters reads by type.
- fragment_genome.py — (brief) fragments genome for testing.
- gen_test_data.py — Legacy unit tests for prediction logic.
- generate_testing_pred.py — Test-data generator: creates intersect BED, GFF predictions and validation CSV.
- generate_validation_from_preds.py — (brief) create validation rows from predictions.
- gene_coverage.py — (brief) coverage summaries for genes.
- genome_split.py / genome_split_new.py — (brief) split genome into pieces for processing.
- rORForise.py — CLI entry that runs `evaluate.read_preds()` and `evaluate()` with arguments.
- testing.py — (brief) test harness.
- validate_generator.py — Validates generated test data by re-running `check_pred`.
- utils.py — (brief) utility helpers.

Descriptions

## adapters.py
Purpose: Convert GFF-style prediction outputs (either genome coordinates or read-relative coordinates) into a canonical `preds` dictionary suitable for `evaluate.evaluate()`, and provide shared typed records for BED/GFF/FASTA parsing.
- Main API:
  - `BedRecord`, `GffPrediction`
  - `read_fasta_sequences(path)` -> {read_name: sequence}
  - `load_bed_map(intersect_bed_path, read_sequences=None)` -> {read_name: [bed_rows...]}
  - `read_bed_records(intersect_bed_path, read_sequences=None)` -> list[`BedRecord`]
  - `read_gff_predictions(gff_path)` -> list[`GffPrediction`]
  - `gff_genome_to_read_preds(gff_path, intersect_bed_path, coordinate_type="auto", read_sequences=None)` -> {read_name: {pred_id: (pstart, pend, strand)}}
  - `prediction_overlaps_cds(bed_record, pred_start, pred_end)` -> bool
- Notes: Supports the sequence-bearing rORForise TSV/BED layout and BED12/GFF intersects when a read FASTA is supplied. Auto-detects genome versus read-relative GFF coordinates, clamps predictions to read boundaries, suffixes duplicate prediction IDs, and supports gzip.

## bam_gff_intersect.py
Purpose: Create the intersect-style TSV consumed by `evaluate.py` from a BAM alignment and a GFF annotation file.
- CLI flags:
  - `-b/--bam`
  - `-g/--gff`
  - `-f/--features` (default `CDS`, comma-separated)
  - `-o/--output`
- Notes: Uses interval overlap between reads and features, writes the 11-column rORForise mapping format, and requires the optional `pysam` dependency (`pip install rORForise[bam]`).

## check_pred.py
Purpose: Given a CDS mapping, read mapping and a prediction on the read, determine which correctness categories the prediction satisfies (start/stop/frame/direction etc.).
- Public API:
  - `check_pred(cds_open, cds_close, cds_direction, read_open, read_close, read_direction, pred_start, pred_end, pred_direction, read_seq)` -> set of (answer_code, codon)
  - `answers` (dict of name->code), `inverse_answers` (code->name)
- Notes: Core logic used by `evaluate.py` to classify predictions.

## evaluate.py
Purpose: Run full evaluation of predictions against intersect BED (read-to-CDS mapping), tally detailed per-prediction answers, write summary/detailed/codon/track CSVs, and compute simple context accuracies.
- Main API:
  - `read_preds(predictions_gff, output_dir, output_prefix, verbose=False)` -> (total_preds, preds dict)
  - `evaluate(intersect_bed_filename, preds, gc_prob, output_dir, output_prefix, overlap_threshold=60, verbose=False, read_sequences=None)`
- Output files: `{prefix}_run_summary.csv`, `{prefix}_detailed_results.csv`, `{prefix}_codon_summary.csv`, `{prefix}_prediction_lengths.csv`, `{prefix}_track_predictions.csv`, and `{prefix}_codons_*.csv`.
- Notes: Uses shared adapter parsing, computes codon enrichment using provided GC probability, keeps output directories by default, and reports read counts separately from CDS-mapping and prediction/CDS-comparison counts.

## metrics.py
Purpose: Provide a reusable metrics layer for reporting and benchmarking.
- Main API:
  - `evaluate_tool(tool_name, bed_path, preds, overlap_threshold=60, read_sequences=None, tool_type="unspecified")` -> (`prediction_rows`, `summary_row`)
  - `stratified_metric_rows(prediction_rows)` -> rows by read context, strand, boundary class, and codon
  - `boundary_diagnostic_rows(prediction_rows)` -> edge/internal boundary success rows
  - `pairwise_read_overlap(tool_preds)` -> pairwise read-level overlap rows
  - `pairwise_prediction_overlap(tool_preds)` -> exact/reciprocal/same-frame overlap rows
  - `rank_tools(metric_rows)` -> composite ranking rows
- Notes: Emits long-form rows with answer flags, offsets, on-target decomposition, reciprocal overlap, read/CDS context, boundary classes, codons, prediction lengths, and Wilson 95% confidence intervals for key rates.

## plots.py / report.py
Purpose: Generate static visual outputs and an interactive HTML report.
- Figures:
  - `context_accuracy.png`
  - `start_stop_accuracy.png`
  - `frame_direction_breakdown.png`
  - `on_target_rate.png`
  - `prediction_lengths.png`
  - `codon_enrichment_heatmap.png`
  - `read_tracks/*.png`
  - `multi_tool_tracks/*.png`
- Notes: Uses Matplotlib with a writable temp cache directory for CLI/CI environments. The report includes searchable/sortable tables, downloadable CSV links, metric glossary cards, tool scorecards, on-target decomposition, boundary diagnostics, stratified metrics, and read-level examples.

## benchmark.py
Purpose: Benchmark multiple prediction tools from one read/CDS mapping file and produce tables, rankings, pairwise agreement, tool scorecards, stratified diagnostics, legacy per-tool outputs, figures, and `report.html`.
- CLI flags:
  - `-b/--bed`
  - `--manifest`
  - `--reads`
  - `-o/--outdir`
  - `--tool Name=predictions.gff[:auto|read|genome][:tool_type]` (repeatable)
  - `--tool-type Name=tool_type`
  - `--run-tools` for manifest tool commands
  - `--gc_prob`
  - `--coords` default coordinate mode
  - `-l/--overlap_threshold`
  - `--report html|none`
- Output files: `benchmark_summary.csv`, `per_tool_metrics.csv`, `per_prediction_results.csv`, `tool_scorecards.csv`, `boundary_diagnostics.csv`, `stratified_metrics.csv`, `pairwise_read_overlap.csv`, `pairwise_prediction_overlap.csv`, `tool_rankings.csv`, `benchmark_run_manifest.json`, optional `tool_execution_summary.csv`, `evaluation_outputs/*`, `figures/*`, and `report.html`.
- Notes: Tool types include `start_aware`, `stop_to_stop`, `partial_orf`, `classifier`, `genome_caller`, `read_caller`, and `unspecified`. JSON manifests and simple YAML manifests are supported.

## generate_testing_pred.py
Purpose: Create synthetic reads, mapping-to-CDS BED, GFF predictions and a validation CSV describing expected answers for unit testing the evaluation pipeline.
- CLI flags:
  - `-n/--num_reads` (default 200)
  - `-c/--num_cds` (default 50)
  - `-o/--output_dir` (default `test_output`)
  - `--output_prefix` (default `test`)
  - `-gc/--gc_content` (default 0.32)
  - `--seed` (default 42)
  - `--gzip` (write gzipped outputs)
  - `--min_pred_len` (default 30)  # newly added option in this branch
  - `--strand_mode` (choices `mixed`, `plus`, `minus`) controls strand generation
  - `--read_length` (default 150) length of one read in a pair, or the single-read length when `--no-paired` is used
  - `--paired/--no-paired`, `--fragment_length_mean`, `--fragment_length_sd`, `--pair_orientation`
- Outputs: `{prefix}_intersect.bed`, `{prefix}_predictions.gff`, `{prefix}_validation.csv` in output dir.
- Notes: `--strand_mode plus/minus/mixed` controls generated CDS strands and deterministic plus/minus modes avoid random read-strand flips.

## validate_generator.py
Purpose: Given a generated intersect BED and validation CSV, re-run `check_pred` on each validation row and ensure expected answers are present. Exits non-zero on mismatch.
- CLI flags: `-b/--bed`, `-v/--validation`
- Notes: Loads `check_pred` from package or from repo `src` if necessary; reports example failures.

## compare_tools.py
Purpose: Run `evaluate` for multiple tools' prediction GFFs, collect context accuracies and compute pairwise overlap/exact-match statistics between tools.
- CLI flags:
  - `-b/--bed` (intersect bed)
  - `-o/--outdir` (output dir)
  - `--pred` (repeatable) e.g. `--pred Pyrodigal=path/to/pyrodigal.gff`
  - `--gc_prob` (float, default 0.32)
  - `--coords` (`auto`, `read`, `genome`)
  - `-l/--overlap_threshold`
- Output files: `compare_summary.csv`, `pairwise_read_overlap.csv`, plus each tool's evaluate outputs present in `outdir`.
- Notes: Uses `adapters.gff_genome_to_read_preds` to map read-relative or genome-coordinate GFFs before calling `evaluate.evaluate()`.

## compute_start_stop_accuracy.py
Purpose: Script to compute start/stop accuracies for one or more tools by mapping GFF predictions, calling `check_pred`, and counting how many predictions capture and correctly predict starts/stops. Helpful for basic tool benchmarking.
- CLI flags:
  - `-b/--bed` (intersect bed)
  - `--pred` (repeatable) e.g. `--pred Pyrodigal=path/to/pyrodigal.gff`
  - `--out` optional CSV output
- Notes: Loads adapters and check_pred from the local `src/rORForise` files to avoid package import ambiguity; prints per-tool breakdown and optional CSV.

## validate_generator.py
(duplicate above)

## rORForise.py
Purpose: Lightweight CLI wrapper to load predictions via `evaluate.read_preds()` and run `evaluate.evaluate()` with command-line arguments controlling overlap_threshold, gc_prob, output directory, etc.
- CLI flags: `-int_bed/--intersect_bed`, `-p_gff/--predictions_gff`, `-o/--output_dir`, `-gc_prob`, `-prefix/--output_prefix`, `-l/--overlap_threshold`, `--verbose`, `--clear-output`, `--report html|none`, `--report-tool-name`, `--max-tracks`
- Notes: Existing output directories are preserved by default; use `--clear-output` for the old destructive behavior.

## Other modules
- constants.py: defines `rORForise_VERSION`.
- generate_validation_from_preds.py: creates validation rows from read-level GFF predictions and all matching read/CDS mappings.
- fragment_genome.py: fragments a genome into rORForise-compatible read/CDS mapping files; `--compress` writes gzip outputs.
- gen_test_data.py/testing.py/examples_checking.py/gene_coverage.py/gene_split*.py/filter_read_type.py/count_Read_Per_Gene.py/utils.py: small helpers, tests or utilities used in processing and data generation. Inspect individual files for usage if you plan to extend or refactor.
