rORForise — module and tool descriptions

This document lists the modules and CLI tools present under `src/rORForise` and gives a brief description, main functions/classes, and CLI flags for each.

Index
- adapters.py — Adapter utilities to map GFF predictions to read-relative coordinates used by evaluation.
- bam_gff_intersect.py — (brief) helper to intersect BAM/GFF (not yet summarized)
- check_pred.py — Core prediction checker: given CDS/read/pred positions, categorises correctness.
- compare_tools.py — Wrapper to run evaluations for multiple prediction tools and produce comparisons.
- compute_start_stop_accuracy.py — Script to compute start/stop accuracy per-tool using adapters and check_pred.
- constants.py — Version string.
- count_Read_Per_Gene.py — (brief) helper to compute read counts per gene.
- evaluate.py — Main evaluation logic; computes per-prediction answers, summaries and CSV writers.
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
Purpose: Convert GFF-style prediction outputs (either genome coordinates or read-relative coordinates) into a canonical `preds` dictionary suitable for `evaluate.evaluate()`.
- Main API:
  - `load_bed_map(intersect_bed_path)` -> {read_name: [bed_rows...]}
  - `gff_genome_to_read_preds(gff_path, intersect_bed_path)` -> {read_name: {pred_id: (pstart, pend, strand)}}
- Notes: Auto-detects whether the GFF coordinates are genome coordinates or already read-relative; clamps predictions to read boundaries; supports gzipped files.

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
  - `evaluate(intersect_bed_filename, preds, gc_prob, output_dir, output_prefix, overlap_threshold=60, verbose=False)`
- Output files: `{prefix}_run_summary.csv`, `{prefix}_detailed_results.csv`, `{prefix}_codon_summary.csv`, `{prefix}_prediction_lengths.csv`, `{prefix}_track_predictions.csv`, and `{prefix}_codons_*.csv`.
- Notes: Uses `check_pred` by loading `check_pred.py` relative to the module; computes codon enrichment using provided GC probability.

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
  - `--read_length` (default 300) read length in bases used for generated reads (default 300 = paired-end 150bp x2)
- Outputs: `{prefix}_intersect.bed`, `{prefix}_predictions.gff`, `{prefix}_validation.csv` in output dir.

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
- Output files: `compare_summary.csv`, `pairwise_read_overlap.csv`, plus each tool's evaluate outputs present in `outdir`.
- Notes: Uses `adapters.gff_genome_to_read_preds` to map GFFs to read-relative preds then calls `evaluate.evaluate()` for each tool. Also computes counts of predictions that overlap start/stop/middle for that tool's predictions.

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
- CLI flags: `-int_bed/--intersect_bed`, `-p_gff/--predictions_gff`, `-o/--output_dir`, `-gc_prob`, `-prefix/--output_prefix`, `-l/--overlap_threshold`, `--verbose`

## Other modules
- constants.py: defines `rORForise_VERSION`.
- gen_test_data.py/testing.py/examples_checking.py/fragment_genome.py/gene_coverage.py/generate_validation_from_preds.py/gene_split*.py/filter_read_type.py/count_Read_Per_Gene.py/utils.py: small helpers, tests or utilities used in processing and data generation. Inspect individual files for usage if you plan to extend or refactor.
