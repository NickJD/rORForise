import argparse
import json
import os
import subprocess
import time
from datetime import datetime
from pathlib import Path

try:
    from . import adapters, evaluate, metrics, report
except Exception:
    import adapters
    import evaluate
    import metrics
    import report


COORD_CHOICES = {"auto", "read", "genome"}
TOOL_TYPE_CHOICES = {
    "unspecified",
    "start_aware",
    "stop_to_stop",
    "partial_orf",
    "classifier",
    "genome_caller",
    "read_caller",
}


def infer_tool_type(name):
    lowered = name.lower()
    if "naivestorf" in lowered or "naive-storf" in lowered:
        return "stop_to_stop"
    if "fraggenescan" in lowered or "pyrodigal" in lowered or "prodigal" in lowered:
        return "start_aware"
    return "unspecified"


def parse_tool_spec(spec, default_coords="auto"):
    if "=" not in spec:
        raise argparse.ArgumentTypeError("Tool specs must look like Name=path/to/predictions.gff[:coords][:tool_type]")
    name, remainder = spec.split("=", 1)
    coords = default_coords
    tool_type = None
    path = remainder
    parts = remainder.rsplit(":", 2)
    if len(parts) == 3 and parts[-1] in TOOL_TYPE_CHOICES and parts[-2] in COORD_CHOICES:
        path, coords, tool_type = parts
    else:
        possible_path, possible_suffix = remainder.rsplit(":", 1) if ":" in remainder else (remainder, None)
        if possible_suffix in COORD_CHOICES:
            path = possible_path
            coords = possible_suffix
        elif possible_suffix in TOOL_TYPE_CHOICES:
            path = possible_path
            tool_type = possible_suffix
    if not name:
        raise argparse.ArgumentTypeError("Tool name cannot be empty")
    return name, path, coords, tool_type or infer_tool_type(name)


def parse_name_value(spec, label):
    if "=" not in spec:
        raise argparse.ArgumentTypeError(f"{label} specs must look like Name=value")
    name, value = spec.split("=", 1)
    if not name or not value:
        raise argparse.ArgumentTypeError(f"{label} specs must look like Name=value")
    return name, value


def _simple_yaml_manifest(text):
    data = {"tools": []}
    current_tool = None
    for raw_line in text.splitlines():
        line = raw_line.split("#", 1)[0].rstrip()
        if not line.strip():
            continue
        stripped = line.strip()
        if stripped.startswith("- "):
            current_tool = {}
            data.setdefault("tools", []).append(current_tool)
            stripped = stripped[2:]
            if ":" in stripped:
                key, value = stripped.split(":", 1)
                current_tool[key.strip()] = value.strip().strip("'\"")
            continue
        if ":" not in stripped:
            continue
        key, value = stripped.split(":", 1)
        key = key.strip()
        value = value.strip().strip("'\"")
        if raw_line.startswith((" ", "\t")) and current_tool is not None:
            current_tool[key] = value
        else:
            data[key] = value
            current_tool = None
    return data


def load_manifest(path):
    path = Path(path)
    text = path.read_text(encoding="utf-8")
    if path.suffix.lower() == ".json":
        data = json.loads(text)
    else:
        data = _simple_yaml_manifest(text)
    base = path.parent
    for key in ("intersect", "reads"):
        if data.get(key):
            data[key] = str((base / data[key]).resolve()) if not os.path.isabs(data[key]) else data[key]
    for tool in data.get("tools", []):
        for key in ("predictions", "prediction_gff", "command_cwd"):
            if tool.get(key) and not os.path.isabs(tool[key]):
                tool[key] = str((base / tool[key]).resolve())
    return data


def tool_type_rows(tool_configs):
    descriptions = {
        "start_aware": "Designed to infer biological start/stop boundaries.",
        "stop_to_stop": "Expected to find uninterrupted stop/edge-delimited ORF segments; start accuracy is mostly boundary coincidence.",
        "partial_orf": "Expected to recover partial coding sequence on reads.",
        "classifier": "Expected to classify coding potential more than exact ORF boundaries.",
        "genome_caller": "Predictions are naturally in genome coordinates.",
        "read_caller": "Predictions are naturally in read coordinates.",
        "unspecified": "No tool behavior metadata supplied.",
    }
    return [
        {
            "tool": config["name"],
            "tool_type": config["tool_type"],
            "interpretation": descriptions.get(config["tool_type"], descriptions["unspecified"]),
            "prediction_gff": config["pred_path"],
            "coordinate_type": config["coords"],
        }
        for config in tool_configs
    ]


def build_tool_configs(tool_specs, default_coords="auto", tool_type_specs=None, manifest_tools=None):
    explicit_types = {}
    for spec in tool_type_specs or []:
        name, value = parse_name_value(spec, "--tool-type")
        explicit_types[name] = value

    configs = []
    for raw_spec in tool_specs or []:
        name, pred_path, coords, tool_type = parse_tool_spec(raw_spec, default_coords=default_coords)
        configs.append({"name": name, "pred_path": pred_path, "coords": coords, "tool_type": explicit_types.get(name, tool_type)})

    for tool in manifest_tools or []:
        name = tool.get("name")
        pred_path = tool.get("predictions") or tool.get("prediction_gff")
        if not name or not pred_path:
            continue
        configs.append(
            {
                "name": name,
                "pred_path": pred_path,
                "coords": tool.get("coords") or default_coords,
                "tool_type": explicit_types.get(name, tool.get("type") or infer_tool_type(name)),
                "command": tool.get("command"),
                "command_cwd": tool.get("command_cwd"),
            }
        )
    return configs


def run_tool_commands(tool_configs, output_dir):
    rows = []
    for config in tool_configs:
        command = config.get("command")
        if not command:
            continue
        start = time.monotonic()
        completed = subprocess.run(
            command,
            shell=True,
            cwd=config.get("command_cwd") or None,
            text=True,
            capture_output=True,
        )
        elapsed = time.monotonic() - start
        rows.append(
            {
                "tool": config["name"],
                "command": command,
                "returncode": completed.returncode,
                "runtime_seconds": elapsed,
                "stdout_tail": completed.stdout[-1000:],
                "stderr_tail": completed.stderr[-1000:],
            }
        )
        if completed.returncode != 0:
            metrics.write_csv_dicts(Path(output_dir) / "tool_execution_summary.csv", rows)
            raise RuntimeError(f"Tool command failed for {config['name']} with exit code {completed.returncode}")
    return rows


def selected_summary(row):
    return {
        "tool": row["tool"],
        "tool_type": row.get("tool_type"),
        "start_accuracy_pct": row.get("start_accuracy_pct"),
        "stop_accuracy_pct": row.get("stop_accuracy_pct"),
        "middle_accuracy_pct": row.get("middle_accuracy_pct"),
        "frame_accuracy_pct": row.get("frame_accuracy_pct"),
        "direction_accuracy_pct": row.get("direction_accuracy_pct"),
        "on_target_rate_pct": row.get("on_target_rate_pct"),
        "reciprocal_overlap_80_pct": row.get("reciprocal_overlap_80_pct"),
        "avg_prediction_overlap_pct": row.get("avg_prediction_overlap_pct"),
        "avg_cds_covered_pct": row.get("avg_cds_covered_pct"),
        "total_predictions": row.get("total_predictions"),
        "inspected_prediction_cds_comparisons": row.get("inspected_prediction_cds_comparisons"),
    }


def run_benchmark(
    bed_path,
    tool_specs,
    output_dir,
    gc_prob=0.32,
    overlap_threshold=60,
    default_coords="auto",
    report_format="html",
    max_tracks=12,
    reads_path=None,
    tool_type_specs=None,
    manifest_path=None,
    run_tools=False,
):
    manifest = load_manifest(manifest_path) if manifest_path else {}
    bed_path = manifest.get("intersect") or bed_path
    reads_path = manifest.get("reads") or reads_path
    if not bed_path:
        raise ValueError("An intersect file is required via --bed or manifest intersect.")
    read_sequences = adapters.read_fasta_sequences(reads_path) if reads_path else None
    tool_configs = build_tool_configs(
        tool_specs,
        default_coords=default_coords,
        tool_type_specs=tool_type_specs,
        manifest_tools=manifest.get("tools"),
    )
    if not tool_configs:
        raise ValueError("At least one --tool or manifest tool entry is required.")

    records = adapters.read_bed_records(bed_path, read_sequences=read_sequences)
    if not records:
        raise ValueError(
            "No readable intersect records were found. "
            "Use the sequence-bearing rORForise TSV/BED layout with read_name, read_start, "
            "read_end, read_strand, feature_type, CDS coordinates, CDS strand, and read_sequence, "
            "or provide --reads FASTA for BED12/GFF intersects."
        )
    if not any(record.feature_type == "CDS" for record in records):
        raise ValueError("No CDS records were found in the intersect input.")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    eval_dir = output_dir / "evaluation_outputs"
    eval_dir.mkdir(exist_ok=True)

    execution_rows = run_tool_commands(tool_configs, output_dir) if run_tools else []

    all_prediction_rows = []
    metric_rows = []
    all_codon_rows = []
    tool_preds = {}

    for config in tool_configs:
        tool_name = config["name"]
        pred_path = config["pred_path"]
        coords = config["coords"]
        tool_type = config["tool_type"]
        prefix = metrics.safe_name(tool_name)
        preds = adapters.gff_genome_to_read_preds(pred_path, bed_path, coordinate_type=coords, read_sequences=read_sequences)
        tool_preds[tool_name] = preds

        evaluate.evaluate(
            bed_path,
            preds,
            gc_prob,
            eval_dir,
            output_prefix=prefix,
            overlap_threshold=overlap_threshold,
            verbose=False,
            read_sequences=read_sequences,
        )

        prediction_rows, metric_row = metrics.evaluate_tool(
            tool_name,
            bed_path,
            preds,
            overlap_threshold=overlap_threshold,
            read_sequences=read_sequences,
            tool_type=tool_type,
        )
        metric_row["prediction_gff"] = pred_path
        metric_row["coordinate_type"] = coords
        all_prediction_rows.extend(prediction_rows)
        metric_rows.append(metric_row)
        all_codon_rows.extend(metrics.load_codon_summary(eval_dir / f"{prefix}_codon_summary.csv", tool_name))

    read_overlap_rows = metrics.pairwise_read_overlap(tool_preds)
    prediction_overlap_rows = metrics.pairwise_prediction_overlap(tool_preds)
    stratified_rows = metrics.stratified_metric_rows(all_prediction_rows)
    boundary_rows = metrics.boundary_diagnostic_rows(all_prediction_rows)
    ranking_rows = metrics.rank_tools(metric_rows)
    summary_rows = [selected_summary(row) for row in metric_rows]
    scorecard_rows = tool_type_rows(tool_configs)

    metrics.write_csv_dicts(output_dir / "per_prediction_results.csv", all_prediction_rows)
    metrics.write_csv_dicts(output_dir / "per_tool_metrics.csv", metric_rows)
    metrics.write_csv_dicts(output_dir / "benchmark_summary.csv", summary_rows)
    metrics.write_csv_dicts(output_dir / "tool_scorecards.csv", scorecard_rows)
    metrics.write_csv_dicts(output_dir / "stratified_metrics.csv", stratified_rows)
    metrics.write_csv_dicts(output_dir / "boundary_diagnostics.csv", boundary_rows)
    metrics.write_csv_dicts(output_dir / "pairwise_read_overlap.csv", read_overlap_rows)
    metrics.write_csv_dicts(output_dir / "pairwise_prediction_overlap.csv", prediction_overlap_rows)
    metrics.write_csv_dicts(output_dir / "tool_rankings.csv", ranking_rows)
    if execution_rows:
        metrics.write_csv_dicts(output_dir / "tool_execution_summary.csv", execution_rows)

    run_manifest = {
        "generated_at": datetime.now().isoformat(),
        "intersect": str(bed_path),
        "reads": str(reads_path) if reads_path else "",
        "gc_prob": gc_prob,
        "overlap_threshold": overlap_threshold,
        "tools": tool_configs,
        "outputs": [
            "per_prediction_results.csv",
            "per_tool_metrics.csv",
            "benchmark_summary.csv",
            "tool_scorecards.csv",
            "stratified_metrics.csv",
            "boundary_diagnostics.csv",
            "pairwise_read_overlap.csv",
            "pairwise_prediction_overlap.csv",
            "tool_rankings.csv",
            "report.html" if report_format == "html" else "",
        ],
    }
    (output_dir / "benchmark_run_manifest.json").write_text(json.dumps(run_manifest, indent=2), encoding="utf-8")

    report_path = None
    if report_format == "html":
        report_path = report.build_report(
            output_dir,
            "rORForise Benchmark Report",
            metric_rows,
            all_prediction_rows,
            codon_rows=all_codon_rows,
            ranking_rows=ranking_rows,
            pairwise_rows=prediction_overlap_rows,
            stratified_rows=stratified_rows,
            boundary_rows=boundary_rows,
            scorecard_rows=scorecard_rows,
            max_tracks=max_tracks,
        )

    return {
        "output_dir": output_dir,
        "report": report_path,
        "metric_rows": metric_rows,
        "prediction_rows": all_prediction_rows,
        "ranking_rows": ranking_rows,
        "stratified_rows": stratified_rows,
        "boundary_rows": boundary_rows,
    }


def parse_args():
    parser = argparse.ArgumentParser(description="Benchmark multiple read-level ORF prediction tools.")
    parser.add_argument("-b", "--bed", help="Intersect BED/TSV file with read-to-CDS mappings.")
    parser.add_argument("--manifest", help="Benchmark manifest JSON/YAML containing intersect, reads, and tools.")
    parser.add_argument("--reads", help="Optional read FASTA/FASTA.GZ used when the intersect file lacks read_sequence.")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory for benchmark tables and report.")
    parser.add_argument(
        "--tool",
        action="append",
        help="Tool spec Name=predictions.gff[:coords][:tool_type], where coords is auto/read/genome.",
    )
    parser.add_argument("--tool-type", action="append", help="Tool behavior spec Name=tool_type.")
    parser.add_argument("--run-tools", action="store_true", help="Run manifest tool commands before benchmarking.")
    parser.add_argument("--gc_prob", type=float, default=0.32, help="GC probability for codon expected counts.")
    parser.add_argument("--coords", choices=sorted(COORD_CHOICES), default="auto", help="Default coordinate mode for tools.")
    parser.add_argument("-l", "--overlap_threshold", type=int, default=60, help="Minimum read/CDS overlap to inspect.")
    parser.add_argument("--report", choices=["html", "none"], default="html", help="Report format to generate.")
    parser.add_argument("--max-tracks", type=int, default=12, help="Maximum read-level track examples in the HTML report.")
    return parser.parse_args()


def main():
    args = parse_args()
    result = run_benchmark(
        args.bed,
        args.tool,
        args.outdir,
        gc_prob=args.gc_prob,
        overlap_threshold=args.overlap_threshold,
        default_coords=args.coords,
        report_format=args.report,
        max_tracks=args.max_tracks,
        reads_path=args.reads,
        tool_type_specs=args.tool_type,
        manifest_path=args.manifest,
        run_tools=args.run_tools,
    )
    print(f"Benchmark complete: {result['output_dir']}")
    if result["report"]:
        print(f"HTML report: {result['report']}")


if __name__ == "__main__":
    main()
