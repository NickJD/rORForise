from html import escape
from pathlib import Path

try:
    from . import metrics
    from . import plots
except Exception:
    import metrics
    import plots


def _relative(path, base):
    return Path(path).resolve().relative_to(Path(base).resolve()).as_posix()


def _fmt(key, value):
    if value is None:
        return "N/A"
    if value == "":
        return ""
    if key.endswith("_pct") or key == "composite_score":
        try:
            return f"{float(value):.2f}"
        except (TypeError, ValueError):
            return str(value)
    if isinstance(value, float):
        return f"{value:.2f}"
    return str(value)


def _table(rows, columns, limit=None):
    chosen = list(rows)
    if limit is not None:
        chosen = chosen[:limit]
    head = "".join(f"<th>{escape(label)}</th>" for _, label in columns)
    body_rows = []
    for row in chosen:
        cells = "".join(f"<td>{escape(_fmt(key, row.get(key, '')))}</td>" for key, _ in columns)
        body_rows.append(f"<tr>{cells}</tr>")
    if not body_rows:
        body_rows.append(f"<tr><td colspan=\"{len(columns)}\">No data</td></tr>")
    return f"<table class=\"data-table\"><thead><tr>{head}</tr></thead><tbody>{''.join(body_rows)}</tbody></table>"


def _figure(title, path, output_dir):
    if not path:
        return ""
    rel = _relative(path, output_dir)
    return f"<figure><img src=\"{escape(rel)}\" alt=\"{escape(title)}\"><figcaption>{escape(title)}</figcaption></figure>"


def _download_links(output_dir):
    files = [
        "benchmark_summary.csv",
        "per_tool_metrics.csv",
        "tool_rankings.csv",
        "tool_scorecards.csv",
        "boundary_diagnostics.csv",
        "stratified_metrics.csv",
        "pairwise_prediction_overlap.csv",
        "per_prediction_results.csv",
        "benchmark_run_manifest.json",
    ]
    links = []
    for name in files:
        path = output_dir / name
        if path.exists():
            links.append(f"<a href=\"{escape(name)}\">{escape(name)}</a>")
    return "<div class=\"downloads\">" + "".join(links) + "</div>" if links else ""


def build_report(
    output_dir,
    title,
    metric_rows,
    prediction_rows,
    codon_rows=None,
    ranking_rows=None,
    pairwise_rows=None,
    stratified_rows=None,
    boundary_rows=None,
    scorecard_rows=None,
    max_tracks=10,
):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir = output_dir / "figures"
    tracks_dir = figures_dir / "read_tracks"
    comparisons_dir = figures_dir / "multi_tool_tracks"
    figures_dir.mkdir(parents=True, exist_ok=True)

    codon_rows = codon_rows or []
    ranking_rows = ranking_rows or []
    pairwise_rows = pairwise_rows or []
    stratified_rows = stratified_rows or []
    boundary_rows = boundary_rows or []
    scorecard_rows = scorecard_rows or []

    figure_paths = [
        ("Context Accuracy", plots.plot_context_accuracy(metric_rows, figures_dir / "context_accuracy.png")),
        ("Start/Stop Boundary Accuracy", plots.plot_start_stop_tolerance(metric_rows, figures_dir / "start_stop_accuracy.png")),
        ("Frame and Direction Accuracy", plots.plot_frame_direction(metric_rows, figures_dir / "frame_direction_breakdown.png")),
        ("On-Target Rate", plots.plot_on_target(metric_rows, figures_dir / "on_target_rate.png")),
        ("Prediction Lengths", plots.plot_prediction_lengths(prediction_rows, figures_dir / "prediction_lengths.png")),
        ("Codon Enrichment Heatmap", plots.plot_codon_enrichment(codon_rows, figures_dir / "codon_enrichment_heatmap.png")),
    ]
    track_paths = plots.plot_read_tracks(prediction_rows, tracks_dir, max_examples=max_tracks)
    comparison_paths = plots.plot_multitool_read_tracks(prediction_rows, comparisons_dir, max_reads=max(4, max_tracks // 3))

    metric_columns = [
        ("tool", "Tool"),
        ("tool_type", "Type"),
        ("start_accuracy_pct", "Start %"),
        ("stop_accuracy_pct", "Stop %"),
        ("middle_accuracy_pct", "Middle %"),
        ("frame_accuracy_pct", "Frame %"),
        ("direction_accuracy_pct", "Direction %"),
        ("on_target_rate_pct", "On-target %"),
        ("reciprocal_overlap_80_pct", "Reciprocal 80% %"),
        ("total_predictions", "Predictions"),
    ]
    ranking_columns = [
        ("rank", "Rank"),
        ("tool", "Tool"),
        ("tool_type", "Type"),
        ("composite_score", "Composite"),
        ("start_accuracy_pct", "Start %"),
        ("stop_accuracy_pct", "Stop %"),
        ("on_target_rate_pct", "On-target %"),
    ]
    scorecard_columns = [
        ("tool", "Tool"),
        ("tool_type", "Type"),
        ("coordinate_type", "Coords"),
        ("interpretation", "Interpretation"),
    ]
    overlap_columns = [
        ("tool", "Tool"),
        ("on_target_rate_pct", "Any CDS overlap %"),
        ("reciprocal_overlap_80_pct", "Reciprocal 80% %"),
        ("prediction_fully_inside_cds_pct", "Prediction inside CDS %"),
        ("cds_fully_inside_prediction_pct", "CDS inside prediction %"),
        ("avg_prediction_overlap_pct", "Avg prediction overlap %"),
        ("avg_cds_covered_pct", "Avg CDS covered %"),
    ]
    boundary_columns = [
        ("tool", "Tool"),
        ("boundary", "Boundary"),
        ("boundary_class", "Class"),
        ("eligible_comparisons", "Eligible"),
        ("correct", "Correct"),
        ("correct_pct", "Correct %"),
    ]
    stratified_columns = [
        ("tool", "Tool"),
        ("stratification", "Stratification"),
        ("value", "Value"),
        ("comparisons", "Comparisons"),
        ("start_accuracy_pct", "Start %"),
        ("stop_accuracy_pct", "Stop %"),
        ("frame_accuracy_pct", "Frame %"),
        ("direction_accuracy_pct", "Direction %"),
        ("on_target_rate_pct", "On-target %"),
    ]
    prediction_columns = [
        ("tool", "Tool"),
        ("read_name", "Read"),
        ("prediction_id", "Prediction"),
        ("pred_start", "Start"),
        ("pred_end", "End"),
        ("pred_strand", "Strand"),
        ("answers", "Answers"),
    ]
    pairwise_columns = [
        ("tool_a", "Tool A"),
        ("tool_b", "Tool B"),
        ("shared_reads", "Shared Reads"),
        ("exact_pred_matches", "Exact Matches"),
        ("reciprocal_overlap_80", "80% Reciprocal Overlap"),
        ("same_frame_overlaps", "Same Frame"),
    ]

    figures_html = "".join(_figure(label, path, output_dir) for label, path in figure_paths)
    tracks_html = "".join(_figure(path.name, path, output_dir) for path in track_paths)
    comparisons_html = "".join(_figure(path.name, path, output_dir) for path in comparison_paths)
    downloads_html = _download_links(output_dir)

    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{escape(title)}</title>
  <style>
    body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; margin: 0; color: #20242a; background: #f7f8fa; }}
    header {{ background: #1f2933; color: white; padding: 28px 36px; }}
    main {{ max-width: 1180px; margin: 0 auto; padding: 28px 24px 48px; }}
    h1 {{ margin: 0; font-size: 30px; }}
    h2 {{ margin-top: 34px; border-bottom: 1px solid #d8dde6; padding-bottom: 8px; }}
    h3 {{ margin-top: 24px; }}
    .grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(360px, 1fr)); gap: 18px; }}
    figure {{ background: white; border: 1px solid #d8dde6; border-radius: 6px; padding: 12px; margin: 0; }}
    figure img {{ display: block; max-width: 100%; height: auto; }}
    figcaption {{ font-size: 13px; color: #52606d; margin-top: 8px; }}
    table {{ border-collapse: collapse; width: 100%; background: white; border: 1px solid #d8dde6; }}
    .table-wrap {{ overflow-x: auto; margin-bottom: 18px; }}
    th, td {{ padding: 8px 10px; border-bottom: 1px solid #e4e7ec; text-align: left; font-size: 13px; }}
    th {{ background: #eef2f7; font-weight: 650; cursor: pointer; user-select: none; }}
    code {{ background: #eef2f7; padding: 2px 5px; border-radius: 4px; }}
    .note {{ color: #52606d; }}
    .toolbar {{ display: flex; gap: 12px; align-items: center; margin: 14px 0; }}
    .toolbar input {{ max-width: 360px; width: 100%; padding: 8px 10px; border: 1px solid #cbd2d9; border-radius: 6px; }}
    .downloads {{ display: flex; flex-wrap: wrap; gap: 8px; margin: 12px 0 20px; }}
    .downloads a {{ background: white; border: 1px solid #cbd2d9; border-radius: 6px; color: #1f4f82; padding: 7px 9px; text-decoration: none; font-size: 13px; }}
    .glossary {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(260px, 1fr)); gap: 10px; }}
    .glossary div {{ background: white; border: 1px solid #d8dde6; border-radius: 6px; padding: 10px; font-size: 13px; }}
  </style>
</head>
<body>
<header>
  <h1>{escape(title)}</h1>
  <p class="note">Generated by rORForise.</p>
</header>
<main>
  <div class="toolbar">
    <input id="table-filter" type="search" placeholder="Filter tables">
  </div>

  <h2>Downloads</h2>
  {downloads_html}

  <h2>Glossary</h2>
  <div class="glossary">
    <div><strong>On-target</strong><br>Prediction overlaps the annotated CDS interval by at least 1 nt.</div>
    <div><strong>Reciprocal 80%</strong><br>Prediction and CDS each cover at least 80% of the other interval.</div>
    <div><strong>Start/stop accuracy</strong><br>Exact biological boundary match for eligible read/CDS/prediction comparisons.</div>
    <div><strong>Boundary class</strong><br>Whether the assessed prediction boundary sits on a read edge, near an edge, or internally.</div>
  </div>

  <h2>Summary</h2>
  <div class="table-wrap">{_table(metric_rows, metric_columns)}</div>

  <h2>Tool Scorecards</h2>
  <div class="table-wrap">{_table(scorecard_rows, scorecard_columns) if scorecard_rows else '<p class="note">No tool metadata supplied.</p>'}</div>

  <h2>On-Target Decomposition</h2>
  <div class="table-wrap">{_table(metric_rows, overlap_columns)}</div>

  <h2>Figures</h2>
  <div class="grid">{figures_html}</div>

  <h2>Tool Rankings</h2>
  <div class="table-wrap">{_table(ranking_rows, ranking_columns) if ranking_rows else '<p class="note">No ranking table for this report.</p>'}</div>

  <h2>Boundary Diagnostics</h2>
  <div class="table-wrap">{_table(boundary_rows, boundary_columns, limit=200) if boundary_rows else '<p class="note">No boundary diagnostics available.</p>'}</div>

  <h2>Stratified Metrics</h2>
  <div class="table-wrap">{_table(stratified_rows, stratified_columns, limit=200) if stratified_rows else '<p class="note">No stratified metrics available.</p>'}</div>

  <h2>Pairwise Prediction Overlap</h2>
  <div class="table-wrap">{_table(pairwise_rows, pairwise_columns) if pairwise_rows else '<p class="note">Pairwise overlap requires at least two tools.</p>'}</div>

  <h2>Read-Level Tracks</h2>
  <div class="grid">{tracks_html or '<p class="note">No track examples available.</p>'}</div>

  <h2>Multi-Tool Read Comparisons</h2>
  <div class="grid">{comparisons_html or '<p class="note">No multi-tool read comparisons available.</p>'}</div>

  <h2>Example Prediction Rows</h2>
  <div class="table-wrap">{_table(prediction_rows, prediction_columns, limit=50)}</div>
</main>
<script>
const filterInput = document.getElementById('table-filter');
filterInput.addEventListener('input', () => {{
  const q = filterInput.value.toLowerCase();
  document.querySelectorAll('tbody tr').forEach(row => {{
    row.style.display = row.textContent.toLowerCase().includes(q) ? '' : 'none';
  }});
}});

document.querySelectorAll('th').forEach((th, columnIndex) => {{
  th.addEventListener('click', () => {{
    const table = th.closest('table');
    const tbody = table.querySelector('tbody');
    const rows = Array.from(tbody.querySelectorAll('tr'));
    const ascending = th.dataset.sort !== 'asc';
    rows.sort((a, b) => {{
      const av = a.children[columnIndex]?.textContent.trim() || '';
      const bv = b.children[columnIndex]?.textContent.trim() || '';
      const an = Number(av);
      const bn = Number(bv);
      const cmp = !Number.isNaN(an) && !Number.isNaN(bn) ? an - bn : av.localeCompare(bv);
      return ascending ? cmp : -cmp;
    }});
    th.dataset.sort = ascending ? 'asc' : 'desc';
    rows.forEach(row => tbody.appendChild(row));
  }});
}});
</script>
</body>
</html>
"""
    report_path = output_dir / "report.html"
    report_path.write_text(html, encoding="utf-8")
    return report_path


def build_single_tool_report(output_dir, prefix, tool_name, bed_path, preds, overlap_threshold=60, max_tracks=10):
    prediction_rows, metric_row = metrics.evaluate_tool(tool_name, bed_path, preds, overlap_threshold=overlap_threshold)
    metrics.write_csv_dicts(output_dir / f"{prefix}_per_prediction_results.csv", prediction_rows)
    metrics.write_csv_dicts(output_dir / f"{prefix}_tool_metrics.csv", [metric_row])
    codon_rows = metrics.load_codon_summary(output_dir / f"{prefix}_codon_summary.csv", tool_name)
    return build_report(
        output_dir,
        f"rORForise Report: {tool_name}",
        [metric_row],
        prediction_rows,
        codon_rows=codon_rows,
        ranking_rows=metrics.rank_tools([metric_row]),
        stratified_rows=metrics.stratified_metric_rows(prediction_rows),
        boundary_rows=metrics.boundary_diagnostic_rows(prediction_rows),
        scorecard_rows=[{"tool": tool_name, "tool_type": metric_row.get("tool_type", "unspecified"), "coordinate_type": "", "interpretation": ""}],
        max_tracks=max_tracks,
    )
