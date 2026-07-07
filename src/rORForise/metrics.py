import csv
import math
import os
from collections import defaultdict
from pathlib import Path

try:
    from . import check_pred as cp
    from .adapters import genome_interval_for_prediction, prediction_overlaps_cds, read_bed_records
except Exception:
    import check_pred as cp
    from adapters import genome_interval_for_prediction, prediction_overlaps_cds, read_bed_records


ANSWER_NAMES = tuple(cp.answers.keys())
CONTEXTS = ("start", "stop", "middle")


def safe_name(name):
    cleaned = []
    for char in str(name):
        cleaned.append(char if char.isalnum() or char in ("-", "_", ".") else "_")
    return "".join(cleaned).strip("_") or "tool"


def safe_div(num, den):
    return float(num) / float(den) if den else None


def pct(num, den):
    value = safe_div(num, den)
    return value * 100.0 if value is not None else None


def wilson_ci_pct(successes, total, z=1.96):
    if not total:
        return None, None
    phat = successes / total
    denom = 1 + z * z / total
    center = (phat + z * z / (2 * total)) / denom
    spread = z * math.sqrt((phat * (1 - phat) + z * z / (4 * total)) / total) / denom
    return max(0.0, (center - spread) * 100.0), min(100.0, (center + spread) * 100.0)


def format_pct(value):
    return "N/A" if value is None else f"{value:.2f}%"


def parse_pct(value):
    if value in (None, "", "N/A"):
        return None
    return float(str(value).rstrip("%"))


def read_csv_dicts(path):
    if not path or not os.path.exists(path):
        return []
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def write_csv_dicts(path, rows, fieldnames=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = []
        for row in rows:
            for key in row.keys():
                if key not in fieldnames:
                    fieldnames.append(key)
    formatted_rows = [_format_csv_row(row) for row in rows]
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(formatted_rows)


def _format_csv_row(row):
    formatted = dict(row)
    for key, value in formatted.items():
        if key.endswith("_pct") or key == "composite_score":
            formatted[key] = _format_decimal(value)
    return formatted


def _format_decimal(value):
    if value in (None, "", "N/A"):
        return ""
    try:
        return f"{float(value):.2f}"
    except (TypeError, ValueError):
        return value


def read_metric_csv(path):
    return {row["metric"]: row["value"] for row in read_csv_dicts(path)}


def _answer_names(answer_details):
    return {cp.inverse_answers[answer] for answer, _ in answer_details}


def _answer_codons(answer_details):
    codons = {}
    for answer, codon in answer_details:
        if codon:
            codons[cp.inverse_answers[answer]] = codon
    return codons


def _cds_interval_on_read(record):
    if record.read_strand == "-":
        start = record.read_end - record.cds_end + 1
        end = record.read_end - record.cds_start + 1
    else:
        start = record.cds_start - record.read_start + 1
        end = record.cds_end - record.read_start + 1
    start = max(1, min(record.read_length, start))
    end = max(1, min(record.read_length, end))
    return min(start, end), max(start, end)


def _boundary_offsets(record, pred_start, pred_end, pred_strand):
    category = (record.cds_strand, record.read_strand, pred_strand)
    if category == ("+", "+", "+"):
        pred_cds_start = record.read_start + (pred_start - 1)
        pred_cds_end = record.read_start + (pred_end - 1)
        return pred_cds_start - record.cds_start, pred_cds_end - record.cds_end
    if category == ("+", "-", "-"):
        pred_cds_start = record.read_end - (pred_end - 1)
        pred_cds_end = record.read_end - (pred_start - 1)
        return pred_cds_start - record.cds_start, pred_cds_end - record.cds_end
    if category == ("-", "+", "-"):
        pred_cds_start = record.read_start + (pred_start - 1)
        pred_cds_end = record.read_start + (pred_end - 1)
        return pred_cds_end - record.cds_end, pred_cds_start - record.cds_start
    if category == ("-", "-", "+"):
        pred_cds_start = record.read_end - (pred_start - 1)
        pred_cds_end = record.read_end - (pred_end - 1)
        return pred_cds_start - record.cds_end, pred_cds_end - record.cds_start

    pred_genome_start, pred_genome_end = genome_interval_for_prediction(record, pred_start, pred_end)
    cds_bio_start = record.cds_start if record.cds_strand == "+" else record.cds_end
    cds_bio_stop = record.cds_end if record.cds_strand == "+" else record.cds_start
    nearest_to_start = min(pred_genome_start - cds_bio_start, pred_genome_end - cds_bio_start, key=abs)
    nearest_to_stop = min(pred_genome_start - cds_bio_stop, pred_genome_end - cds_bio_stop, key=abs)
    return nearest_to_start, nearest_to_stop


def _boundary_on_read(record, pred_start, pred_end, pred_strand, boundary):
    category = (record.cds_strand, record.read_strand, pred_strand)
    if boundary == "start":
        if category in {("+", "+", "+"), ("-", "-", "+")}:
            return pred_start
        if category in {("+", "-", "-"), ("-", "+", "-")}:
            return pred_end
    if boundary == "stop":
        if category in {("+", "+", "+"), ("-", "-", "+")}:
            return pred_end
        if category in {("+", "-", "-"), ("-", "+", "-")}:
            return pred_start
    return None


def _edge_class(boundary_on_read, read_length):
    if boundary_on_read is None:
        return "wrong_direction"
    if boundary_on_read in (1, read_length):
        return "exact_read_edge"
    if boundary_on_read <= 3 or boundary_on_read >= read_length - 2:
        return "near_read_edge_3nt"
    return "internal_boundary"


def _cds_biological_start(record):
    return record.cds_start if record.cds_strand == "+" else record.cds_end


def _cds_biological_stop(record):
    return record.cds_end if record.cds_strand == "+" else record.cds_start


def _overlap_metrics(record, pred_start, pred_end):
    pred_genome_start, pred_genome_end = genome_interval_for_prediction(record, pred_start, pred_end)
    overlap_bp = max(0, min(pred_genome_end, record.cds_end) - max(pred_genome_start, record.cds_start) + 1)
    pred_len = pred_genome_end - pred_genome_start + 1
    cds_len = record.cds_end - record.cds_start + 1
    return {
        "pred_genome_start": pred_genome_start,
        "pred_genome_end": pred_genome_end,
        "prediction_cds_overlap_bp": overlap_bp,
        "prediction_overlap_pct": pct(overlap_bp, pred_len),
        "cds_covered_pct": pct(overlap_bp, cds_len),
        "reciprocal_overlap_80": int(pred_len > 0 and cds_len > 0 and overlap_bp / pred_len >= 0.8 and overlap_bp / cds_len >= 0.8),
        "prediction_fully_inside_cds": int(pred_genome_start >= record.cds_start and pred_genome_end <= record.cds_end),
        "cds_fully_inside_prediction": int(pred_genome_start <= record.cds_start and pred_genome_end >= record.cds_end),
        "extends_before_cds": int(pred_genome_start < record.cds_start),
        "extends_after_cds": int(pred_genome_end > record.cds_end),
    }


def _read_context(record):
    if record.captures_cds_start and record.captures_cds_end:
        return "captures_start_and_stop"
    if record.captures_cds_start:
        return "captures_start"
    if record.captures_cds_end:
        return "captures_stop"
    return "middle_only"


def build_prediction_rows(tool_name, bed_path, preds, overlap_threshold=60, read_sequences=None, tool_type="unspecified"):
    rows = []
    records = [record for record in read_bed_records(bed_path, read_sequences=read_sequences) if record.feature_type == "CDS"]

    for record in records:
        if record.overlap_bp < overlap_threshold or record.read_name not in preds:
            continue
        cds_read_start, cds_read_end = _cds_interval_on_read(record)
        for pred_id, (pred_start, pred_end, pred_strand) in preds[record.read_name].items():
            answers = cp.check_pred(
                record.cds_start,
                record.cds_end,
                record.cds_strand,
                record.read_start,
                record.read_end,
                record.read_strand,
                pred_start,
                pred_end,
                pred_strand,
                record.read_seq,
            )
            answer_names = _answer_names(answers)
            answer_codons = _answer_codons(answers)
            start_offset, stop_offset = _boundary_offsets(record, pred_start, pred_end, pred_strand)
            middle_only = not record.captures_cds_start and not record.captures_cds_end
            on_target = prediction_overlaps_cds(record, pred_start, pred_end)
            start_boundary = _boundary_on_read(record, pred_start, pred_end, pred_strand, "start")
            stop_boundary = _boundary_on_read(record, pred_start, pred_end, pred_strand, "stop")
            overlap_details = _overlap_metrics(record, pred_start, pred_end)

            row = {
                "tool": tool_name,
                "tool_type": tool_type,
                "read_name": record.read_name,
                "prediction_id": pred_id,
                "pred_start": pred_start,
                "pred_end": pred_end,
                "pred_strand": pred_strand,
                "pred_length": pred_end - pred_start + 1,
                "read_start": record.read_start,
                "read_end": record.read_end,
                "read_length": record.read_length,
                "read_strand": record.read_strand,
                "cds_start": record.cds_start,
                "cds_end": record.cds_end,
                "cds_strand": record.cds_strand,
                "cds_read_start": cds_read_start,
                "cds_read_end": cds_read_end,
                "read_cds_overlap_bp": record.overlap_bp,
                "captures_start": int(record.captures_cds_start),
                "captures_stop": int(record.captures_cds_end),
                "middle_only": int(middle_only),
                "read_context": _read_context(record),
                "on_target": int(on_target),
                "on_target_category": "overlaps_cds" if on_target else "no_cds_overlap",
                "start_offset": start_offset,
                "stop_offset": stop_offset,
                "start_boundary_on_read": start_boundary if start_boundary is not None else "",
                "stop_boundary_on_read": stop_boundary if stop_boundary is not None else "",
                "start_boundary_class": _edge_class(start_boundary, record.read_length),
                "stop_boundary_class": _edge_class(stop_boundary, record.read_length),
                "cds_start_at_read_edge": int(_cds_biological_start(record) in {record.read_start, record.read_end}),
                "cds_stop_at_read_edge": int(_cds_biological_stop(record) in {record.read_start, record.read_end}),
                "answers": "; ".join(sorted(answer_names)),
                "start_codon": answer_codons.get("correct start")
                or answer_codons.get("alternative start")
                or answer_codons.get("incorrect start")
                or answer_codons.get("middle incorrect start")
                or answer_codons.get("middle or alternative start")
                or "",
                "stop_codon": answer_codons.get("correct stop")
                or answer_codons.get("alternative stop")
                or answer_codons.get("incorrect stop")
                or answer_codons.get("middle incorrect stop")
                or answer_codons.get("middle or alternative stop")
                or "",
            }
            row.update(overlap_details)
            for answer_name in ANSWER_NAMES:
                row[answer_name.replace(" ", "_")] = int(answer_name in answer_names)
            rows.append(row)
    return rows


def summarize_tool(tool_name, bed_path, preds, rows, overlap_threshold=60, read_sequences=None, tool_type="unspecified"):
    records = read_bed_records(bed_path, read_sequences=read_sequences)
    cds_records = [record for record in records if record.feature_type == "CDS"]
    unique_reads = {record.read_name for record in records}
    cds_reads = {record.read_name for record in cds_records}
    good_records = [record for record in cds_records if record.overlap_bp >= overlap_threshold]
    good_reads = {record.read_name for record in good_records}
    reads_without_predictions = {record.read_name for record in good_records if record.read_name not in preds}

    inspected = len(rows)
    start_rows = [row for row in rows if int(row["captures_start"])]
    stop_rows = [row for row in rows if int(row["captures_stop"])]
    middle_rows = [row for row in rows if int(row["middle_only"])]
    both_rows = [row for row in rows if int(row["captures_start"]) and int(row["captures_stop"])]

    def count(flag, chosen_rows=rows):
        return sum(int(row.get(flag, 0)) for row in chosen_rows)

    exact_orf = sum(int(row.get("correct_start", 0)) and int(row.get("correct_stop", 0)) for row in both_rows)
    start_within_3 = sum(abs(int(row["start_offset"])) <= 3 for row in start_rows)
    start_within_9 = sum(abs(int(row["start_offset"])) <= 9 for row in start_rows)
    stop_within_3 = sum(abs(int(row["stop_offset"])) <= 3 for row in stop_rows)
    stop_within_9 = sum(abs(int(row["stop_offset"])) <= 9 for row in stop_rows)
    on_target = sum(int(row["on_target"]) for row in rows)
    reciprocal_80 = count("reciprocal_overlap_80")
    fully_inside_cds = count("prediction_fully_inside_cds")
    cds_inside_prediction = count("cds_fully_inside_prediction")
    correct_start_exact_edge = sum(int(row.get("correct_start", 0)) for row in start_rows if row.get("start_boundary_class") == "exact_read_edge")
    correct_start_near_edge = sum(int(row.get("correct_start", 0)) for row in start_rows if row.get("start_boundary_class") == "near_read_edge_3nt")
    correct_start_internal = sum(int(row.get("correct_start", 0)) for row in start_rows if row.get("start_boundary_class") == "internal_boundary")

    def average(field):
        values = [float(row[field]) for row in rows if row.get(field) not in (None, "")]
        return sum(values) / len(values) if values else None

    summary = {
        "tool": tool_name,
        "tool_type": tool_type,
        "total_reads": len(unique_reads),
        "cds_aligned_reads": len(cds_reads),
        "reads_good_overlap": len(good_reads),
        "reads_without_predictions": len(reads_without_predictions),
        "total_predictions": sum(len(read_preds) for read_preds in preds.values()),
        "inspected_prediction_cds_comparisons": inspected,
        "on_target_predictions": on_target,
        "off_target_predictions": inspected - on_target,
        "reciprocal_overlap_80": reciprocal_80,
        "prediction_fully_inside_cds": fully_inside_cds,
        "cds_fully_inside_prediction": cds_inside_prediction,
        "start_eligible": len(start_rows),
        "stop_eligible": len(stop_rows),
        "middle_eligible": len(middle_rows),
        "both_boundary_eligible": len(both_rows),
        "correct_start": count("correct_start", start_rows),
        "correct_stop": count("correct_stop", stop_rows),
        "correct_middle": count("middle", middle_rows),
        "alternative_start": count("alternative_start"),
        "alternative_stop": count("alternative_stop"),
        "incorrect_start": count("incorrect_start"),
        "incorrect_stop": count("incorrect_stop"),
        "correct_frame": count("correct_frame"),
        "incorrect_frame": count("incorrect_frame"),
        "correct_direction": count("correct_direction"),
        "incorrect_direction": count("incorrect_direction"),
        "exact_orf_boundary": exact_orf,
        "start_within_3nt": start_within_3,
        "start_within_9nt": start_within_9,
        "stop_within_3nt": stop_within_3,
        "stop_within_9nt": stop_within_9,
        "correct_start_exact_read_edge": correct_start_exact_edge,
        "correct_start_near_read_edge_3nt": correct_start_near_edge,
        "correct_start_internal_boundary": correct_start_internal,
    }

    summary.update(
        {
            "prediction_coverage_pct": pct(len(good_reads) - len(reads_without_predictions), len(good_reads)),
            "on_target_rate_pct": pct(on_target, inspected),
            "off_target_rate_pct": pct(inspected - on_target, inspected),
            "reciprocal_overlap_80_pct": pct(reciprocal_80, inspected),
            "prediction_fully_inside_cds_pct": pct(fully_inside_cds, inspected),
            "cds_fully_inside_prediction_pct": pct(cds_inside_prediction, inspected),
            "avg_prediction_overlap_pct": average("prediction_overlap_pct"),
            "avg_cds_covered_pct": average("cds_covered_pct"),
            "start_accuracy_pct": pct(summary["correct_start"], len(start_rows)),
            "stop_accuracy_pct": pct(summary["correct_stop"], len(stop_rows)),
            "middle_accuracy_pct": pct(summary["correct_middle"], len(middle_rows)),
            "frame_accuracy_pct": pct(summary["correct_frame"], inspected),
            "direction_accuracy_pct": pct(summary["correct_direction"], inspected),
            "exact_orf_boundary_accuracy_pct": pct(exact_orf, len(both_rows)),
            "start_within_3nt_pct": pct(start_within_3, len(start_rows)),
            "start_within_9nt_pct": pct(start_within_9, len(start_rows)),
            "stop_within_3nt_pct": pct(stop_within_3, len(stop_rows)),
            "stop_within_9nt_pct": pct(stop_within_9, len(stop_rows)),
        }
    )
    _add_ci(summary, "start_accuracy", summary["correct_start"], len(start_rows))
    _add_ci(summary, "stop_accuracy", summary["correct_stop"], len(stop_rows))
    _add_ci(summary, "frame_accuracy", summary["correct_frame"], inspected)
    _add_ci(summary, "direction_accuracy", summary["correct_direction"], inspected)
    _add_ci(summary, "on_target_rate", on_target, inspected)
    return summary


def _add_ci(summary, prefix, successes, total):
    low, high = wilson_ci_pct(successes, total)
    summary[f"{prefix}_ci95_low_pct"] = low
    summary[f"{prefix}_ci95_high_pct"] = high


def evaluate_tool(tool_name, bed_path, preds, overlap_threshold=60, read_sequences=None, tool_type="unspecified"):
    rows = build_prediction_rows(
        tool_name,
        bed_path,
        preds,
        overlap_threshold=overlap_threshold,
        read_sequences=read_sequences,
        tool_type=tool_type,
    )
    summary = summarize_tool(
        tool_name,
        bed_path,
        preds,
        rows,
        overlap_threshold=overlap_threshold,
        read_sequences=read_sequences,
        tool_type=tool_type,
    )
    return rows, summary


def stratified_metric_rows(prediction_rows):
    rows = []
    strata = [
        ("read_context", "read_context"),
        ("read_strand", "read_strand"),
        ("cds_strand", "cds_strand"),
        ("pred_strand", "pred_strand"),
        ("start_boundary_class", "start_boundary_class"),
        ("stop_boundary_class", "stop_boundary_class"),
        ("start_codon", "start_codon"),
    ]
    for tool in sorted({row["tool"] for row in prediction_rows}):
        tool_rows = [row for row in prediction_rows if row["tool"] == tool]
        for field, label in strata:
            values = sorted({row.get(field, "") or "blank" for row in tool_rows})
            for value in values:
                chosen = [row for row in tool_rows if (row.get(field, "") or "blank") == value]
                if not chosen:
                    continue
                start_rows = [row for row in chosen if int(row.get("captures_start", 0))]
                stop_rows = [row for row in chosen if int(row.get("captures_stop", 0))]
                rows.append(_summarize_prediction_subset(tool, label, value, chosen, start_rows, stop_rows))
    return rows


def _summarize_prediction_subset(tool, stratification, value, rows, start_rows, stop_rows):
    total = len(rows)
    correct_start = sum(int(row.get("correct_start", 0)) for row in start_rows)
    correct_stop = sum(int(row.get("correct_stop", 0)) for row in stop_rows)
    correct_frame = sum(int(row.get("correct_frame", 0)) for row in rows)
    correct_direction = sum(int(row.get("correct_direction", 0)) for row in rows)
    on_target = sum(int(row.get("on_target", 0)) for row in rows)
    return {
        "tool": tool,
        "stratification": stratification,
        "value": value,
        "comparisons": total,
        "start_eligible": len(start_rows),
        "stop_eligible": len(stop_rows),
        "correct_start": correct_start,
        "correct_stop": correct_stop,
        "correct_frame": correct_frame,
        "correct_direction": correct_direction,
        "on_target": on_target,
        "start_accuracy_pct": pct(correct_start, len(start_rows)),
        "stop_accuracy_pct": pct(correct_stop, len(stop_rows)),
        "frame_accuracy_pct": pct(correct_frame, total),
        "direction_accuracy_pct": pct(correct_direction, total),
        "on_target_rate_pct": pct(on_target, total),
    }


def boundary_diagnostic_rows(prediction_rows):
    rows = []
    for tool in sorted({row["tool"] for row in prediction_rows}):
        tool_rows = [row for row in prediction_rows if row["tool"] == tool]
        for boundary in ("start", "stop"):
            class_field = f"{boundary}_boundary_class"
            correct_field = f"correct_{boundary}"
            eligible_field = "captures_start" if boundary == "start" else "captures_stop"
            for boundary_class in sorted({row.get(class_field, "unknown") for row in tool_rows}):
                chosen = [row for row in tool_rows if row.get(class_field, "unknown") == boundary_class and int(row.get(eligible_field, 0))]
                correct = sum(int(row.get(correct_field, 0)) for row in chosen)
                rows.append(
                    {
                        "tool": tool,
                        "boundary": boundary,
                        "boundary_class": boundary_class,
                        "eligible_comparisons": len(chosen),
                        "correct": correct,
                        "correct_pct": pct(correct, len(chosen)),
                    }
                )
    return rows


def pairwise_read_overlap(tool_preds):
    rows = []
    names = list(tool_preds.keys())
    for i, tool_a in enumerate(names):
        for tool_b in names[i + 1:]:
            reads_a = set(tool_preds[tool_a].keys())
            reads_b = set(tool_preds[tool_b].keys())
            shared_reads = reads_a & reads_b
            rows.append(
                {
                    "tool_a": tool_a,
                    "tool_b": tool_b,
                    "reads_with_preds_a": len(reads_a),
                    "reads_with_preds_b": len(reads_b),
                    "shared_reads": len(shared_reads),
                    "shared_reads_frac_a": safe_div(len(shared_reads), len(reads_a)) or 0,
                    "shared_reads_frac_b": safe_div(len(shared_reads), len(reads_b)) or 0,
                }
            )
    return rows


def _interval_overlap(a_start, a_end, b_start, b_end):
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


def pairwise_prediction_overlap(tool_preds):
    rows = []
    names = list(tool_preds.keys())
    for i, tool_a in enumerate(names):
        for tool_b in names[i + 1:]:
            preds_a = tool_preds[tool_a]
            preds_b = tool_preds[tool_b]
            shared_reads = set(preds_a) & set(preds_b)
            exact_matches = 0
            reciprocal_80 = 0
            same_frame = 0
            same_strand_overlap = 0

            for read_name in shared_reads:
                values_a = list(preds_a[read_name].values())
                values_b = list(preds_b[read_name].values())
                exact_matches += len(set(values_a) & set(values_b))
                for a_start, a_end, a_strand in values_a:
                    a_len = a_end - a_start + 1
                    for b_start, b_end, b_strand in values_b:
                        overlap = _interval_overlap(a_start, a_end, b_start, b_end)
                        if overlap <= 0:
                            continue
                        if a_strand == b_strand:
                            same_strand_overlap += 1
                        b_len = b_end - b_start + 1
                        if overlap / a_len >= 0.8 and overlap / b_len >= 0.8:
                            reciprocal_80 += 1
                        if a_strand == b_strand and (a_start - b_start) % 3 == 0:
                            same_frame += 1

            rows.append(
                {
                    "tool_a": tool_a,
                    "tool_b": tool_b,
                    "shared_reads": len(shared_reads),
                    "exact_pred_matches": exact_matches,
                    "reciprocal_overlap_80": reciprocal_80,
                    "same_strand_overlaps": same_strand_overlap,
                    "same_frame_overlaps": same_frame,
                    "total_preds_a": sum(len(v) for v in preds_a.values()),
                    "total_preds_b": sum(len(v) for v in preds_b.values()),
                }
            )
    return rows


def rank_tools(metric_rows):
    ranked = []
    score_fields = (
        "start_accuracy_pct",
        "stop_accuracy_pct",
        "frame_accuracy_pct",
        "direction_accuracy_pct",
        "on_target_rate_pct",
    )
    for row in metric_rows:
        values = [row.get(field) for field in score_fields if row.get(field) is not None]
        score = sum(values) / len(values) if values else None
        ranked.append(
            {
                "rank": 0,
                "tool": row["tool"],
                "tool_type": row.get("tool_type", "unspecified"),
                "composite_score": score,
                "start_accuracy_pct": row.get("start_accuracy_pct"),
                "stop_accuracy_pct": row.get("stop_accuracy_pct"),
                "frame_accuracy_pct": row.get("frame_accuracy_pct"),
                "direction_accuracy_pct": row.get("direction_accuracy_pct"),
                "on_target_rate_pct": row.get("on_target_rate_pct"),
            }
        )
    ranked.sort(key=lambda row: -1 if row["composite_score"] is None else row["composite_score"], reverse=True)
    for index, row in enumerate(ranked, 1):
        row["rank"] = index
    return ranked


def load_codon_summary(path, tool_name):
    rows = []
    for row in read_csv_dicts(path):
        row = dict(row)
        row["tool"] = tool_name
        rows.append(row)
    return rows
