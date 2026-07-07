from pathlib import Path
from collections import defaultdict, deque
import os


def _pyplot():
    cache_dir = Path(os.environ.get("MPLCONFIGDIR", "/private/tmp/rorforise_matplotlib"))
    cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir))
    xdg_cache = cache_dir / "xdg"
    xdg_cache.mkdir(exist_ok=True)
    os.environ.setdefault("XDG_CACHE_HOME", str(xdg_cache))

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def _float(row, key):
    value = row.get(key)
    if value in (None, "", "N/A"):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _save_placeholder(path, title, message="No data available"):
    plt = _pyplot()
    fig, ax = plt.subplots(figsize=(8, 3))
    ax.axis("off")
    ax.set_title(title)
    ax.text(0.5, 0.5, message, ha="center", va="center", fontsize=12)
    fig.tight_layout()
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path


def plot_context_accuracy(metric_rows, path):
    plt = _pyplot()
    contexts = [
        ("start_accuracy_pct", "Start"),
        ("stop_accuracy_pct", "Stop"),
        ("middle_accuracy_pct", "Middle"),
    ]
    rows = list(metric_rows)
    if not rows:
        return _save_placeholder(path, "Context Accuracy")

    labels = [row["tool"] for row in rows]
    x = range(len(labels))
    width = 0.25
    fig, ax = plt.subplots(figsize=(max(7, len(labels) * 1.2), 4.5))
    colors = ["#2864a6", "#2d8b57", "#b45f06"]

    for offset, (field, label) in enumerate(contexts):
        values = [_float(row, field) or 0 for row in rows]
        ax.bar([i + (offset - 1) * width for i in x], values, width=width, label=label, color=colors[offset])

    ax.set_ylabel("Accuracy (%)")
    ax.set_title("Context Accuracy")
    ax.set_ylim(0, 100)
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_start_stop_tolerance(metric_rows, path):
    plt = _pyplot()
    rows = list(metric_rows)
    if not rows:
        return _save_placeholder(path, "Start/Stop Boundary Accuracy")

    labels = [row["tool"] for row in rows]
    fields = [
        ("start_accuracy_pct", "Start exact", "#2864a6"),
        ("start_within_9nt_pct", "Start <=9 nt", "#7aa6d8"),
        ("stop_accuracy_pct", "Stop exact", "#2d8b57"),
        ("stop_within_9nt_pct", "Stop <=9 nt", "#8ccf9f"),
    ]
    x = range(len(labels))
    width = 0.2
    fig, ax = plt.subplots(figsize=(max(8, len(labels) * 1.4), 4.7))
    for offset, (field, label, color) in enumerate(fields):
        values = [_float(row, field) or 0 for row in rows]
        ax.bar([i + (offset - 1.5) * width for i in x], values, width=width, label=label, color=color)
    ax.set_ylabel("Predictions (%)")
    ax.set_title("Boundary Accuracy and Tolerance")
    ax.set_ylim(0, 100)
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.legend(frameon=False, ncols=2)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_frame_direction(metric_rows, path):
    plt = _pyplot()
    rows = list(metric_rows)
    if not rows:
        return _save_placeholder(path, "Frame and Direction Accuracy")

    labels = [row["tool"] for row in rows]
    x = range(len(labels))
    width = 0.35
    fig, ax = plt.subplots(figsize=(max(7, len(labels) * 1.2), 4.5))
    ax.bar([i - width / 2 for i in x], [_float(row, "frame_accuracy_pct") or 0 for row in rows], width=width, label="Frame", color="#5f6caf")
    ax.bar([i + width / 2 for i in x], [_float(row, "direction_accuracy_pct") or 0 for row in rows], width=width, label="Direction", color="#8a5a83")
    ax.set_ylabel("Accuracy (%)")
    ax.set_title("Frame and Direction Accuracy")
    ax.set_ylim(0, 100)
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_on_target(metric_rows, path):
    plt = _pyplot()
    rows = list(metric_rows)
    if not rows:
        return _save_placeholder(path, "On-Target Rate")

    labels = [row["tool"] for row in rows]
    values = [_float(row, "on_target_rate_pct") or 0 for row in rows]
    fig, ax = plt.subplots(figsize=(max(7, len(labels) * 1.1), 4.2))
    ax.bar(labels, values, color="#3d7c78")
    ax.set_ylabel("On-target comparisons (%)")
    ax.set_title("On-Target Prediction Rate")
    ax.set_ylim(0, 100)
    ax.tick_params(axis="x", rotation=30)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_prediction_lengths(prediction_rows, path):
    plt = _pyplot()
    rows = list(prediction_rows)
    if not rows:
        return _save_placeholder(path, "Prediction Lengths")

    tools = sorted({row["tool"] for row in rows})
    fig, ax = plt.subplots(figsize=(8, 4.8))
    for tool in tools:
        lengths = [int(row["pred_length"]) for row in rows if row["tool"] == tool]
        if lengths:
            ax.hist(lengths, bins=min(30, max(5, len(set(lengths)))), alpha=0.45, label=tool)
    ax.set_xlabel("Prediction length (nt)")
    ax.set_ylabel("Count")
    ax.set_title("Prediction Length Distribution")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_codon_enrichment(codon_rows, path):
    plt = _pyplot()
    rows = [
        row
        for row in codon_rows
        if row.get("answer_type") in {"incorrect start", "incorrect stop", "middle incorrect start", "middle incorrect stop"}
        and row.get("obs_exp_ratio") not in ("", "N/A", None)
    ]
    if not rows:
        return _save_placeholder(path, "Codon Enrichment Heatmap")

    codons = sorted({row["codon"] for row in rows})
    labels = sorted({f"{row.get('tool', 'tool')} / {row['answer_type']}" for row in rows})
    matrix = []
    lookup = {(f"{row.get('tool', 'tool')} / {row['answer_type']}", row["codon"]): float(row["obs_exp_ratio"]) for row in rows}
    for label in labels:
        matrix.append([lookup.get((label, codon), 0.0) for codon in codons])

    fig, ax = plt.subplots(figsize=(max(10, len(codons) * 0.18), max(4, len(labels) * 0.35)))
    image = ax.imshow(matrix, aspect="auto", cmap="viridis")
    ax.set_xticks(range(len(codons)))
    ax.set_xticklabels(codons, rotation=90, fontsize=7)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_title("Incorrect Boundary Codon Enrichment")
    cbar = fig.colorbar(image, ax=ax)
    cbar.set_label("Observed / expected")
    fig.tight_layout()
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def _track_interest_score(row):
    return (
        int(row.get("correct_start", 0))
        + int(row.get("correct_stop", 0))
        + int(row.get("correct_direction", 0))
    )


def _select_track_rows(prediction_rows, max_examples):
    candidates = sorted(
        prediction_rows,
        key=lambda row: (
            _track_interest_score(row),
            row["tool"],
            row["read_name"],
        ),
    )
    by_tool = defaultdict(deque)
    for row in candidates:
        by_tool[row["tool"]].append(row)

    tools = sorted(by_tool)
    chosen = []
    while len(chosen) < max_examples and any(by_tool.values()):
        for tool in tools:
            if by_tool[tool]:
                chosen.append(by_tool[tool].popleft())
                if len(chosen) >= max_examples:
                    break
    return chosen


def plot_read_tracks(prediction_rows, output_dir, max_examples=10):
    plt = _pyplot()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    for old_track in output_dir.glob("track_*.png"):
        old_track.unlink()
    chosen = _select_track_rows(prediction_rows, max_examples)
    paths = []
    for index, row in enumerate(chosen, 1):
        read_length = int(row["read_length"])
        pred_start = int(row["pred_start"])
        pred_end = int(row["pred_end"])
        cds_start = int(row["cds_read_start"])
        cds_end = int(row["cds_read_end"])

        fig, ax = plt.subplots(figsize=(8, 1.9))
        ax.hlines(1.0, 1, read_length, color="#333333", linewidth=2)
        ax.plot([cds_start, cds_end], [1.18, 1.18], color="#2d8b57", linewidth=8, solid_capstyle="butt")
        ax.plot([pred_start, pred_end], [0.82, 0.82], color="#2864a6", linewidth=8, solid_capstyle="butt")
        ax.set_xlim(1, read_length)
        ax.set_ylim(0.55, 1.45)
        ax.set_yticks([1.18, 1.0, 0.82])
        ax.set_yticklabels(["CDS", "Read", "Prediction"], fontsize=8)
        ax.set_xlabel("Read coordinate (nt)")
        ax.set_title(f"{row['tool']} | {row['read_name']} | {row['answers']}", fontsize=9)
        fig.tight_layout()
        path = output_dir / f"track_{index:03d}_{row['tool']}_{row['read_name'].replace('/', '_')}.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        paths.append(path)
    return paths


def plot_multitool_read_tracks(prediction_rows, output_dir, max_reads=8):
    plt = _pyplot()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    for old_track in output_dir.glob("comparison_*.png"):
        old_track.unlink()

    grouped = defaultdict(list)
    for row in prediction_rows:
        grouped[row["read_name"]].append(row)

    def group_score(item):
        read_name, rows = item
        tools = {row["tool"] for row in rows}
        correctness = sum(_track_interest_score(row) for row in rows) / max(1, len(rows))
        return (-len(tools), correctness, read_name)

    chosen_groups = [
        (read_name, rows)
        for read_name, rows in sorted(grouped.items(), key=group_score)
        if len({row["tool"] for row in rows}) > 1
    ][:max_reads]

    paths = []
    for index, (read_name, rows) in enumerate(chosen_groups, 1):
        first = rows[0]
        read_length = int(first["read_length"])
        cds_start = int(first["cds_read_start"])
        cds_end = int(first["cds_read_end"])
        tools = sorted({row["tool"] for row in rows})
        y_by_tool = {tool: i + 1 for i, tool in enumerate(reversed(tools))}
        cds_y = len(tools) + 1

        fig, ax = plt.subplots(figsize=(9, max(2.4, 0.36 * len(tools) + 1.2)))
        ax.hlines(cds_y - 0.35, 1, read_length, color="#333333", linewidth=1.5)
        ax.plot([cds_start, cds_end], [cds_y, cds_y], color="#2d8b57", linewidth=7, solid_capstyle="butt")

        for row in rows:
            y = y_by_tool[row["tool"]]
            ax.plot(
                [int(row["pred_start"]), int(row["pred_end"])],
                [y, y],
                color="#2864a6" if int(row.get("correct_direction", 0)) else "#a6423a",
                linewidth=6,
                alpha=0.85,
                solid_capstyle="butt",
            )

        ax.set_xlim(1, read_length)
        ax.set_ylim(0.4, cds_y + 0.6)
        ax.set_yticks([cds_y] + [y_by_tool[tool] for tool in tools])
        ax.set_yticklabels(["CDS"] + tools, fontsize=8)
        ax.set_xlabel("Read coordinate (nt)")
        ax.set_title(f"Multi-tool comparison | {read_name}", fontsize=10)
        ax.plot([], [], color="#2864a6", linewidth=6, label="Correct direction")
        ax.plot([], [], color="#a6423a", linewidth=6, label="Incorrect direction")
        ax.legend(loc="upper right", frameon=False, fontsize=8, ncols=2)
        fig.tight_layout()
        safe_read = "".join(char if char.isalnum() or char in "-_." else "_" for char in read_name)
        path = output_dir / f"comparison_{index:03d}_{safe_read}.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        paths.append(path)
    return paths
