import csv
from dataclasses import dataclass
import gzip


@dataclass(frozen=True)
class BedRecord:
    chrom: str
    read_name: str
    read_start: int
    read_end: int
    read_strand: str
    feature_type: str
    cds_start: int
    cds_end: int
    cds_strand: str
    read_seq: str
    score: str = "."
    gene_id: str | None = None

    @property
    def read_length(self):
        return self.read_end - self.read_start + 1

    @property
    def overlap_bp(self):
        return max(0, min(self.read_end, self.cds_end) - max(self.read_start, self.cds_start) + 1)

    @property
    def captures_cds_start(self):
        if self.cds_strand == "-":
            return self.read_end >= self.cds_end
        return self.read_start <= self.cds_start

    @property
    def captures_cds_end(self):
        if self.cds_strand == "-":
            return self.read_start <= self.cds_start
        return self.read_end >= self.cds_end


@dataclass(frozen=True)
class GffPrediction:
    read_name: str
    pred_id: str
    start: int
    end: int
    strand: str
    seqid: str
    source: str
    feature_type: str
    score: str
    phase: str
    attributes: dict[str, str]


def open_text(path):
    return gzip.open(path, "rt", encoding="utf-8") if str(path).endswith(".gz") else open(path, "r", encoding="utf-8")


def parse_gff_attributes(attributes):
    parsed = {}
    for field in attributes.split(";"):
        field = field.strip()
        if not field:
            continue
        if "=" in field:
            key, value = field.split("=", 1)
        elif " " in field:
            key, value = field.split(" ", 1)
        else:
            parsed.setdefault("ID", field)
            continue
        parsed[key.strip()] = value.strip().strip('"')
    return parsed


def _as_int(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _clean_read_name(name):
    return str(name).lstrip("@")


def read_fasta_sequences(path):
    sequences = {}
    if not path:
        return sequences
    current = None
    chunks = []
    with open_text(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current is not None:
                    sequences[_clean_read_name(current)] = "".join(chunks).upper()
                current = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if current is not None:
            sequences[_clean_read_name(current)] = "".join(chunks).upper()
    return sequences


def _lookup_read_sequence(read_sequences, read_name):
    if not read_sequences:
        return None
    return read_sequences.get(read_name) or read_sequences.get(_clean_read_name(read_name)) or read_sequences.get(f"@{read_name}")


def parse_bed_row(row, read_sequences=None):
    if len(row) < 11:
        return None

    compact = _parse_compact_row(row)
    if compact is not None:
        return compact

    return _parse_bed12_intersect_row(row, read_sequences=read_sequences)


def _parse_compact_row(row):
    read_start = _as_int(row[2])
    read_end = _as_int(row[3])
    cds_start = _as_int(row[7])
    cds_end = _as_int(row[8])
    if None in (read_start, read_end, cds_start, cds_end):
        return None

    gene_id = None
    if len(row) > 11 and row[11]:
        gene_id = row[11]

    return BedRecord(
        chrom=row[0],
        read_name=_clean_read_name(row[1]),
        read_start=read_start,
        read_end=read_end,
        read_strand=row[4],
        score=row[5],
        feature_type=row[6],
        cds_start=cds_start,
        cds_end=cds_end,
        cds_strand=row[9],
        read_seq=row[10],
        gene_id=gene_id,
    )


def _parse_bed12_intersect_row(row, read_sequences=None):
    if len(row) < 21:
        return None

    read_start_0 = _as_int(row[1])
    read_end = _as_int(row[2])
    cds_start = _as_int(row[15])
    cds_end = _as_int(row[16])
    if None in (read_start_0, read_end, cds_start, cds_end):
        return None

    read_name = _clean_read_name(row[3])
    read_seq = _lookup_read_sequence(read_sequences, read_name)
    if not read_seq:
        return None

    gene_id = None
    if len(row) > 20 and row[20]:
        attrs = parse_gff_attributes(row[20])
        gene_id = attrs.get("ID") or attrs.get("Name") or row[20]

    return BedRecord(
        chrom=row[0],
        read_name=read_name,
        read_start=read_start_0 + 1,
        read_end=read_end,
        read_strand=row[5],
        score=row[4],
        feature_type=row[14],
        cds_start=cds_start,
        cds_end=cds_end,
        cds_strand=row[18],
        read_seq=read_seq,
        gene_id=gene_id,
    )


def read_bed_records(intersect_bed_path, read_sequences=None):
    records = []
    with open_text(intersect_bed_path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            record = parse_bed_row(row, read_sequences=read_sequences)
            if record is not None:
                records.append(record)
    return records


def load_bed_map(intersect_bed_path, read_sequences=None):
    bed_map = {}
    for record in read_bed_records(intersect_bed_path, read_sequences=read_sequences):
        bed_map.setdefault(record.read_name, []).append(record)
    return bed_map


def genome_interval_for_prediction(bed_record, pred_start, pred_end):
    if bed_record.read_strand == "-":
        start = bed_record.read_end - pred_end + 1
        end = bed_record.read_end - pred_start + 1
    else:
        start = bed_record.read_start + pred_start - 1
        end = bed_record.read_start + pred_end - 1
    return min(start, end), max(start, end)


def prediction_overlaps_cds(bed_record, pred_start, pred_end):
    pred_genome_start, pred_genome_end = genome_interval_for_prediction(bed_record, pred_start, pred_end)
    return min(pred_genome_end, bed_record.cds_end) - max(pred_genome_start, bed_record.cds_start) + 1 > 0


def _unique_id(existing, pred_id):
    if pred_id not in existing:
        return pred_id
    index = 2
    while f"{pred_id}_{index}" in existing:
        index += 1
    return f"{pred_id}_{index}"


def read_gff_predictions(gff_path):
    predictions = []
    with open_text(gff_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue

            try:
                start = int(cols[3])
                end = int(cols[4])
            except ValueError:
                continue

            attrs = parse_gff_attributes(cols[8])
            pred_id = attrs.get("ID") or attrs.get("Name") or f"{cols[0]}:{start}-{end}:{cols[6]}"
            read_name = cols[0].lstrip("@")
            predictions.append(
                GffPrediction(
                    read_name=read_name,
                    pred_id=pred_id.replace("@", ""),
                    start=start,
                    end=end,
                    strand=cols[6],
                    seqid=cols[0],
                    source=cols[1],
                    feature_type=cols[2],
                    score=cols[5],
                    phase=cols[7],
                    attributes=attrs,
                )
            )
    return predictions


def predictions_to_dict(predictions):
    preds = {}
    for prediction in predictions:
        read_preds = preds.setdefault(prediction.read_name, {})
        pred_id = _unique_id(read_preds, prediction.pred_id)
        read_preds[pred_id] = (prediction.start, prediction.end, prediction.strand)
    return preds


def _detect_coordinate_type(start, end, bed_record):
    read_len = bed_record.read_length
    if 1 <= start <= read_len and 1 <= end <= read_len:
        return "read"
    if bed_record.read_start <= start <= bed_record.read_end or bed_record.read_start <= end <= bed_record.read_end:
        return "genome"
    return "read" if end <= int(read_len * 1.5) else "genome"


def _flip_strand(strand):
    if strand == "+":
        return "-"
    if strand == "-":
        return "+"
    return strand


def _map_prediction_to_read(prediction, bed_record, coord_type):
    read_len = bed_record.read_length
    strand = prediction.strand
    if coord_type == "genome":
        if bed_record.read_strand == "-":
            pstart = bed_record.read_end - prediction.end + 1
            pend = bed_record.read_end - prediction.start + 1
            # The read sequence is stored in the read's own (reverse-complemented)
            # orientation, so a genome-coordinate prediction's strand must be
            # flipped to be meaningful relative to that read-space orientation.
            strand = _flip_strand(strand)
        else:
            pstart = prediction.start - bed_record.read_start + 1
            pend = prediction.end - bed_record.read_start + 1
    else:
        pstart = prediction.start
        pend = prediction.end

    pstart = max(1, min(read_len, pstart))
    pend = max(1, min(read_len, pend))
    return min(pstart, pend), max(pstart, pend), strand


def gff_genome_to_read_preds(gff_path, intersect_bed_path, coordinate_type="auto", read_sequences=None):
    bed_map = load_bed_map(intersect_bed_path, read_sequences=read_sequences)
    preds = {}

    for prediction in read_gff_predictions(gff_path):
        if prediction.read_name not in bed_map:
            continue
        bed_record = bed_map[prediction.read_name][0]
        chosen_type = coordinate_type
        if coordinate_type == "auto":
            chosen_type = _detect_coordinate_type(prediction.start, prediction.end, bed_record)

        read_preds = preds.setdefault(prediction.read_name, {})
        pred_id = _unique_id(read_preds, prediction.pred_id)
        read_preds[pred_id] = _map_prediction_to_read(prediction, bed_record, chosen_type)
    return preds
