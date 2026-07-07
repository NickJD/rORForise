
import argparse
import csv
from pathlib import Path

# robust import of check_pred
try:
    from rORForise import check_pred as cp
except Exception:
    try:
        import check_pred as cp
    except Exception:
        import importlib.util, os
        repo_root = Path(__file__).resolve().parents[2]
        candidate = repo_root / 'src' / 'rORForise' / 'check_pred.py'
        if candidate.exists():
            spec = importlib.util.spec_from_file_location('check_pred', str(candidate))
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            cp = mod
        else:
            raise

try:
    from .adapters import read_bed_records as parse_bed_records, read_gff_predictions
except Exception:
    try:
        from adapters import read_bed_records as parse_bed_records, read_gff_predictions
    except Exception:
        import importlib.util
        candidate = Path(__file__).resolve().parent / 'adapters.py'
        spec = importlib.util.spec_from_file_location('rORForise_adapters', str(candidate))
        adapters = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(adapters)
        parse_bed_records = adapters.read_bed_records
        read_gff_predictions = adapters.read_gff_predictions


def group_bed_records(bed_path):
    bed = {}
    for record in parse_bed_records(bed_path):
        bed.setdefault(record.read_name, []).append(record)
    return bed


def write_validation(out_path, rows):
    header = ['read_name', 'read_start', 'read_end', 'read_strand',
              'pred_id', 'pred_start', 'pred_end', 'pred_strand',
              'cds_id', 'cds_start', 'cds_end', 'cds_strand',
              'pred_scenario', 'should_overlap_start', 'should_overlap_stop',
              'is_middle_only', 'expected_answers', 'pred_start_codon', 'pred_stop_codon']
    with open(out_path, 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for r in rows:
            w.writerow(r)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('-b', '--bed', required=True)
    p.add_argument('-g', '--gff', required=True)
    p.add_argument('-o', '--out', required=True)
    args = p.parse_args()

    bed = group_bed_records(args.bed)
    preds = read_gff_predictions(args.gff)

    rows = []
    for prediction in preds:
        entries = [entry for entry in bed.get(prediction.read_name, []) if entry.feature_type == 'CDS']
        if not entries:
            continue
        for chosen in entries:
            answers = cp.check_pred(chosen.cds_start, chosen.cds_end, chosen.cds_strand,
                                    chosen.read_start, chosen.read_end, chosen.read_strand,
                                    prediction.start, prediction.end, prediction.strand, chosen.read_seq)
            answer_names = '; '.join(sorted(set(cp.inverse_answers[a] for a, _ in answers)))
            pred_start_codon = ''
            pred_stop_codon = ''
            for a, cod in answers:
                name = cp.inverse_answers.get(a)
                if name and 'start' in name and cod:
                    pred_start_codon = cod
                if name and 'stop' in name and cod:
                    pred_stop_codon = cod

            row = [prediction.read_name, chosen.read_start, chosen.read_end, chosen.read_strand,
                   prediction.pred_id, prediction.start, prediction.end, prediction.strand,
                   chosen.gene_id or f"cds_{chosen.cds_start}_{chosen.cds_end}",
                   chosen.cds_start, chosen.cds_end, chosen.cds_strand,
                   'generated', chosen.captures_cds_start, chosen.captures_cds_end,
                   not chosen.captures_cds_start and not chosen.captures_cds_end,
                   answer_names, pred_start_codon, pred_stop_codon]
            rows.append(row)

    write_validation(args.out, rows)

if __name__ == '__main__':
    main()
