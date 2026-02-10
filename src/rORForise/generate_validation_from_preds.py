
import argparse
import csv
import gzip
from pathlib import Path

# robust import of check_pred
try:
    from rORForise import check_pred as cp
except Exception:
    try:
        import check_pred as cp
    except Exception:
        import importlib.util, os
        repo_root = Path(__file__).resolve().parents[1]
        candidate = repo_root / 'src' / 'rORForise' / 'check_pred.py'
        if candidate.exists():
            spec = importlib.util.spec_from_file_location('check_pred', str(candidate))
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            cp = mod
        else:
            raise


def read_bed_records(bed_path):
    opener = (gzip.open if str(bed_path).endswith('.gz') else open)
    bed = {}
    with opener(bed_path, 'rt', encoding='utf-8') as f:
        r = csv.reader(f, delimiter='\t')
        header = next(r, None)
        for row in r:
            if len(row) < 11:
                continue
            read_name = row[1]
            entry = {
                'read_start': int(row[2]),
                'read_end': int(row[3]),
                'read_strand': row[4],
                'cds_start': int(row[7]),
                'cds_end': int(row[8]),
                'cds_strand': row[9],
                'read_seq': row[10]
            }
            bed.setdefault(read_name, []).append(entry)
    return bed


def read_gff(gff_path):
    opener = (gzip.open if str(gff_path).endswith('.gz') else open)
    preds = []
    with opener(gff_path, 'rt', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.strip().split('\t')
            if len(cols) < 9:
                continue
            read_name = cols[0].lstrip('@')
            start = int(cols[3])
            end = int(cols[4])
            strand = cols[6]
            attrs = cols[8]
            pred_id = attrs.replace('ID=', '').split(';')[0]
            preds.append((read_name, pred_id, start, end, strand))
    return preds


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

    bed = read_bed_records(args.bed)
    preds = read_gff(args.gff)

    rows = []
    for read_name, pred_id, pstart, pend, pstr in preds:
        entries = bed.get(read_name, [])
        # find the bed entry whose cds range matches the pred gff's cds in attributes if possible
        chosen = None
        for e in entries:
            # choose nearest matching mapping by overlapping coordinates
            if e['read_start'] <= pstart <= e['read_end'] or e['read_start'] <= pend <= e['read_end']:
                chosen = e
                break
        if chosen is None and entries:
            chosen = entries[0]
        if chosen is None:
            continue
        read_seq = chosen['read_seq']
        cds_start = chosen['cds_start']
        cds_end = chosen['cds_end']
        cds_strand = chosen['cds_strand']
        read_start = chosen['read_start']
        read_end = chosen['read_end']
        read_strand = chosen['read_strand']

        answers = cp.check_pred(cds_start, cds_end, cds_strand,
                                read_start, read_end, read_strand,
                                pstart, pend, pstr, read_seq)
        answer_names = '; '.join(sorted(set(cp.inverse_answers[a] for a, _ in answers)))
        # codon fields: find any codons in answer tuples
        pred_start_codon = ''
        pred_stop_codon = ''
        for a, cod in answers:
            name = cp.inverse_answers.get(a)
            if name and 'start' in name and cod:
                pred_start_codon = cod
            if name and 'stop' in name and cod:
                pred_stop_codon = cod

        row = [read_name, read_start, read_end, read_strand,
               pred_id, pstart, pend, pstr,
               f"cds_{cds_start}_{cds_end}", cds_start, cds_end, cds_strand,
               'generated', False, False, False,
               answer_names, pred_start_codon, pred_stop_codon]
        rows.append(row)

    write_validation(args.out, rows)

if __name__ == '__main__':
    main()

