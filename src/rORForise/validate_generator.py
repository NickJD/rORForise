"""Validate generated test data against check_pred logic.

This script loads the intersect BED, predictions GFF and validation CSV produced by
`generate_testing_pred.py` and for each validation row re-runs the
`check_pred` function to confirm the expected answers were predicted.

Exits with code 0 if all validation rows pass, non-zero otherwise.
"""

import argparse
import csv
import gzip
import sys
from pathlib import Path

# try to import package modules
try:
    from rORForise import check_pred as cp
except Exception:
    try:
        import check_pred as cp
    except Exception:
        # Try to load from src path using Path resolution
        import importlib.util, os
        repo_root = Path(__file__).resolve().parents[1]  # project root (rORForise/)
        candidate = repo_root / 'src' / 'rORForise' / 'check_pred.py'
        if candidate.exists():
            spec = importlib.util.spec_from_file_location('check_pred', str(candidate))
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            cp = mod
        else:
            raise ModuleNotFoundError('Could not locate check_pred module')


def read_bed(bed_path):
    opener = (gzip.open if str(bed_path).endswith('.gz') else open)
    bed_map = {}
    with opener(bed_path, 'rt', encoding='utf-8') as f:
        r = csv.reader(f, delimiter='\t')
        try:
            header = next(r)
        except StopIteration:
            return bed_map
        for row in r:
            if len(row) < 11:
                continue
            chrom, read_name, read_start, read_end, read_strand = row[0], row[1], int(row[2]), int(row[3]), row[4]
            cds_start, cds_end, cds_strand = int(row[7]), int(row[8]), row[9]
            read_seq = row[10]
            bed_map.setdefault(read_name, []).append({
                'chrom': chrom,
                'read_start': read_start,
                'read_end': read_end,
                'read_strand': read_strand,
                'cds_start': cds_start,
                'cds_end': cds_end,
                'cds_strand': cds_strand,
                'read_seq': read_seq
            })
    return bed_map


def main():
    p = argparse.ArgumentParser()
    p.add_argument('-b', '--bed', required=True)
    p.add_argument('-v', '--validation', required=True)
    args = p.parse_args()

    bed_map = read_bed(args.bed)

    failures = []
    total = 0

    with open(args.validation, 'r', encoding='utf-8') as f:
        rdr = csv.DictReader(f)
        for row in rdr:
            total += 1
            read_name = row['read_name']
            pred_start = int(row['pred_start'])
            pred_end = int(row['pred_end'])
            pred_strand = row['pred_strand']
            cds_start = int(row['cds_start'])
            cds_end = int(row['cds_end'])
            cds_strand = row['cds_strand']
            expected_answers = [s.strip() for s in row['expected_answers'].split(';') if s.strip()]

            # find matching bed row — prefer matching cds coordinates
            candidates = bed_map.get(read_name, [])
            match = None
            for c in candidates:
                if c['cds_start'] == cds_start and c['cds_end'] == cds_end:
                    match = c
                    break
            if match is None and candidates:
                match = candidates[0]

            if match is None:
                failures.append((read_name, 'no_bed_entry'))
                continue

            read_start = match['read_start']
            read_end = match['read_end']
            read_strand = match['read_strand']
            read_seq = match['read_seq']

            # run check_pred
            answers = cp.check_pred(cds_start, cds_end, cds_strand,
                                     read_start, read_end, read_strand,
                                     pred_start, pred_end, pred_strand,
                                     read_seq)
            actual = set(cp.inverse_answers.get(code) for code, _ in answers)

            # expected must be subset of actual
            missing = [e for e in expected_answers if e not in actual]
            if missing:
                failures.append((read_name, pred_start, pred_end, expected_answers, sorted(actual)))

    print(f"Validated {total} rows; failures: {len(failures)}")
    if failures:
        print('\nExample failures:')
        for ex in failures[:10]:
            print(ex)
        sys.exit(2)
    else:
        print('All validation rows matched check_pred outputs')
        sys.exit(0)

if __name__ == '__main__':
    main()

