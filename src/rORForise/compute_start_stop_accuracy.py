import argparse
import csv
import gzip
import importlib.util
import importlib
from pathlib import Path


# Load check_pred directly from the repo's src/rORForise directory when present.
repo_root = Path(__file__).resolve().parents[2]
candidate_check = repo_root / 'src' / 'rORForise' / 'check_pred.py'

cp = None

if candidate_check.exists():
    # load check_pred.py from repo
    spec2 = importlib.util.spec_from_file_location('rORForise_check_pred_local', str(candidate_check))
    cp = importlib.util.module_from_spec(spec2)
    spec2.loader.exec_module(cp)
    print(f"Loaded check_pred from repo: {candidate_check.name}")
else:
    # fallback: try importing helper module named 'check_pred' from sys.path
    try:
        cp = importlib.import_module('check_pred')
        print('Loaded check_pred from top-level module check_pred')
    except Exception as e:
        raise ImportError('Could not import check_pred module from repo or top-level: ' + str(e))


def _open(path):
    return gzip.open(path, 'rt', encoding='utf-8') if str(path).endswith('.gz') else open(path, 'r', encoding='utf-8')


def load_bed(bed_path):
    bed_map = {}
    with _open(bed_path) as fh:
        r = csv.reader(fh, delimiter='\t')
        for row in r:
            if len(row) < 11:
                continue
            # Skip header rows where numeric fields are non-numeric (e.g., 'Start')
            try:
                read_start = int(row[2])
                read_end = int(row[3])
                cds_start = int(row[7])
                cds_end = int(row[8])
            except ValueError:
                # header or malformed line -> skip
                continue
            read_name = row[1]
            bed_map[read_name] = {
                'read_start': read_start,
                'read_end': read_end,
                'read_strand': row[4],
                'cds_start': cds_start,
                'cds_end': cds_end,
                'cds_strand': row[9],
                'read_seq': row[10]
            }
    return bed_map


def read_readlevel_gff(gff_path):
    preds = {}
    opener = (gzip.open if str(gff_path).endswith('.gz') else open)
    with opener(gff_path, 'rt', encoding='utf-8') as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            cols = line.split('\t')
            if len(cols) < 9:
                continue
            read_name = cols[0].lstrip('@')
            try:
                pred_start = int(cols[3])
                pred_end = int(cols[4])
            except ValueError:
                continue
            pred_strand = cols[6]
            attr = cols[8]
            pred_id = attr.replace('ID=', '').split(';')[0].replace('@', '')
            preds.setdefault(read_name, {})[pred_id] = (pred_start, pred_end, pred_strand)
    return preds


def analyse_one(bed_map, gff_path):
    preds = read_readlevel_gff(gff_path)
    stats = {
        'total_preds': 0,
        'reads_with_preds': len(preds),
        'preds_capture_start': 0,
        'preds_correct_start': 0,
        'preds_capture_stop': 0,
        'preds_correct_stop': 0,
        # breakdown
        'both': 0,
        'start_only': 0,
        'stop_only': 0,
        'middle_only': 0,
    }

    for read, pdict in preds.items():
        if read not in bed_map:
            continue
        be = bed_map[read]
        rs = be['read_start']; re = be['read_end']
        cds_s = be['cds_start']; cds_e = be['cds_end']
        read_dir = be['read_strand']; cds_dir = be['cds_strand']
        read_seq = be['read_seq']

        for pid, (ps, pe, pdir) in pdict.items():
            stats['total_preds'] += 1
            # determine whether read captures CDS start/stop using same logic as evaluate
            if cds_dir == '+' and read_dir == '+':
                read_captures_start = rs <= cds_s
                read_captures_stop = re >= cds_e
            elif cds_dir == '+' and read_dir == '-':
                read_captures_start = rs <= cds_s
                read_captures_stop = re >= cds_e
            elif cds_dir == '-' and read_dir == '+':
                read_captures_start = re >= cds_e
                read_captures_stop = rs <= cds_s
            elif cds_dir == '-' and read_dir == '-':
                read_captures_start = re >= cds_e
                read_captures_stop = rs <= cds_s
            else:
                read_captures_start = False
                read_captures_stop = False

            # classify into breakdown
            if read_captures_start and read_captures_stop:
                stats['both'] += 1
            elif read_captures_start:
                stats['start_only'] += 1
            elif read_captures_stop:
                stats['stop_only'] += 1
            else:
                stats['middle_only'] += 1

            # call check_pred with read_seq (it expects oriented read_seq as in BED)
            answers = cp.check_pred(cds_s, cds_e, cds_dir, rs, re, read_dir, ps, pe, pdir, read_seq)
            ans_names = {cp.inverse_answers[a] for a, _ in answers}

            if read_captures_start:
                stats['preds_capture_start'] += 1
                if 'correct start' in ans_names:
                    stats['preds_correct_start'] += 1
            if read_captures_stop:
                stats['preds_capture_stop'] += 1
                if 'correct stop' in ans_names:
                    stats['preds_correct_stop'] += 1

    return stats


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-b','--bed', required=True, help='Intersect BED file (read mappings, read_seq in col 11)')
    parser.add_argument('--pred', action='append', required=True, help='toolname=predictions.gff (read-level GFF required)')
    parser.add_argument('--out', type=str, default=None, help='Optional CSV output file to write breakdown per tool')
    args = parser.parse_args()

    bed_map = load_bed(args.bed)
    print(f"Loaded BED map with {len(bed_map)} reads")

    results = []
    for entry in args.pred:
        name, path = entry.split('=', 1)
        print(f"Processing tool {name} -> {path}")
        try:
            preds_probe = read_readlevel_gff(path)
            total_mapped = sum(len(v) for v in preds_probe.values())
            print(f"  Loaded predictions: {total_mapped} across {len(preds_probe)} reads (sample 5 reads: {list(preds_probe.keys())[:5]})")
        except Exception as e:
            print(f"  Error loading predictions for {name}: {e}")
            preds_probe = {}

        stats = analyse_one(bed_map, path)
        results.append((name, stats))

    # print summary
    for name, s in results:
        print('==', name, '==')
        print(' total mapped preds:', s['total_preds'])
        print(' reads_with_preds:', s['reads_with_preds'])
        cs = s['preds_capture_start']
        cpv = s['preds_correct_start']
        start_pct = (cpv/cs*100) if cs>0 else 0.0
        cs2 = s['preds_capture_stop']
        cpv2 = s['preds_correct_stop']
        stop_pct = (cpv2/cs2*100) if cs2>0 else 0.0
        print(' preds where read captures CDS start:', cs, f'({start_pct:.2f}% correct start)')
        print(' preds where read captures CDS stop :', cs2, f'({stop_pct:.2f}% correct stop)')
        print(' breakdown: both(start+stop) =', s['both'], f"({s['both']/s['total_preds']*100:.2f}% of preds)" )
        print('            start_only       =', s['start_only'], f"({s['start_only']/s['total_preds']*100:.2f}% of preds)" )
        print('            stop_only        =', s['stop_only'], f"({s['stop_only']/s['total_preds']*100:.2f}% of preds)" )
        print('            middle_only      =', s['middle_only'], f"({s['middle_only']/s['total_preds']*100:.2f}% of preds)" )
        print()

    # write CSV if requested
    if args.out:
        with open(args.out, 'w', newline='') as ofh:
            w = csv.writer(ofh)
            w.writerow(['tool','total_preds','reads_with_preds','preds_capture_start','preds_correct_start','preds_capture_stop','preds_correct_stop','both','start_only','stop_only','middle_only'])
            for name, s in results:
                w.writerow([name, s['total_preds'], s['reads_with_preds'], s['preds_capture_start'], s['preds_correct_start'], s['preds_capture_stop'], s['preds_correct_stop'], s['both'], s['start_only'], s['stop_only'], s['middle_only']])
        print('Wrote breakdown CSV to', args.out)
    print('compute_start_stop_accuracy.py finished')
