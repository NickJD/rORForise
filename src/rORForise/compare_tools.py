# Aggregator to compare multiple prediction GFFs using the evaluate pipeline.
import argparse
import os
import importlib.util
from pathlib import Path
import csv

# robust adapters import: only used here for loading BED map (not for GFF->read mapping)
try:
    from rORForise import adapters
except Exception:
    try:
        import adapters
    except Exception:
        candidate = Path(__file__).resolve().parent / 'adapters.py'
        if candidate.exists():
            spec = importlib.util.spec_from_file_location('rORForise_adapters', str(candidate))
            adapters = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(adapters)
        else:
            raise

# robust evaluate module import: prefer package submodule 'rORForise.evaluate',
# but fall back to top-level evaluate.py if necessary.
import importlib
try:
    ev = importlib.import_module('rORForise.evaluate')
except Exception:
    try:
        ev = importlib.import_module('evaluate')
    except Exception:
        candidate = Path(__file__).resolve().parent / 'evaluate.py'
        if candidate.exists():
            spec = importlib.util.spec_from_file_location('rORForise_evaluate', str(candidate))
            ev = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(ev)
        else:
            raise


def run_eval_for_preds(bed, pred_gff, outdir, prefix, gc_prob):
    # read_preds returns (total_preds, preds_dict)
    total_preds, preds = ev.read_preds(pred_gff, outdir, prefix, verbose=False)

    # run evaluate (expects preds to be the read->preds dict)
    ev.evaluate(bed, preds, gc_prob, outdir, output_prefix=prefix, overlap_threshold=60, verbose=False)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('-b','--bed', required=True)
    p.add_argument('-o','--outdir', required=True)
    p.add_argument('--pred', action='append', help='toolname=predictions.gff (read-level GFF required)', required=True)
    p.add_argument('--gc_prob', type=float, default=0.32, help='GC probability for evaluation (default: 0.32)')
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    # collect per-tool data for summary and pairwise comparisons
    bed_map = adapters.load_bed_map(args.bed)
    per_tool = {}
    for pred in args.pred:
        name, path = pred.split('=', 1)
        prefix = f"cmp_{name}"
        print(f"Running evaluation for tool '{name}' -> {path}")
        # This will call ev.read_preds internally and run evaluate
        run_eval_for_preds(args.bed, path, args.outdir, prefix, args.gc_prob)

        # read the context accuracy file produced by evaluate
        metrics_file = os.path.join(args.outdir, f"{prefix}_context_accuracy.csv")
        metrics = {}
        if os.path.exists(metrics_file):
            with open(metrics_file) as fh:
                r = csv.DictReader(fh)
                metrics = {row['context']: row for row in r}
        else:
            print(f"Warning: expected metrics file not found for {name}: {metrics_file}")

        # Load predictions from read-level GFF using evaluate.read_preds
        try:
            total_loaded, preds = ev.read_preds(path, args.outdir, prefix, verbose=False)
        except Exception as e:
            print(f"Error reading predictions for {name}: {e}")
            preds = {}
            total_loaded = 0

        total_preds_mapped = sum(len(v) for v in preds.values())
        both = start_only = stop_only = middle_only = 0
        for read, pd in preds.items():
            if read not in bed_map:
                continue
            be = bed_map[read][0]
            rs = be['read_start']; re = be['read_end']
            cds_s = be['cds_start']; cds_e = be['cds_end']
            read_dir = be['read_strand']; cds_dir = be['cds_strand']
            for pid, (ps, pe, pstr) in pd.items():
                if cds_dir == '+' and read_dir == '+':
                    cap_start = rs <= cds_s
                    cap_stop = re >= cds_e
                elif cds_dir == '+' and read_dir == '-':
                    cap_start = rs <= cds_s
                    cap_stop = re >= cds_e
                elif cds_dir == '-' and read_dir == '+':
                    cap_start = re >= cds_e
                    cap_stop = rs <= cds_s
                elif cds_dir == '-' and read_dir == '-':
                    cap_start = re >= cds_e
                    cap_stop = rs <= cds_s
                else:
                    cap_start = cap_stop = False

                if cap_start and cap_stop:
                    both += 1
                elif cap_start:
                    start_only += 1
                elif cap_stop:
                    stop_only += 1
                else:
                    middle_only += 1

        per_tool[name] = {
            'metrics': metrics,
            'preds': preds,
            'total_preds_mapped': total_preds_mapped,
            'both': both,
            'start_only': start_only,
            'stop_only': stop_only,
            'middle_only': middle_only,
        }

    # write a per-tool summary
    outsum = os.path.join(args.outdir, 'compare_summary.csv')
    with open(outsum, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['tool', 'start_acc', 'stop_acc', 'middle_acc', 'total_preds_mapped', 'both', 'start_only', 'stop_only', 'middle_only'])
        for name, data in per_tool.items():
            metrics = data['metrics']
            w.writerow([
                name,
                metrics.get('start', {}).get('accuracy_pct', 'N/A'),
                metrics.get('stop',{}).get('accuracy_pct','N/A'),
                metrics.get('middle',{}).get('accuracy_pct','N/A'),
                data['total_preds_mapped'],
                data['both'],
                data['start_only'],
                data['stop_only'],
                data['middle_only']
            ])
    print('Wrote', outsum)

    # compute pairwise read/pred overlaps and exact matches
    pairfile = os.path.join(args.outdir, 'pairwise_read_overlap.csv')
    with open(pairfile, 'w', newline='') as pf:
        w = csv.writer(pf)
        w.writerow(['tool_a','tool_b','reads_with_preds_a','reads_with_preds_b','shared_reads','shared_reads_frac_a','shared_reads_frac_b','exact_pred_matches','total_preds_a','total_preds_b'])
        names = list(per_tool.keys())
        for i in range(len(names)):
            for j in range(i+1, len(names)):
                a = names[i]; b = names[j]
                preds_a = per_tool[a]['preds']
                preds_b = per_tool[b]['preds']
                reads_a = set(preds_a.keys())
                reads_b = set(preds_b.keys())
                shared_reads = reads_a & reads_b
                shared_count = len(shared_reads)
                reads_with_a = len(reads_a)
                reads_with_b = len(reads_b)
                shared_frac_a = (shared_count / reads_with_a) if reads_with_a>0 else 0
                shared_frac_b = (shared_count / reads_with_b) if reads_with_b>0 else 0

                # exact_pred_matches: count (per read) how many identical (ps,pe,str) tuples are present in both
                exact_matches = 0
                for read in shared_reads:
                    set_a = set((ps,pe,strd) for (ps,pe,strd) in preds_a[read].values())
                    set_b = set((ps,pe,strd) for (ps,pe,strd) in preds_b[read].values())
                    exact_matches += len(set_a & set_b)

                total_preds_a = per_tool[a]['total_preds_mapped']
                total_preds_b = per_tool[b]['total_preds_mapped']

                w.writerow([a,b,reads_with_a,reads_with_b,shared_count,f"{shared_frac_a:.4f}",f"{shared_frac_b:.4f}",exact_matches,total_preds_a,total_preds_b])
    print('Wrote pairwise overlaps to', pairfile)


if __name__ == '__main__':
     main()

