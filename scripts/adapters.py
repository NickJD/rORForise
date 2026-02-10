"""Adapters for converting prediction outputs into a canonical preds dict

Functions:
- gff_genome_to_read_preds(gff_path, intersect_bed_path)
    Map predictions given in genome coordinates (seqname is read_name) into
    read-relative preds usable by evaluate.evaluate.

The adapters are intentionally minimal: they assume the GFF's seqname column
identifies the read name as used in the intersect BED produced by the
generator. They will return a preds dict: {read_name: {pred_id: (ps,pe,dir)}}
where ps/pe are 1-based positions on the read sequence as stored in the BED.
"""

import csv
import gzip


def _open(path):
    return gzip.open(path, 'rt', encoding='utf-8') if str(path).endswith('.gz') else open(path, 'r', encoding='utf-8')


def load_bed_map(intersect_bed_path):
    """Return a map of read_name -> list of bed entries (using the first when needed)
    bed entry columns expected: chrom, read_name, read_start, read_end, read_strand, score, feature_type, cds_start, cds_end, cds_strand, read_sequence
    """
    bed_map = {}
    with _open(intersect_bed_path) as fh:
        r = csv.reader(fh, delimiter='\t')
        header = next(r, None)
        for row in r:
            if len(row) < 11:
                continue
            read_name = row[1]
            bed_map.setdefault(read_name, []).append({
                'read_start': int(row[2]),
                'read_end': int(row[3]),
                'read_strand': row[4],
                'cds_start': int(row[7]),
                'cds_end': int(row[8]),
                'cds_strand': row[9],
                'read_seq': row[10]
            })
    return bed_map


def gff_genome_to_read_preds(gff_path, intersect_bed_path):
    """Convert GFF with either genome coordinates or read-relative coordinates to preds dict.

    For each GFF line, we expect columns like: seqname, source, feature, start, end, score, strand, frame, attributes
    seqname must be the read_name found in the intersect BED. If multiple BED rows exist for a read,
    the first is used.

    The adapter auto-detects whether the coordinates in the GFF are genome coordinates (absolute)
    or read-relative coordinates (1-based positions within the read). This is done per-record by
    comparing the GFF start/end to the read mapping in the BED.

    Returns: preds dict {read_name: {pred_id: (pstart_on_read, pend_on_read, pred_strand)}}
    """
    bed_map = load_bed_map(intersect_bed_path)
    preds = {}

    with _open(gff_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            cols = line.strip().split('\t')
            if len(cols) < 9:
                continue
            seqname = cols[0].lstrip('@')
            try:
                gstart = int(cols[3])
                gend = int(cols[4])
            except ValueError:
                # skip malformed numeric fields
                continue
            gstrand = cols[6]
            attrs = cols[8]
            pred_id = attrs.replace('ID=', '').split(';')[0].replace('@', '')

            if seqname not in bed_map:
                # cannot map this prediction to any read mapping
                continue
            bed = bed_map[seqname][0]
            read_start = bed['read_start']
            read_end = bed['read_end']
            read_strand = bed['read_strand']
            read_len = read_end - read_start + 1

            # Heuristic detection of coordinate type:
            # - If gstart/gend fall within genome read interval [read_start, read_end], assume genome coords
            # - Else if gstart/gend fall within [1, read_len], assume read-relative coords
            # - Else if gend <= read_len * 1.5, prefer read-relative
            # - Else assume genome coords
            coord_type = 'genome'
            if (gstart >= read_start and gstart <= read_end) or (gend >= read_start and gend <= read_end):
                coord_type = 'genome'
            elif (gstart >= 1 and gstart <= read_len) and (gend >= 1 and gend <= read_len):
                coord_type = 'read'
            else:
                if gend <= int(read_len * 1.5):
                    coord_type = 'read'
                else:
                    coord_type = 'genome'

            if coord_type == 'genome':
                # map genome coords to read-relative positions
                if read_strand == '+':
                    pstart = gstart - read_start + 1
                    pend = gend - read_start + 1
                else:
                    # read seq first base corresponds to genome read_end
                    pstart = read_end - gend + 1
                    pend = read_end - gstart + 1
            else:
                # GFF coordinates are already read-relative
                pstart = gstart
                pend = gend

            # clamp to read boundaries
            pstart = max(1, pstart)
            pend = max(1, min(read_len, pend))

            if seqname not in preds:
                preds[seqname] = {}
            preds[seqname][pred_id] = (pstart, pend, gstrand)

    return preds

