import unittest

from rORForise import check_pred as cp
from rORForise.adapters import BedRecord, GffPrediction, _map_prediction_to_read


def _revcomp(seq):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement[base] for base in reversed(seq))


class GenomeCoordinateMappingForReverseReadsTest(unittest.TestCase):
    """A genome-coordinate GFF prediction's strand is relative to the genome,
    but check_pred expects prediction strand relative to the read's own
    (possibly reverse-complemented) orientation. `_map_prediction_to_read`
    must flip the strand for reads aligned to the genome in reverse,
    otherwise correct predictions are reported as wrong direction and vice
    versa (see history: this used to silently invert every genome-coordinate
    prediction on a minus-strand-aligned read).
    """

    def setUp(self):
        genome_fwd = ("N" * 5) + "ATG" + ("N" * 3) + "TAA" + ("N" * 4)
        self.read_seq = _revcomp(genome_fwd)
        self.bed_record = BedRecord(
            chrom="chr", read_name="r1", read_start=95, read_end=112,
            read_strand="-", feature_type="CDS", cds_start=100, cds_end=108,
            cds_strand="+", read_seq=self.read_seq,
        )

    def _check(self, pred_strand):
        prediction = GffPrediction(
            read_name="r1", pred_id="p1", start=100, end=108, strand=pred_strand,
            seqid="chr", source="x", feature_type="CDS", score=".", phase=".",
            attributes={},
        )
        pstart, pend, pdir = _map_prediction_to_read(prediction, self.bed_record, "genome")
        answers = cp.check_pred(
            self.bed_record.cds_start, self.bed_record.cds_end, self.bed_record.cds_strand,
            self.bed_record.read_start, self.bed_record.read_end, self.bed_record.read_strand,
            pstart, pend, pdir, self.bed_record.read_seq,
        )
        return {cp.inverse_answers[a] for a, _ in answers}

    def test_genome_plus_prediction_on_reverse_read_is_correct(self):
        names = self._check("+")
        self.assertIn("correct direction", names)
        self.assertIn("correct start", names)
        self.assertIn("correct stop", names)

    def test_genome_minus_prediction_on_reverse_read_is_incorrect_direction(self):
        names = self._check("-")
        self.assertEqual(names, {"incorrect direction"})


if __name__ == "__main__":
    unittest.main()
