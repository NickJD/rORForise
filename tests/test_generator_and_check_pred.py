import csv
import tempfile
import unittest
from pathlib import Path

from rORForise import check_pred as cp
from rORForise import evaluate as ev
from rORForise import generate_testing_pred as gen


def answer_names(answer_set):
    return {cp.inverse_answers[answer] for answer, _ in answer_set}


class GeneratorAndEvaluationTest(unittest.TestCase):
    def test_strand_mode_controls_cds_strands(self):
        plus = gen.generate_cds_features(num_cds=10, strand_mode="plus")
        minus = gen.generate_cds_features(num_cds=10, strand_mode="minus")
        self.assertEqual({strand for _, _, strand in plus}, {"+"})
        self.assertEqual({strand for _, _, strand in minus}, {"-"})

    def test_minus_strand_correct_start_generation(self):
        cds_start, cds_end, cds_strand = 100, 199, "-"
        read_start, read_end, read_strand = 150, 220, "-"
        read_seq = gen.generate_sequence(read_end - read_start + 1)
        pred_start, pred_end, pred_strand, read_seq, _ = gen.generate_prediction(
            read_start,
            read_end,
            read_strand,
            cds_start,
            cds_end,
            cds_strand,
            read_seq,
            scenario="correct_start",
        )

        actual = cp.check_pred(
            cds_start, cds_end, cds_strand,
            read_start, read_end, read_strand,
            pred_start, pred_end, pred_strand,
            read_seq,
        )
        self.assertIn("correct start", answer_names(actual))

    def test_evaluation_percentages_do_not_exceed_100_for_simple_case(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            bed = tmp_path / "reads.bed"
            gff = tmp_path / "preds.gff"
            out = tmp_path / "out"

            read_seq = "C" * 18
            bed.write_text(
                "\t".join([
                    "chrom", "read_name", "read_start", "read_end", "read_strand",
                    "score", "feature_type", "cds_start", "cds_end", "cds_strand",
                    "read_sequence",
                ]) + "\n" +
                "\t".join(["chr1", "read1", "95", "112", "+", ".", "CDS", "100", "108", "+", read_seq]) + "\n",
                encoding="utf-8",
            )
            gff.write_text("##gff-version 3\nread1\tTest\tCDS\t9\t14\t.\t+\t0\tID=pred1\n", encoding="utf-8")

            _, preds = ev.read_preds(gff, out, "case")
            ev.evaluate(bed, preds, 0.5, out, "case", overlap_threshold=1)

            with (out / "case_detailed_results.csv").open(newline="", encoding="utf-8") as fh:
                for row in csv.DictReader(fh):
                    self.assertLessEqual(float(row["percentage"]), 100.0, row)


if __name__ == "__main__":
    unittest.main()
