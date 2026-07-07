import unittest

from rORForise import check_pred as cp


def names(answer_set):
    return {cp.inverse_answers[answer] for answer, _ in answer_set}


def seq_with_codons(length, inserts):
    seq = ["C"] * length
    for start_1based, codon in inserts.items():
        idx = start_1based - 1
        seq[idx:idx + 3] = list(codon)
    return "".join(seq)


class CheckPredVectorsTest(unittest.TestCase):
    def assert_answers_include(self, actual, expected):
        missing = set(expected) - names(actual)
        self.assertFalse(missing, f"Missing answers: {sorted(missing)} from {sorted(names(actual))}")

    def test_plus_read_plus_prediction_exact_boundaries(self):
        read_seq = seq_with_codons(18, {6: "ATG", 12: "TAA"})
        actual = cp.check_pred(100, 108, "+", 95, 112, "+", 6, 14, "+", read_seq)
        self.assert_answers_include(actual, ["correct direction", "correct start", "correct stop", "correct frame"])

    def test_plus_cds_reverse_read_exact_boundaries(self):
        actual = cp.check_pred(100, 108, "+", 95, 112, "-", 5, 13, "-", "C" * 18)
        self.assert_answers_include(actual, ["correct direction", "correct start", "correct stop", "correct frame"])

    def test_minus_cds_forward_read_exact_boundaries(self):
        read_seq = seq_with_codons(18, {6: "TTA", 12: "CAT"})
        actual = cp.check_pred(100, 108, "-", 95, 112, "+", 6, 14, "-", read_seq)
        self.assert_answers_include(actual, ["correct direction", "correct start", "correct stop", "correct frame"])

    def test_minus_cds_reverse_read_exact_boundaries(self):
        actual = cp.check_pred(100, 108, "-", 95, 112, "-", 5, 13, "+", "C" * 18)
        self.assert_answers_include(actual, ["correct direction", "correct start", "correct stop", "correct frame"])

    def test_wrong_prediction_strand_is_incorrect_direction(self):
        actual = cp.check_pred(100, 108, "+", 95, 112, "+", 6, 14, "-", "C" * 18)
        self.assertEqual(names(actual), {"incorrect direction"})

    def test_in_frame_noncanonical_start_is_incorrect_start(self):
        read_seq = seq_with_codons(18, {9: "CCC"})
        actual = cp.check_pred(100, 108, "+", 95, 112, "+", 9, 14, "+", read_seq)
        self.assert_answers_include(actual, ["correct direction", "correct frame", "incorrect start"])


if __name__ == "__main__":
    unittest.main()
