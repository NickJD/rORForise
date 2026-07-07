import csv
import json
import tempfile
import unittest
from pathlib import Path

from rORForise import benchmark, metrics


class BenchmarkAndReportTest(unittest.TestCase):
    def test_benchmark_writes_tables_figures_and_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bed = root / "reads.bed"
            tool_a = root / "tool_a.gff"
            tool_b = root / "tool_b.gff"
            outdir = root / "benchmark"

            read_seq = "C" * 5 + "ATG" + "C" * 3 + "TAA" + "C" * 4
            bed.write_text(
                "\t".join(
                    [
                        "chrom",
                        "read_name",
                        "read_start",
                        "read_end",
                        "read_strand",
                        "score",
                        "feature_type",
                        "cds_start",
                        "cds_end",
                        "cds_strand",
                        "read_sequence",
                    ]
                )
                + "\n"
                + "\t".join(["chr1", "read1", "95", "112", "+", ".", "CDS", "100", "108", "+", read_seq])
                + "\n",
                encoding="utf-8",
            )
            tool_a.write_text("##gff-version 3\nread1\tA\tCDS\t6\t14\t.\t+\t0\tID=exact\n", encoding="utf-8")
            tool_b.write_text("##gff-version 3\nread1\tB\tCDS\t9\t14\t.\t+\t0\tID=shifted\n", encoding="utf-8")

            result = benchmark.run_benchmark(
                bed,
                [f"Exact={tool_a}:read", f"Shifted={tool_b}:read"],
                outdir,
                gc_prob=0.5,
                overlap_threshold=1,
                report_format="html",
                max_tracks=2,
            )

            self.assertTrue((outdir / "per_tool_metrics.csv").exists())
            self.assertTrue((outdir / "per_prediction_results.csv").exists())
            self.assertTrue((outdir / "pairwise_prediction_overlap.csv").exists())
            self.assertTrue((outdir / "tool_rankings.csv").exists())
            self.assertTrue(result["report"].exists())
            self.assertTrue((outdir / "figures" / "context_accuracy.png").exists())
            self.assertTrue((outdir / "figures" / "prediction_lengths.png").exists())

            with (outdir / "per_tool_metrics.csv").open(newline="", encoding="utf-8") as fh:
                rows = {row["tool"]: row for row in csv.DictReader(fh)}
            self.assertEqual(rows["Exact"]["start_accuracy_pct"], "100.00")
            self.assertEqual(rows["Exact"]["stop_accuracy_pct"], "100.00")
            self.assertEqual(rows["Shifted"]["start_accuracy_pct"], "0.00")

            with (outdir / "tool_rankings.csv").open(newline="", encoding="utf-8") as fh:
                rankings = list(csv.DictReader(fh))
            self.assertRegex(rankings[0]["composite_score"], r"^\d+\.\d{2}$")

    def test_benchmark_rejects_unreadable_intersect_layout(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bed = root / "bed12_intersect_without_sequence.bed"
            gff = root / "predictions.gff"
            outdir = root / "benchmark"

            bed.write_text(
                "\t".join(
                    [
                        "chr1",
                        "10",
                        "160",
                        "read1",
                        "42",
                        "+",
                        "10",
                        "160",
                        "0,0,0",
                        "1",
                        "150,",
                        "0,",
                        "chr1",
                        "source",
                        "CDS",
                        "30",
                        "120",
                        ".",
                        "+",
                        "0",
                        "ID=cds1",
                        "90",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            gff.write_text("read1\tA\tCDS\t1\t90\t.\t+\t0\tID=pred1\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "No readable intersect records"):
                benchmark.run_benchmark(bed, [f"Tool={gff}:read"], outdir, report_format="none")

    def test_manifest_can_use_bed12_intersect_with_read_fasta(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bed = root / "bed12_intersect.bed"
            reads = root / "reads.fasta"
            gff = root / "predictions.gff"
            manifest = root / "benchmark.json"
            outdir = root / "benchmark"

            read_seq = "C" * 5 + "ATG" + "C" * 3 + "TAA" + "C" * 4
            reads.write_text(">read1\n" + read_seq + "\n", encoding="utf-8")
            bed.write_text(
                "\t".join(
                    [
                        "chr1",
                        "94",
                        "112",
                        "read1",
                        "42",
                        "+",
                        "94",
                        "112",
                        "0,0,0",
                        "1",
                        "18,",
                        "0,",
                        "chr1",
                        "source",
                        "CDS",
                        "100",
                        "108",
                        ".",
                        "+",
                        "0",
                        "ID=cds1",
                        "9",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            gff.write_text("read1\tA\tCDS\t6\t14\t.\t+\t0\tID=exact\n", encoding="utf-8")
            manifest.write_text(
                json.dumps(
                    {
                        "intersect": bed.name,
                        "reads": reads.name,
                        "tools": [
                            {
                                "name": "Exact",
                                "predictions": gff.name,
                                "coords": "read",
                                "type": "start_aware",
                            }
                        ],
                    }
                ),
                encoding="utf-8",
            )

            benchmark.run_benchmark(
                None,
                None,
                outdir,
                manifest_path=manifest,
                overlap_threshold=1,
                report_format="none",
            )

            with (outdir / "per_tool_metrics.csv").open(newline="", encoding="utf-8") as fh:
                row = next(csv.DictReader(fh))
            self.assertEqual(row["start_accuracy_pct"], "100.00")
            self.assertEqual(row["tool_type"], "start_aware")
            self.assertTrue((outdir / "benchmark_run_manifest.json").exists())

    def test_reverse_strand_start_offset_matches_check_pred_orientation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bed = root / "reads.bed"

            read_seq = list("C" * 150)
            read_seq[118:121] = list("CAT")
            bed.write_text(
                "\t".join(
                    [
                        "chrom",
                        "read_name",
                        "read_start",
                        "read_end",
                        "read_strand",
                        "score",
                        "feature_type",
                        "cds_start",
                        "cds_end",
                        "cds_strand",
                        "read_sequence",
                    ]
                )
                + "\n"
                + "\t".join(["chr1", "read1", "100", "249", "+", ".", "CDS", "80", "220", "-", "".join(read_seq)])
                + "\n",
                encoding="utf-8",
            )

            rows, summary = metrics.evaluate_tool("ReverseTool", bed, {"read1": {"pred1": (2, 121, "-")}}, overlap_threshold=1)

            self.assertEqual(rows[0]["correct_start"], 1)
            self.assertEqual(rows[0]["start_offset"], 0)
            self.assertEqual(summary["correct_start"], 1)


if __name__ == "__main__":
    unittest.main()
