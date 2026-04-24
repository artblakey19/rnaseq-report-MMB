import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))

from l2s2_query import build_signature, parse_results  # noqa: E402


class NullLogger:
    def info(self, *args, **kwargs):
        pass

    def warning(self, *args, **kwargs):
        pass


class L2S2QueryTests(unittest.TestCase):
    def test_build_signature_ranks_by_padj_and_direction(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            de_path = Path(tmpdir) / "de.csv"
            pd.DataFrame(
                [
                    {"gene_name": "UP2", "log2FoldChange": 2.0, "padj": 0.02},
                    {"gene_name": "DOWN1", "log2FoldChange": -1.0, "padj": 0.001},
                    {"gene_name": "UP1", "log2FoldChange": 1.0, "padj": 0.01},
                    {"gene_name": "UP1", "log2FoldChange": 1.5, "padj": 0.015},
                    {"gene_name": "NA_PADJ", "log2FoldChange": 1.0, "padj": None},
                ]
            ).to_csv(de_path, index=False)

            up, down = build_signature(str(de_path), top_up=2, top_down=2, logger=NullLogger())

        self.assertEqual(up, ["UP1", "UP2"])
        self.assertEqual(down, ["DOWN1"])

    def test_parse_results_sorts_by_fdr_then_score(self):
        response = {
            "data": {
                "currentBackground": {
                    "pairedEnrich": {
                        "consensus": [
                            {"drug": "B", "oddsRatio": 5.0, "pvalue": 0.02, "adjPvalue": 0.05},
                            {"drug": "A", "oddsRatio": 2.0, "pvalue": 0.01, "adjPvalue": 0.01},
                            {"drug": "C", "oddsRatio": 9.0, "pvalue": 0.03, "adjPvalue": 0.05},
                        ],
                        "moas": [],
                    }
                }
            }
        }

        rows = parse_results(response, logger=NullLogger())

        self.assertEqual([row["perturbagen"] for row in rows], ["A", "C", "B"])
        self.assertEqual([row["rank"] for row in rows], [1, 2, 3])


if __name__ == "__main__":
    unittest.main()
