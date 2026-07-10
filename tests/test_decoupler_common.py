import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))

from _decoupler_common import ulm_contribution_table  # noqa: E402


GENES = list("ABCDEFGH")
STAT = [5.0, -3.0, 2.0, 0.5, -1.5, 4.0, -0.25, 1.0]
MAT = pd.DataFrame([STAT], index=["condition"], columns=GENES)
NET = pd.DataFrame(
    {
        "source": ["T1", "T1", "T1", "T2", "T2", "T2"],
        "target": ["A", "B", "C", "D", "E", "F"],
        "weight": [1.0, -1.0, 1.0, -1.0, 1.0, -1.0],
    }
)


def weight_vector(source):
    sub = NET[NET["source"] == source].set_index("target")["weight"]
    return np.array([sub.get(g, 0.0) for g in GENES])


def pearson(source):
    return float(np.corrcoef(weight_vector(source), np.array(STAT))[0, 1])


SCORES = pd.Series({s: pearson(s) for s in ("T1", "T2")})


class UlmContributionTableTests(unittest.TestCase):
    def test_contributions_reconstruct_pearson_r(self):
        table = ulm_contribution_table(MAT, NET, SCORES)
        y = np.array(STAT)
        n = len(GENES)

        for source in ("T1", "T2"):
            x = weight_vector(source)
            # The emitted rows alone are the whole covariance numerator: non-targets
            # carry weight 0, so no residual is left behind.
            cov_numerator = table[table["source"] == source]["contribution"].sum()
            r = cov_numerator / ((n - 1) * x.std(ddof=1) * y.std(ddof=1))
            self.assertAlmostEqual(r, pearson(source), places=10)

    def test_contribution_is_weight_times_centred_stat(self):
        table = ulm_contribution_table(MAT, NET, SCORES)
        centred = pd.Series(STAT, index=GENES) - np.mean(STAT)
        expected = table["weight"].to_numpy() * centred.loc[table["target"]].to_numpy()
        np.testing.assert_allclose(table["contribution"].to_numpy(), expected)

    def test_rank_is_oriented_by_score_sign(self):
        table = ulm_contribution_table(MAT, NET, SCORES)
        for source in ("T1", "T2"):
            sub = table[table["source"] == source].sort_values("rank")
            first = sub.iloc[0]["contribution"]
            if SCORES[source] >= 0:
                self.assertAlmostEqual(first, sub["contribution"].max())
            else:
                self.assertAlmostEqual(first, sub["contribution"].min())
            self.assertEqual(list(sub["rank"]), list(range(1, len(sub) + 1)))

    def test_emits_only_scored_sources(self):
        table = ulm_contribution_table(MAT, NET, SCORES[["T1"]])
        self.assertEqual(set(table["source"]), {"T1"})

    def test_targets_absent_from_matrix_are_dropped(self):
        net = pd.concat(
            [NET, pd.DataFrame({"source": ["T1"], "target": ["ZZZ"], "weight": [1.0]})],
            ignore_index=True,
        )
        table = ulm_contribution_table(MAT, net, SCORES)
        self.assertNotIn("ZZZ", set(table["target"]))


if __name__ == "__main__":
    unittest.main()
