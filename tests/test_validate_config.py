import copy
import sys
import unittest
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))

from validate_config import load_config_tables, validate_config  # noqa: E402


FIXTURE_CONFIG = ROOT / "tests" / "test_data" / "config_test.yaml"


def load_fixture():
    with FIXTURE_CONFIG.open() as handle:
        return yaml.safe_load(handle)


class ValidateConfigTests(unittest.TestCase):
    def test_fixture_config_validates(self):
        cfg = load_fixture()
        samples, contrasts, counts_header = load_config_tables(cfg)
        validate_config(cfg, samples, contrasts, counts_header)

    def test_unknown_config_key_fails(self):
        cfg = load_fixture()
        cfg["resources"] = {"threads_default": 4}
        samples, contrasts, counts_header = load_config_tables(cfg)
        with self.assertRaisesRegex(ValueError, "Unsupported config key: resources"):
            validate_config(cfg, samples, contrasts, counts_header)

    def test_duplicate_sample_fails(self):
        cfg = load_fixture()
        samples, contrasts, counts_header = load_config_tables(cfg)
        samples = copy.deepcopy(samples)
        samples.loc[1, "sample"] = samples.loc[0, "sample"]
        with self.assertRaisesRegex(ValueError, "Duplicate sample ids"):
            validate_config(cfg, samples, contrasts, counts_header)

    def test_unsupported_contrast_factor_fails(self):
        cfg = load_fixture()
        samples, contrasts, counts_header = load_config_tables(cfg)
        contrasts = copy.deepcopy(contrasts)
        contrasts.loc[0, "factor"] = "batch"
        with self.assertRaisesRegex(ValueError, "only factor 'condition' is supported"):
            validate_config(cfg, samples, contrasts, counts_header)

    def test_invalid_contrast_id_fails(self):
        cfg = load_fixture()
        samples, contrasts, counts_header = load_config_tables(cfg)
        contrasts = copy.deepcopy(contrasts)
        contrasts.loc[0, "contrast_id"] = "Treatment vs Control"
        with self.assertRaisesRegex(ValueError, "contrast_id values must match"):
            validate_config(cfg, samples, contrasts, counts_header)


if __name__ == "__main__":
    unittest.main()
