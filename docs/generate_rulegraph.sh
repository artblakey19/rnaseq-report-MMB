#!/usr/bin/env bash
# Regenerate docs/rulegraph.{dot,svg} from the current Snakefile.
# Requires: snakemake (project venv) + graphviz (system `dot`).
set -euo pipefail

cd "$(dirname "$0")/.."

.venv/bin/snakemake \
    --snakefile workflow/Snakefile \
    --rulegraph > docs/rulegraph.dot

dot -Tsvg docs/rulegraph.dot > docs/rulegraph.svg
echo "wrote docs/rulegraph.dot, docs/rulegraph.svg"
