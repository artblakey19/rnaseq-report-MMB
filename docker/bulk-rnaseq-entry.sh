#!/bin/bash
# Pipeline image sub-command dispatch.
#
#   init   → runs workflow/scripts/init_project.py interactively, writing
#            config/config.yaml, config/samples.tsv, config/contrasts.tsv
#            inside /project (cwd) based on host-mounted counts files.
#   *      → default: runs the Snakemake pipeline. Positional args are
#            forwarded to snakemake (e.g. -c1, --configfile ..., --dry-run).
set -euo pipefail

case "${1:-}" in
    init)
        shift
        cd /project
        exec micromamba run -n base \
            python /app/workflow/scripts/init_project.py "$@"
        ;;
    *)
        exec micromamba run -n base \
            snakemake --snakefile /app/workflow/Snakefile --use-conda "$@"
        ;;
esac
