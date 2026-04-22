#!/bin/bash
# Generic UID/HOME shim used by both the pipeline and the Jupyter image.
#
# Reads host UID/GID off the bind-mounted /project (so outputs land with host
# ownership), overridable via HOST_UID / HOST_GID env vars. Points HOME at
# /tmp so conda/mamba/jupyter have a writable cache dir. Then drops
# privileges with gosu and exec's whatever inner command the image supplies.
set -euo pipefail

HOST_UID="${HOST_UID:-$(stat -c %u /project 2>/dev/null || echo 0)}"
HOST_GID="${HOST_GID:-$(stat -c %g /project 2>/dev/null || echo 0)}"

export HOME=/tmp
export PROJECT_ROOT=/project

exec gosu "${HOST_UID}:${HOST_GID}" "$@"
