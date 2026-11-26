#!/usr/bin/env bash
set -eu
exec "${CLEAN_PYTHON:-python3}" -B "$(dirname "$0")/clean_stage.py" "$@"
