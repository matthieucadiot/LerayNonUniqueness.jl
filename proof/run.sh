#!/usr/bin/env bash
# Submit or inspect a clean Della certification.
set -euo pipefail
root=$(cd "$(dirname "$0")" && pwd -P)
action=${1:-submit}
run=${2:-"$(dirname "$root")/clean_run_20261123"}
if [ "$#" -gt 2 ]; then
    echo 'Usage: bash run.sh [plan|submit|status] [absolute scratch run directory]' >&2
    exit 2
fi
py=
for candidate in python3 python3.13 python3.12 python3.11; do
    if command -v "$candidate" >/dev/null && "$candidate" -c 'import sys; sys.exit(sys.version_info < (3, 11))' 2>/dev/null; then
        py=$(command -v "$candidate")
        break
    fi
done
if [ -z "$py" ]; then
    echo 'Activate the Della Python/Conda environment with Python 3.11 or newer and repeat this command.' >&2
    exit 1
fi
export CLEAN_PYTHON="$py" PYTHONDONTWRITEBYTECODE=1
case "$action" in
    plan) exec "$py" -B "$root/scripts/weekend.py" plan ;;
    status) exec "$py" -B "$root/scripts/weekend.py" status "$run" ;;
    submit) ;;
    *) echo 'Usage: bash run.sh [plan|submit|status] [absolute scratch run directory]' >&2; exit 2 ;;
esac
old_runtime="$HOME/Leray_NoSwirl_Proof/Certification_2/Proof_Float64_20261107/.runtime"
runtime="$(dirname "$root")/clean_runtime_20261120"
if [ -z "${JULIA_BIN:-}" ] || [ -z "${JULIA_DEPOT_PATH:-}" ]; then
    (
        flock -x 9
        if [ ! -f "$runtime/READY" ]; then
            mkdir -p "$runtime"
            rsync -a "$old_runtime/julia-1.12.3" "$old_runtime/depot" "$runtime/"
            test -x "$runtime/julia-1.12.3/bin/julia"
            test -d "$runtime/depot/packages"
            touch "$runtime/READY"
        fi
    ) 9>"$runtime.lock"
fi
export JULIA_BIN=${JULIA_BIN:-"$runtime/julia-1.12.3/bin/julia"}
export JULIA_DEPOT_PATH=${JULIA_DEPOT_PATH:-"$runtime/depot:"}
export JULIA_LOAD_PATH='@:@stdlib'
if [ ! -f "$run/.weekend/run.json" ]; then
    for group in eigenfunctions residual coercivity; do
        "$JULIA_BIN" --startup-file=no --threads=1 --project="$root/$group" -e '
            VERSION == v"1.12.3" || error("Use Julia 1.12.3")
            using TOML
            for pkg in keys(TOML.parsefile(Base.active_project())["deps"])
                isnothing(Base.find_package(pkg)) && error("Missing package $pkg in $(Base.active_project()). Check JULIA_DEPOT_PATH.")
            end'
    done
fi
exec "$py" -B "$root/scripts/weekend.py" submit "$run"
