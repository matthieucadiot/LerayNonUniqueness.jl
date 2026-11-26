# Check the published results or rerun the certification

## Check the results on a laptop

Download `certification-evidence.zip` and the release’s plain-text `SHA256SUMS` into the same directory. Add `Proof_Float64_Clean_2611.zip` if you also want the fixed splines or the numerical tests. Choose the checksum command for your system. Each downloaded release ZIP should appear with `OK`.

```sh
sha256sum --ignore-missing -c SHA256SUMS           # Linux (GNU coreutils)
# shasum -a 256 --ignore-missing -c SHA256SUMS     # macOS
```

Then extract the evidence and run the checker with Python 3.11 or newer. It uses only the standard library.

```sh
unzip -q certification-evidence.zip
python3 certification-evidence/verify_evidence.py
```

The accepted run has ID `2d9cec66690641fa88a3aac70215c5cb`. The checker verifies the original final receipt, the linked aggregate files, and the recorded inequalities. See [Tables 1–5](TABLES.md) for the correspondence with the paper. For the rest of this guide, extract the clean proof beside the evidence and enter its directory.

```sh
unzip -q Proof_Float64_Clean_2611.zip
cd Proof_Float64_Clean_2611
python3 ../certification-evidence/verify_evidence.py --proof .
```

The clean ZIP includes a release guide in `START-HERE.md`. It maps the public archive name to the original Della instructions included with the proof.

## Run a numerical check

With Julia **1.12.3** available as `julia`, this recomputes the seven closing estimates using exact rational arithmetic and compares their values and targets to the final receipt. It runs offline using Julia’s standard library.

```sh
julia --startup-file=no -e '
    using TOML
    include("scripts/clean_closure.jl")
    values = closure()
    @assert length(values) == 7
    @assert values == TOML.parsefile(ARGS[1])["closure"]
    println("Seven exact closing estimates match the accepted run.")
' ../certification-evidence/verification.toml
```

For the interval calculations, install the three **pinned** Julia environments once with registry and package access. On a cluster choose a writable depot visible inside compute jobs and keep the same depot for the run.

```sh
export JULIA_BIN="$(command -v julia)"
julia --startup-file=no -e '@assert VERSION == v"1.12.3"'
(
    set -e
    for project in eigenfunctions residual coercivity; do
        julia --startup-file=no --project="$project" -e 'using Pkg; Pkg.instantiate()'
    done
)
```

The following small offline test evaluates actual odd and even response cells, exercises the refined integration grid, compares threaded and serial interval endpoints, and saves and reloads a temporary column checkpoint.

```sh
JULIA_PKG_OFFLINE=true julia --startup-file=no --project=coercivity \
    --threads=4 coercivity/table2/clean_test.jl
```

The expected test summary is **38 passed**. Twelve threaded repetitions must match the serial interval endpoints.

## Recompute the whole certification on a cluster

The accepted calculation used up to 12 simultaneous workers, each allocated 32 cores and 48 GiB in two-hour jobs. Checkpoints and results require a shared, writable scratch directory. Keep the extracted proof directory intact. `bash run.sh plan` prints the original pass schedule offline. The bundled `run.sh submit` launcher is specific to Princeton Della and names the authors’ account, partition and QoS.

For another scheduler, its jobs can call the same numerical entry point, `scripts/clean_stage.sh`. Use Bash, Python 3.11 or newer, the verified source tree and Julia 1.12.3. Set up these paths once from the extracted proof directory, using your site's shared scratch path.

```sh
export CLEAN_SOURCE_ROOT="$(pwd -P)"
export CLEAN_RUN_DIR="$(mktemp -d "${SCRATCH:?set SCRATCH to a shared scratch directory}/leray.XXXXXXXX")"
python3 -c 'import uuid; print(uuid.uuid4().hex)' > "$CLEAN_RUN_DIR/run-id.txt"
```

In each allocated job, export the same absolute `CLEAN_SOURCE_ROOT` and `CLEAN_RUN_DIR`, the same shared Julia depot, and the absolute `JULIA_BIN` path. The job adapter sets `stage`, `lane` and `lanes` according to the table below. With 32 cores allocated to the job, run the worker command.

```sh
export CLEAN_RUN_ID="$(cat "$CLEAN_RUN_DIR/run-id.txt")"
export SLURM_CPUS_PER_TASK=32
bash "$CLEAN_SOURCE_ROOT/scripts/clean_stage.sh" \
    "$stage" "$CLEAN_RUN_DIR" "$lane" "$lanes" 6600
```

| In order | Stage | Assignments | What completes |
| --- | --- | --- | --- |
| 1 | `base` | `lanes=12`, one job per `lane=1` through `12` | Eigenmodes, boxes, residual columns |
| 2 | `assemble_base` | `lane=1`, `lanes=1` | Upstream certificates |
| 3 | `main` | `lanes=12`, one job per `lane=1` through `12` | Coefficient bounds, 15 matrix cases, 58 responses |
| 4 | `final` | `lane=1`, `lanes=1` | Numerical acceptance and final receipt |

Wait for every lane of each stage to exit successfully before advancing. Run at most one copy of a given stage and lane at once. A worker reaching its checkpoint budget can exit 75. Schedule that same lane again with its original identity, lane count and output directory. Investigate other nonzero exits in the job log. The final stage writes `$CLEAN_RUN_DIR/verification.toml` with `all_checks_passed = true` after acceptance. A fresh run receives a new identity, timestamps and file digests.

The original Princeton scheduler wrapper also checks the source and input manifest at submission and at final verification. For another scheduler, run `sha256sum -c CLEAN_SHA256SUMS` (or `shasum -a 256 -c CLEAN_SHA256SUMS` on macOS) from the extracted source directory before submission and again before final assembly.
