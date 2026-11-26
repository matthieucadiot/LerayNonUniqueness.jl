"""Run one disjoint lane of a prebooked Slurm pass. Numerical programs validate the results at assembly."""
import json
import os
import subprocess
import sys
import time
from pathlib import Path

SOURCE = Path(__file__).resolve().parents[1]


def span(total, lane, count):
    return (1 + (lane - 1) * total // count, lane * total // count)


def stage_main(stage, run, lane, lanes, seconds):
    run = run.resolve(strict=True)
    run_id = os.environ["CLEAN_RUN_ID"]
    if run != Path(os.environ["CLEAN_RUN_DIR"]).resolve() or not 1 <= lane <= lanes:
        raise ValueError("Unexpected run directory or lane")
    end = time.monotonic() + seconds
    threads = int(os.environ.get("SLURM_CPUS_PER_TASK", "32"))
    if threads < 1 or seconds < 1:
        raise ValueError("Expected positive thread count and time budget")
    results = run / "results"
    marks = run / ".science" / stage / str(lane)
    marks.mkdir(parents=True, exist_ok=True)
    base_env = dict(os.environ, PROOF_DATA=str(SOURCE / "data"), COERCIVITY_RESULTS=str(results / "coercivity"),
                    COERCIVITY_DATA=str(SOURCE / "data/Coercivity"), JULIA_LOAD_PATH="@:@stdlib",
                    EIGENFUNCTIONS_CERTIFICATE=str(results / "eigenfunctions/Eigenfunctions.jld2"),
                    RESIDUAL_CERTIFICATE=str(results / "residual/Residual_Table1.jld2"),
                    OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")

    def execute(key, component, script, args=(), limit=None):
        stamp = marks / (key + ".json")
        stamp_data = {"run_id": run_id, "stage": stage, "lane": lane, "lanes": lanes,
                      "component": component, "script": script, "args": list(args)}
        if stamp.is_file():
            if json.loads(stamp.read_text()) != stamp_data:
                raise RuntimeError(f"Unexpected progress marker {stamp}")
            return
        remaining = int(end - time.monotonic())
        # Leave time to flush an atomic prefix before Slurm's hard deadline.
        if remaining < 30:
            raise SystemExit(75)
        budget = min(remaining, limit) if limit else remaining
        env = dict(base_env, PROOF_OUTPUT=str(results / component),
                   PROOF_JOB_SECONDS=str(budget), TABLE2_JOB_SECONDS=str(budget))
        folder = SOURCE / component
        project = SOURCE / ("coercivity" if component == "scripts" else component)
        call = [os.environ["JULIA_BIN"], "--startup-file=no", f"--threads={threads},0",
                f"--project={project}", str(folder / script), *map(str,args)]
        if key == "responses":
            call.append(str(budget))
        print(f"{stage} lane {lane}/{lanes}: {key} ({budget}s soft budget)", flush=True)
        result = subprocess.run(call, cwd=folder, env=env)
        if result.returncode:
            raise SystemExit(result.returncode if result.returncode > 0 else 128 - result.returncode)
        tmp = stamp.with_suffix(".tmp")
        tmp.write_text(json.dumps(stamp_data,sort_keys=True))
        tmp.replace(stamp)

    if stage == "base":
        execute("modes","eigenfunctions","Run_modes.jl",span(207,lane,lanes))
        execute("boxes","coercivity","proposition/Run_boxes.jl",span(2482,lane,lanes))
        execute("residual","residual","Run_columns.jl",span(804,lane,lanes))
    elif stage == "assemble_base":
        if (lane,lanes) != (1,1): raise ValueError("Assembly requires one lane")
        execute("modes","eigenfunctions","Assemble.jl")
        execute("residual","residual","Assemble.jl")
        execute("boxes","coercivity","proposition/Assemble_boxes.jl")
    elif stage == "main":
        execute("farfield","coercivity","proposition/Run_farfield_columns.jl",span(804,lane,lanes))
        execute("bounds","coercivity","proposition/Run_bounds_columns.jl",span(2521,lane,lanes))
        for k in range(15):
            row = 1+(lane-1+k)%15
            execute("matrix_"+str(row),"coercivity","proposition/clean_matrices.jl",("work",row,lane,lanes))
        execute("responses","coercivity","table2/clean_response.jl",(lane,lanes))
    elif stage == "final":
        if (lane,lanes) != (1,1): raise ValueError("Final assembly requires one lane")
        execute("bounds","coercivity","proposition/Assemble_bounds.jl")
        execute("table2","coercivity","table2/clean_assemble.jl")
        execute("table4","coercivity","proposition/clean_matrices.jl",("assemble",))
        execute("verification","scripts","clean_final.jl")
    else:
        raise ValueError(f"Unknown stage {stage}")


if __name__ == "__main__":
    if len(sys.argv) != 6:
        raise SystemExit("usage: clean_stage.sh STAGE RUN_ROOT LANE LANES SOFT_SECONDS")
    stage_main(sys.argv[1], Path(sys.argv[2]), *map(int,sys.argv[3:]))
