#!/usr/bin/env python3
"""Prebook a finite Della run. All paths under a run are private to that run."""
import argparse
from collections import Counter
from contextlib import contextmanager
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import signal
import stat
import subprocess
import sys
import tempfile
import tomllib
import uuid

ROOT = Path(__file__).resolve().parents[1]
# Prior 32-core response benchmark: 490–550 worker-hours for an unrefined sweep.
# The 4,289-column repair pilot budget was another 90–130 worker-hours. Fresh
# matrices, lane imbalance, transient failures and setup require additional room.
# 40 main passes allow 880 productive worker-hours at the 110-minute soft budget.
STAGES = (("base", 12, 6), ("assemble_base", 1, 2), ("main", 12, 40), ("final", 1, 2))
SOFT = 6600
CPUS = 32
MEMORY = "48G"
MANIFEST = "CLEAN_SHA256SUMS"
STAGE_PROGRAM = "scripts/clean_stage.sh"
WRAPPER = "scripts/weekend_job.sh"

def require(value, message):
    if not value:
        raise RuntimeError(message)

def now():
    return datetime.now(timezone.utc).isoformat()

def atomic_text(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, delete=False) as f:
        temporary = Path(f.name)
        try:
            f.write(content)
            f.flush()
            os.fsync(f.fileno())
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
    os.replace(temporary, path)

def write_json(path, value):
    atomic_text(path, json.dumps(value, sort_keys=True, indent=2) + "\n")

def read(path):
    return json.loads(path.read_text())

@contextmanager
def lock(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        yield

def digest(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()

def manifest(root, expected=None, full=False):
    listing = root / MANIFEST
    actual = digest(listing)
    if expected is not None:
        require(actual == expected, "Clean package manifest changed after submission")
    seen = set()
    for line in listing.read_text().splitlines():
        match = re.fullmatch(r"([a-f0-9]{64}) [ *](.+)", line)
        require(match, f"Malformed line in {MANIFEST}: {line[:120]}")
        sha, name = match.groups()
        relative = Path(name)
        require(not relative.is_absolute() and ".." not in relative.parts and relative not in seen,
                f"Unsafe or duplicate manifest entry: {name}")
        file = root / relative
        require(file.is_relative_to(root) and file.resolve() == file and file.is_file(),
                f"Missing or indirect input: {name}")
        seen.add(relative)
        if full or file.suffix in (".jl", ".sh", ".py", ".toml", ".tsv"):
            require(digest(file) == sha, f"Changed input: {name}")
    for name in (STAGE_PROGRAM, WRAPPER, "scripts/weekend.py"):
        require(Path(name) in seen, f"{name} missing from {MANIFEST}")
    return actual

def plan():
    jobs = []
    previous_stage = []
    for stage, lanes, passes in STAGES:
        current = []
        for p in range(1, passes + 1):
            key = f"{stage}.{p:02d}"
            dependency = list(previous_stage) if p == 1 else [current[-1]]
            jobs.append({"key": key, "stage": stage, "pass": p, "lanes": lanes, "parents": dependency})
            current.append(key)
        previous_stage = current
    return jobs

def command(args):
    result = subprocess.run(args, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode:
        raise RuntimeError(f"{shlex.join(map(str, args))}: {result.stderr.strip() or result.stdout.strip()}")
    return result.stdout.strip()

def metadata(run):
    return run / ".weekend"

def owner(directory):
    require(directory.is_dir() and not directory.is_symlink() and directory.stat().st_uid == os.getuid(),
            f"Run directory is missing or belongs to another user: {directory}")
    require(not (stat.S_IMODE(directory.stat().st_mode) & 0o022),
            f"Run directory must not be group/world writable: {directory}")

def stored(run):
    owner(run)
    f = metadata(run) / "run.json"
    state = read(f)
    require(state["root"] == str(ROOT) and state["run"] == str(run), "Run paths differ from submission")
    return state

def software():
    julia = shutil.which(os.environ.get("JULIA_BIN") or "julia")
    require(julia and Path(julia).is_file(), "Set JULIA_BIN to the Julia 1.12.3 executable")
    require(command([julia, "--version"]) == "julia version 1.12.3", "Julia 1.12.3 is required")
    depot_env = os.environ.get("JULIA_DEPOT_PATH") or str(ROOT / ".runtime/depot")
    # A trailing separator adds Julia's bundled system depots. Never rely
    # on an implicit user depot as the first writable cache.
    depots = depot_env.split(os.pathsep)
    if depots[-1] == "":
        depots.pop()
    require(depots, "Set an explicit primary Julia depot directory")
    for p in depots:
        require(p and Path(p).is_absolute() and Path(p).is_dir(),
                "Set JULIA_DEPOT_PATH to existing absolute directories visible from all compute nodes")
    return str(Path(julia).resolve()), depot_env

def init_run(run):
    require(not (run == ROOT or run.is_relative_to(ROOT) or ROOT.is_relative_to(run)),
            "Choose a scratch run directory outside the clean package")
    if not run.exists():
        run.mkdir(mode=0o700, parents=True, exist_ok=True)
    owner(run)
    m = metadata(run)
    require(m.is_dir() or not any(run.iterdir()), "Choose a NEW or empty run directory")
    m.mkdir(mode=0o700, exist_ok=True)
    owner(m)

def make_state(run):
    require(set(run.iterdir()) == {metadata(run)}, "An uninitialized run contains existing scientific files")
    package = manifest(ROOT, full=True)
    julia, depot = software()
    tag = uuid.uuid4().hex
    rows = plan()
    for index, row in enumerate(rows, start=1):
        row["name"] = f"lc_{tag[:12]}_{index:02d}"
    state = {"schema": 1, "root": str(ROOT), "run": str(run), "run_id": tag,
             "manifest": package, "julia": julia, "depot": depot, "jobs": rows,
             "created": now(), "committed": False, "released": False}
    (run / "logs").mkdir(exist_ok=True)
    write_json(metadata(run) / "run.json", state)
    return state

def queued():
    rows = command(["squeue", "-h", "--me", "-o", "%F|%j|%T"]).splitlines()
    by_name = {}
    for line in rows:
        job, name, state = line.split("|")
        if name.startswith("lc_"):
            by_name.setdefault(name, set()).add((job, state))
    return by_name

def job_id(snapshot, name):
    # A lost sbatch acknowledgement is recoverable by the unique held job name.
    ids = {job for job, _ in snapshot.get(name, ())}
    require(len(ids) <= 1 and all(i.isdecimal() for i in ids), f"Ambiguous existing job {name}")
    return next(iter(ids), None)

def submit(run):
    init_run(run)
    m = metadata(run)
    with lock(m / "submit.lock"):
        state = stored(run) if (m / "run.json").exists() else make_state(run)
        rows = state["jobs"]
        if state["released"]:
            print((run / "BOOST_REQUEST.txt").read_text(), end="")
            return
        if not state["committed"]:
            # Only the first job needs to be held. All other jobs are behind it
            # until every new job ID is durably committed.
            known = {}
            existing = queued()
            for row in rows:
                key = row["key"]
                if "id" not in row:
                    jobid = job_id(existing, row["name"])
                    if not jobid:
                        args = ["sbatch", "--parsable", "--partition=cpu", "--qos=short",
                                "--account=aionescu", "--constraint=amd", "--nodes=1", "--ntasks=1",
                                f"--cpus-per-task={CPUS}", f"--mem={MEMORY}", "--time=02:00:00",
                                "--no-requeue", "--kill-on-invalid-dep=no",
                                f"--job-name={row['name']}", f"--chdir={run}",
                                f"--output={run}/logs/{key}_%A_%a.out", "--open-mode=append"]
                        if row["lanes"] > 1:
                            args.append(f"--array=1-{row['lanes']}%{row['lanes']}")
                        if row["parents"]:
                            args.append("--dependency=afterany:" + ":".join(known[k] for k in row["parents"]))
                        else:
                            args.append("--hold")
                        args.extend([str(ROOT / WRAPPER), sys.executable, str(ROOT / "scripts/weekend.py"),
                                     str(run), key])
                        jobid = command(args).split(";", 1)[0]
                    require(jobid.isdecimal() and jobid not in known.values(), f"Invalid Slurm ID for {key}: {jobid}")
                    row["id"] = jobid
                    write_json(m / "run.json", state)
                known[key] = row["id"]
            ids = [r["id"] for r in rows]
            lines = ["Subject: Weekend Della priority request for aionescu", "",
                     "Hi Prentice,", "",
                     "Thank you for helping us finish the certification. We are now running it "
                     "from the fixed inputs in a clean directory. Could you raise the priority "
                     "of the following Della job and array IDs for the weekend?", ""]
            lines.extend(", ".join(ids[i:i+10]) for i in range(0, len(ids), 10))
            lines.extend(["", "All passes, retries and final checks have been submitted. "
                          f"At most 12 workers run together, each requesting {CPUS} AMD CPU cores, "
                          f"{MEMORY} of memory and two hours on the short QoS "
                          f"({12*CPUS} cores at peak). Completed work is checkpointed across passes.",
                          "", "Thank you again,", "Alex"])
            atomic_text(run / "BOOST_REQUEST.txt", "\n".join(lines) + "\n")
            inventory = ["stage_pass\tjob_id\tarray_tasks\tworker_cores\tworker_memory\twall_limit"]
            inventory.extend(f"{r['key']}\t{r['id']}\t{r['lanes']}\t{CPUS}\t{MEMORY}\t02:00:00" for r in rows)
            atomic_text(run / "job_ids.tsv", "\n".join(inventory) + "\n")
            state["committed"] = True
            write_json(m / "run.json", state)
        first = rows[0]["id"]
        info = subprocess.run(["scontrol", "-o", "show", "job", first], text=True, capture_output=True)
        if info.returncode == 0 and re.search(r"\bReason=JobHeld(?:User|Admin)\b", info.stdout):
            command(["scontrol", "release", first])
        elif info.returncode:
            # If a prior release succeeded and the process was interrupted, an
            # accounting record must confirm the first job actually started.
            records = command(["sacct", "-n", "-X", "-j", first, "-P", "-o", "JobID,Start"])
            require(any(re.search(r"\|\d{4}-", line) for line in records.splitlines()),
                    "Cannot confirm release of first job. Inspect the recorded ID before retrying submit.")
        else:
            # Never report a cancelled or permanently held chain as released.
            require("JobHeld" not in info.stdout and "JobState=CANCELLED" not in info.stdout,
                    "First job is held or cancelled; inspect it before proceeding")
        state["released"] = True
        write_json(m / "run.json", state)
        print((run / "BOOST_REQUEST.txt").read_text(), end="")

def marker(run, stage, lane):
    return metadata(run) / "done" / f"{stage}.{lane}.json"

def complete(run, state, stage, lane):
    path = marker(run, stage, lane)
    if not path.is_file():
        return False
    v = read(path)
    return v.get("run_id") == state["run_id"] and v.get("manifest") == state["manifest"] and v.get("stage") == stage and v.get("lane") == lane

def final_receipt(run, state):
    with (run / "verification.toml").open("rb") as f:
        receipt = tomllib.load(f)
    require(receipt.get("all_checks_passed") is True and receipt.get("run_id") == state["run_id"],
            "Final scientific verification does not certify this run identity")
    return digest(run / "verification.toml")

def prune(run, state, row):
    stage = row["stage"]
    if not all(complete(run, state, stage, n) for n in range(1, row["lanes"] + 1)):
        return
    stamp = metadata(run) / f"pruned.{stage}.json"
    with lock(metadata(run) / f"prune.{stage}.lock"):
        if stamp.exists():
            return
        later = [r for r in state["jobs"] if r["stage"] == stage and r["pass"] > row["pass"]]
        cancelled = []
        if later:
            try:
                snapshot = queued()
                for candidate in later:
                    entries = snapshot.get(candidate["name"], set())
                    if entries == {(candidate["id"], "PENDING")}:
                        cancelled.append(candidate["id"])
                if cancelled:
                    command(["scancel", *cancelled])
            except RuntimeError as error:
                print("Optional cancellation of a finished stage's spare passes:", error, flush=True)
                return
        write_json(stamp, {"run_id": state["run_id"], "cancelled": cancelled})

def execute_stage(run, state, row, lane):
    env = os.environ.copy()
    env.update({"CLEAN_RUN_DIR": str(run), "CLEAN_RUN_ID": state["run_id"], "CLEAN_SOURCE_ROOT": str(ROOT),
                "JULIA_BIN": state["julia"], "JULIA_DEPOT_PATH": state["depot"], "CLEAN_PYTHON": sys.executable,
                "PYTHONDONTWRITEBYTECODE": "1",
                "JULIA_NUM_THREADS": str(CPUS), "OPENBLAS_NUM_THREADS": "1"})
    args = ["bash", str(ROOT / STAGE_PROGRAM), row["stage"], str(run), str(lane), str(row["lanes"]), str(SOFT)]
    with subprocess.Popen(args, cwd=run, env=env, start_new_session=True) as child:
        try:
            return child.wait(timeout=SOFT + 240)
        except subprocess.TimeoutExpired:
            os.killpg(child.pid, signal.SIGTERM)
            try:
                child.wait(timeout=45)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid, signal.SIGKILL)
                child.wait()
            print("Checkpoint window elapsed; retry will resume durable outputs.", flush=True)
            return 124

def worker(run, key):
    state = stored(run)
    require(state["committed"], "Run was never committed")
    bykey = {r["key"]: r for r in state["jobs"]}
    row = bykey[key]
    lane = int(os.environ.get("SLURM_ARRAY_TASK_ID", "1"))
    require(os.environ.get("SLURM_ARRAY_JOB_ID", os.environ.get("SLURM_JOB_ID")) == row["id"]
            and 1 <= lane <= row["lanes"], "Unexpected Slurm allocation ID or array lane")
    require(int(os.environ.get("SLURM_CPUS_PER_TASK", "0")) == CPUS, "Incorrect CPU allocation")
    stage = row["stage"]
    m = metadata(run)
    event = m / "attempts" / f"{key}.{lane}.json"
    receipt = {"run_id": state["run_id"], "manifest": state["manifest"], "stage": stage, "lane": lane,
               "job": row["id"], "pass": row["pass"], "started": now()}
    with lock(m / "locks" / f"{stage}.{lane}"):
        try:
            receipt["outcome"] = "running"
            write_json(event, receipt)
            prior = [p for p in bykey.values() if p["key"] in row["parents"] and p["stage"] != stage]
            if not prior:
                # Every retry has the same scientific prerequisites even when
                # its immediate Slurm dependency is an earlier pass of its own stage.
                first = next(x for x in bykey.values() if x["stage"] == stage)
                prior = [bykey[k] for k in first["parents"]]
            for previous_stage, lanes in {(p["stage"], p["lanes"]) for p in prior}:
                missing = [i for i in range(1, lanes + 1) if not complete(run, state, previous_stage, i)]
                require(not missing, f"Prerequisite {previous_stage} has unfinished lanes {missing}")
            already = complete(run, state, stage, lane)
            if already and stage == "final":
                require(final_receipt(run, state) == read(marker(run, stage, lane))["receipt_sha256"],
                        "Final verification receipt changed")
            if already:
                receipt["outcome"] = "already_complete"
                code = 0
            else:
                manifest(ROOT, state["manifest"], full=(stage == "final"))
                result = execute_stage(run, state, row, lane)
                receipt["science_exit"] = result
                if result == 0:
                    if stage == "final":
                        receipt["receipt_sha256"] = final_receipt(run, state)
                    receipt["outcome"] = "complete"
                    write_json(marker(run, stage, lane), receipt)
                    code = 0
                elif result in (75, 124):
                    receipt["outcome"] = "partial" if result == 75 else "timeout"
                    code = result
                else:
                    receipt["outcome"] = "science_error"
                    code = 1
            if code == 0:
                prune(run, state, row)
        except Exception as error:
            receipt.update(outcome="error", error=str(error))
            code = 1
        receipt["ended"] = now()
        write_json(event, receipt)
        print(json.dumps(receipt, sort_keys=True), flush=True)
        return code

def status(run):
    state = stored(run)
    rows = state["jobs"]
    print(f"Run {run}\nIdentity {state['run_id']}\nChain released: {state['released']}")
    for stage, lanes, _ in STAGES:
        done = sum(complete(run, state, stage, i) for i in range(1, lanes + 1))
        print(f"{stage}: {done}/{lanes} lanes complete")
    events = [read(p) for p in (metadata(run) / "attempts").glob("*.json")]
    output = run / "results"
    counts = (("Eigenfunction mode files", output / "eigenfunctions/modes", "*.jld2", 207),
              ("Residual column files", output / "residual/columns", "column_*.jld2", 804),
              ("Box column files", output / "coercivity/proposition/box_columns", "boxes_*.jld2", 2482),
              ("Response column files", output / "coercivity/table2_columns", "*/column_*.jld2", 67682))
    for label, folder, pattern, expected in counts:
        print(f"{label}: {sum(1 for _ in folder.glob(pattern))}/{expected}")
    response = output / "coercivity/table2_columns"
    print("Partial response checkpoints:", sum(1 for _ in response.glob("*/partial_*.jld2")))
    matrix = output / "coercivity/proposition/clean_matrix_blocks"
    print("Matrix strip files:", sum(1 for _ in matrix.glob("*/strip_*.jld2")))
    print("File counts measure saved work. Final assembly validates contents and coverage.")
    print("Recorded attempts:", dict(sorted(Counter(e["outcome"] for e in events).items())))
    for e in sorted(events, key=lambda x: x.get("started", "")):
        if e["outcome"] in ("error", "science_error", "timeout"):
            recovered = complete(run, state, e["stage"], e["lane"])
            print(f"  {e['job']} lane {e['lane']}, pass {e['pass']}, {e['outcome']}; lane recovered: {recovered}")
    ids = ",".join(r["id"] for r in rows if "id" in r)
    if ids:
        print("Live queue (array, task, state, CPUs, elapsed, priority, node or reason):")
        sq = subprocess.run(["squeue", "-h", "--me", "-o", "%F|%K|%T|%C|%M|%Q|%R"], text=True, capture_output=True)
        idset = set(ids.split(","))
        print("\n".join(l for l in sq.stdout.splitlines() if l.partition("|")[0] in idset) or sq.stderr.strip() or "(none)")
        print("Slurm accounting (array tasks, including hard Slurm timeouts and unused cancelled spares):")
        sa = subprocess.run(["sacct", "-n", "-X", "--array", "-j", ids, "-P", "-o", "JobID,State,ExitCode,Elapsed"], text=True, capture_output=True)
        counts = Counter(l.split("|")[1].strip() for l in sa.stdout.splitlines() if "|" in l)
        print(dict(sorted(counts.items())) if counts else sa.stderr.strip() or "(none reported)")
        for line in sa.stdout.splitlines():
            fields = line.split("|")
            if len(fields) < 4:
                continue
            job, outcome, exit_code, elapsed = fields[:4]
            if outcome.strip().partition(" ")[0] in {"FAILED", "TIMEOUT", "NODE_FAIL", "OUT_OF_MEMORY", "PREEMPTED", "BOOT_FAIL"}:
                # A worker which saved a prefix deliberately returns 75.
                if exit_code == "75:0":
                    continue
                print(f"  {job}: {outcome}, exit {exit_code}, elapsed {elapsed}")
    if complete(run, state, "final", 1):
        sha = final_receipt(run, state)
        require(sha == read(marker(run, "final", 1))["receipt_sha256"], "Final receipt digest changed")
        print(f"Verified final receipt {run / 'verification.toml'}\nSHA-256 {sha}")
    else:
        print("Final certification has no verified receipt yet.")

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    sub.add_parser("plan", help="offline scheduling plan and reserved capacity")
    for action in ("submit", "status", "worker"):
        p = sub.add_parser(action)
        p.add_argument("run", type=Path, help="absolute private run directory on scratch")
        if action == "worker":
            p.add_argument("key")
    args = parser.parse_args()
    if args.action == "plan":
        for stage, lanes, passes in STAGES:
            print(f"{stage}: {passes} passes of {lanes} x {CPUS}-CPU, two-hour allocations")
        print("Up to 12 simultaneous workers for this run. Productive checkpoint budget 110 minutes per worker/pass.")
        return 0
    run = args.run
    require(run.is_absolute() and not run.is_symlink(), "Use an absolute run directory; its final component must not be a symlink")
    run = run.resolve()
    if args.action == "submit":
        try:
            submit(run)
        except Exception:
            if (metadata(run) / "run.json").exists():
                print(f"Submission/release state is recorded in {metadata(run) / 'run.json'}. "
                      "After resolving the reported issue, repeat submit with this same run directory.",
                      file=sys.stderr)
            raise
    elif args.action == "status":
        status(run)
    else:
        return worker(run, args.key)
    return 0

if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, RuntimeError, ValueError, KeyError) as error:
        sys.exit(str(error))
