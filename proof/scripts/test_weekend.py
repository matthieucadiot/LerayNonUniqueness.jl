"""Exercise the submission chain, worker retries and acceptance using a fake Slurm."""
import contextlib
import importlib.util
import io
import os
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("weekend", HERE / "weekend.py")
weekend = importlib.util.module_from_spec(spec)
spec.loader.exec_module(weekend)

# Different stages/lanes write into the owned run only. Changing answer files
# simulates soft continuation, failed computation and an invalid final receipt.
FAKE_SCIENCE = """#!/bin/bash
set -eu
test "$2" = "$CLEAN_RUN_DIR"
test "$5" = 6600
test "$JULIA_NUM_THREADS" = 32
test -x "$JULIA_BIN"
test -d "$JULIA_DEPOT_PATH"
key="$1.$3"
bump="$2/visited.$key"
echo x >> "$bump"
if [ -f "$2/answer.$key" ]; then
    code=$(cat "$2/answer.$key")
    if [ "$code" != 0 ]; then exit "$code"; fi
fi
if [ "$1" = final ]; then
    if [ -e "$2/invalid" ]; then identity=previous; else identity="$CLEAN_RUN_ID"; fi
    printf 'all_checks_passed = true\\nrun_id = "%s"\\n' "$identity" > "$2/verification.toml"
fi
"""
SLURM = """#!/usr/bin/env python3
import json, os, sys
from pathlib import Path
db = Path(os.environ['FAKE_SLURM'])
state = json.loads(db.read_text())
args = sys.argv[1:]
prog = Path(sys.argv[0]).name
def flag(key):
    return next(i.partition('=')[2] for i in args if i.startswith(key+'='))
if prog == 'sbatch':
    job = str(10000 + len(state['jobs']))
    data = {'id':job, 'name':flag('--job-name'), 'args':args, 'state':'PENDING', 'held':'--hold' in args}
    state['jobs'][job] = data
    print(job + ';della')
    if state.pop('lose', 0) and len(state['jobs']) == 1:
        db.write_text(json.dumps(state))
        sys.exit(1)
elif prog == 'scontrol':
    data = state['jobs'][args[-1]]
    if 'release' in args:
        data['held'] = False
        state['releases'].append((data['id'],len(state['jobs'])))
    else:
        reason = 'JobHeldUser' if data['held'] else 'Priority'
        print('JobState=PENDING Reason='+reason)
elif prog == 'squeue':
    for id,data in state['jobs'].items():
        if data['state'] == 'PENDING':
            print('|'.join((id,data['name'],data['state'])))
elif prog == 'scancel':
    for id in args:
        state['jobs'][id]['state'] = 'CANCELLED'
elif prog == 'sacct':
    pass
else:
    sys.exit('unsupported '+prog)
db.write_text(json.dumps(state))
"""

class DellaTests(unittest.TestCase):
    def setUp(self):
        self.sandbox = tempfile.TemporaryDirectory()
        self.addCleanup(self.sandbox.cleanup)
        base = Path(self.sandbox.name).resolve()
        source = base / "source"
        self.source = source
        self.run = base / "run"
        (source / "scripts").mkdir(parents=True)
        shutil.copyfile(HERE / "weekend.py", source / "scripts/weekend.py")
        shutil.copyfile(HERE / "weekend_job.sh", source / "scripts/weekend_job.sh")
        (source / weekend.STAGE_PROGRAM).write_text(FAKE_SCIENCE)
        (source / "fixed.jld2").write_text("reviewed input")
        files = sorted(f for f in source.rglob("*") if f.is_file())
        (source / weekend.MANIFEST).write_text("".join(f"{weekend.digest(f)}  {f.relative_to(source)}\n" for f in files))
        bindir = base / "bin"
        bindir.mkdir()
        for cmd in ("scontrol", "sbatch", "squeue", "sacct", "scancel", "julia"):
            f = bindir / cmd
            f.write_text("#!/bin/bash\necho 'julia version 1.12.3'\n" if cmd == "julia" else SLURM)
            f.chmod(0o700)
        self.db = base / "slurm.json"
        self.db.write_text('{"jobs": {}, "releases": []}')
        self.patches = [patch.object(weekend, "ROOT", source), patch.dict(os.environ, {
            "FAKE_SLURM": str(self.db), "JULIA_BIN": str(bindir / "julia"), "JULIA_DEPOT_PATH": str(base),
            "PATH": str(bindir) + os.pathsep + os.environ["PATH"]}),
            patch.object(weekend, "STAGES", (("base", 2, 3), ("assemble_base", 1, 1), ("main", 2, 2), ("final", 1, 2)))]
        for p in self.patches:
            p.start()
            self.addCleanup(p.stop)

    def submit(self):
        with contextlib.redirect_stdout(io.StringIO()):
            weekend.submit(self.run)
        return weekend.stored(self.run)

    def worker(self, state, key, lane):
        row = next(r for r in state["jobs"] if r["key"] == key)
        with patch.dict(os.environ, {"SLURM_ARRAY_JOB_ID": row["id"], "SLURM_ARRAY_TASK_ID": str(lane),
                                     "SLURM_CPUS_PER_TASK": "32"}), contextlib.redirect_stdout(io.StringIO()):
            return weekend.worker(self.run, key)

    def event(self, key, lane):
        return weekend.read(weekend.metadata(self.run) / "attempts" / f"{key}.{lane}.json")

    def test_all_ids_committed_and_retrying_submit_is_idempotent(self):
        self.db.write_text('{"jobs": {}, "releases": [], "lose": 1}')
        with self.assertRaisesRegex(RuntimeError, "sbatch"):
            self.submit()
        self.assertEqual(len(weekend.read(self.db)["jobs"]),1)
        state = self.submit()
        self.submit()
        sim = weekend.read(self.db)
        jobs = state["jobs"]
        self.assertEqual(len(jobs),8)
        self.assertEqual(sim["releases"], [[jobs[0]["id"],len(jobs)]])
        for r in jobs:
            args = sim["jobs"][r["id"]]["args"]
            self.assertIn("--cpus-per-task=32", args)
            self.assertIn("--time=02:00:00", args)
            self.assertIn("--qos=short", args)
            if r["lanes"] > 1:
                self.assertIn("--array=1-2%2", args)
            parent_ids = [x["id"] for x in jobs if x["key"] in r["parents"]]
            if parent_ids:
                self.assertIn("--dependency=afterany:" + ":".join(parent_ids), args)
            self.assertIn(r["id"], (self.run/"BOOST_REQUEST.txt").read_text())
        self.assertEqual(jobs[3]["parents"], [r["key"] for r in jobs[:3]])

    def test_failure_resume_prune_and_scientific_acceptance(self):
        state = self.submit()
        self.assertEqual(self.worker(state, "assemble_base.01", 1),1)
        self.assertEqual(self.worker(state, "base.01", 1),0)
        self.assertEqual(self.worker(state, "base.02", 1),0)
        self.assertEqual((self.run/"visited.base.1").read_text(),"x\n")
        answer = self.run/"answer.base.2"
        answer.write_text("75")
        self.worker(state,"base.01",2)
        self.assertEqual(self.event("base.01",2)["outcome"], "partial")
        self.assertFalse(weekend.complete(self.run, state, "base", 2))
        answer.unlink()
        self.assertEqual(self.worker(state,"base.02",2),0)
        self.assertEqual([r["id"] for r in weekend.read(self.db)["jobs"].values() if r["state"] == "CANCELLED"],
                         [state["jobs"][2]["id"]])
        self.assertEqual(self.worker(state,"assemble_base.01",1),0)
        (self.run/"answer.main.1").write_text("1")
        self.assertEqual(self.worker(state,"main.01",1),1)
        self.assertEqual(self.event("main.01",1)["outcome"], "science_error")
        (self.run/"answer.main.1").unlink()
        self.assertEqual(self.worker(state,"main.02",1),0)
        self.assertEqual(self.worker(state,"main.02",2),0)
        (self.run/"invalid").touch()
        self.assertEqual(self.worker(state,"final.01",1),1)
        self.assertFalse(weekend.complete(self.run,state,"final",1))
        (self.run/"invalid").unlink()
        self.assertEqual(self.worker(state,"final.02",1),0)
        self.assertTrue(weekend.complete(self.run,state,"final",1))
        (self.run/"verification.toml").write_text("all_checks_passed=false\n")
        with self.assertRaisesRegex(RuntimeError,"does not certify"), contextlib.redirect_stdout(io.StringIO()):
            weekend.status(self.run)

    def test_status_counts_native_results_and_allows_standard_julia_depot(self):
        with patch.dict(os.environ, {"JULIA_DEPOT_PATH": os.environ["JULIA_DEPOT_PATH"] + os.pathsep}):
            state = self.submit()
        self.assertTrue(state["depot"].endswith(os.pathsep))
        for name in ("coercivity/proposition/box_columns/boxes_1.jld2",
                     "coercivity/proposition/clean_matrix_blocks/odd_T0U0/strip_2.jld2",
                     "coercivity/table2_columns/odd_1/column_2.jld2",
                     "residual/columns/column_1.jld2"):
            f = self.run / "results" / name
            f.parent.mkdir(parents=True, exist_ok=True)
            f.touch()
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            weekend.status(self.run)
        for line in ("Box column files: 1/2482", "Matrix strip files: 1",
                     "Response column files: 1/67682", "Residual column files: 1/804"):
            self.assertIn(line, output.getvalue())

    def test_integrity_and_misrouted_work_fail_before_computation(self):
        state = self.submit()
        (self.source/weekend.STAGE_PROGRAM).write_text("exit 0")
        self.assertEqual(self.worker(state,"base.01",1),1)
        self.assertIn("Changed input", self.event("base.01",1)["error"])
        self.assertFalse((self.run/"visited.base.1").exists())
        with patch.dict(os.environ, {"SLURM_ARRAY_JOB_ID":"999","SLURM_ARRAY_TASK_ID":"1","SLURM_CPUS_PER_TASK":"32"}):
            with self.assertRaisesRegex(RuntimeError,"Unexpected Slurm allocation"):
                weekend.worker(self.run, "base.01")

if __name__ == "__main__":
    unittest.main()
