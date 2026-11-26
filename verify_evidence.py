"""Check the published November 26 records locally using Python 3.11 or later."""

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
import sys
import tomllib


RECEIPT = "42fc3d336755a2c10b9d1338ba298e66bc86b7b1cbe8a3be32451b2703292104"
RUN = "2d9cec66690641fa88a3aac70215c5cb"
TABLE2 = {
    "odd": {"q": "19.2", "q2": "19.86", "qe": "0.0005", "Φq": "15.25"},
    "even": {"q": "34", "q2": "61", "qe": "0.001", "Φq": "146"},
}


def check(condition, message):
    if not condition:
        raise ValueError(message)


def check_hash(path, expected):
    with path.open("rb") as source:
        actual = hashlib.file_digest(source, "sha256").hexdigest()
    check(actual == expected, f"SHA-256 mismatch for {path}")


def read_toml(path):
    with path.open("rb") as source:
        return tomllib.load(source)


def verify(root, proof):
    path = root / "verification.toml"
    check_hash(path, RECEIPT)
    v = read_toml(path)
    check(v["run_id"] == RUN and v["all_checks_passed"] is True, "Run identity or verdict differs")
    for relative in [".weekend/done/final.1.json", ".weekend/attempts/final.01.1.json"]:
        marker = json.loads((root / relative).read_text())
        check(marker["run_id"] == RUN and marker["receipt_sha256"] == RECEIPT
              and marker["manifest"] == v["package_sha256"] and marker["stage"] == "final"
              and marker["job"] == "18207736" and marker["outcome"] == "complete"
              and marker["science_exit"] == 0, f"Unexpected final job record {relative}")

    results = root / "results"
    co = results / "coercivity"
    prop = co / "proposition"
    check_hash(results / "eigenfunctions/Eigenfunctions.jld2", v["upstream"]["eigenfunctions_sha256"])
    check_hash(results / "residual/Residual_Table1.jld2", v["upstream"]["residual_sha256"])
    check_hash(co / "Table2_bounds.jld2", v["table2_data_sha256"])
    check_hash(prop / "Boxes.jld2", v["proposition"]["boxes_sha256"])
    check_hash(prop / "Bounds.jld2", v["proposition"]["bounds_sha256"])

    for name, file in [("table2", co / "table2_checks.toml"), ("table4", prop / "table4_checks.toml")]:
        check_hash(file, v[name + "_sha256"])
        report = read_toml(file)
        check(report == v[name] and report["run_id"] == RUN and report["all_checks_passed"] is True,
              f"Component report differs: {name}")

    coverage = v["table2"]["coverage"]
    check((coverage["columns"], coverage["parent_cells"], coverage["selected_columns"])
          == (67682, 57725138, 4289), "Unexpected response coverage")
    for parity, thresholds in TABLE2.items():
        for field, target in thresholds.items():
            entry = v["table2"][parity][field]
            check(entry["passed"] is True and 0 <= Fraction(entry["upper"]) <= Fraction(target),
                  f"Table 2 bound failed for {parity} {field}")

    matrices = v["table4"]["rows"]
    check(len(matrices) == v["table4"]["matrix_count"] == 15, "Incomplete Table 4")
    for name, row in matrices.items():
        check_hash(prop / "matrices" / (name + ".jld2"), row["result_sha256"])
        check(row["verified"] is True and 0 <= row["norm_upper"] <= 30
              and 0 <= Fraction(row["loss_upper"]) <= Fraction("0.00000000009")
              and row["minimum_margin_lower"] >= 0, f"Table 4 bound failed for {name}")

    check(len(v["closure"]) == 7, "Missing closing inequalities")
    for name, entry in v["closure"].items():
        value = Fraction(entry["exact_value"].replace("//", "/"))
        target = Fraction(entry["exact_target"].replace("//", "/"))
        check(value < target, f"Closing inequality failed for {name}")

    if proof:
        manifest = proof / "CLEAN_SHA256SUMS"
        check_hash(manifest, v["package_sha256"])
        for line in manifest.read_text().splitlines():
            sha, relative = line.split("  ", 1)
            check_hash(proof / relative, sha)
        check_hash(proof / "data/Coercivity/Fixed_weighted_responses.jld2", v["fixed_weights_sha256"])
        print("Verified the original code and fixed input manifest.")
    print("Verified the final receipt, 15 matrix files, component reports and aggregate certificates.")
    print("Table 2 bounds, Table 4 recorded margins and seven exact closing comparisons passed.")
    print(f"Run {RUN}; {coverage['parent_cells']:,} parent cells in {coverage['columns']:,} columns.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("evidence", type=Path, nargs="?", default=Path(__file__).resolve().parent,
                        help="extracted certification-evidence directory")
    parser.add_argument("--proof", type=Path, help="also verify the extracted full proof package")
    args = parser.parse_args()
    try:
        verify(args.evidence, args.proof)
    except (OSError, KeyError, TypeError, ValueError) as error:
        sys.exit(f"Evidence check failed: {error}")
