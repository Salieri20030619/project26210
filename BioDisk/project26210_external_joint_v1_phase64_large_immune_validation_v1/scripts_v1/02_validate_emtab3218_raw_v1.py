#!/usr/bin/env python3
import csv
import json
import os
import sys


def walk(obj):
    if isinstance(obj, dict):
        yield obj
        for value in obj.values():
            yield from walk(value)
    elif isinstance(obj, list):
        for value in obj:
            yield from walk(value)


def main(manifest_path, extraction_dir, output_dir):
    with open(manifest_path, encoding="utf-8") as handle:
        manifest = json.load(handle)
    expected = {
        node["path"]: int(node["size"])
        for node in walk(manifest)
        if node.get("type") == "file" and str(node.get("path", "")).lower().endswith(".cel")
    }
    observed = {}
    duplicate_names = []
    for root, _, files in os.walk(extraction_dir):
        for name in files:
            if not name.lower().endswith(".cel"):
                continue
            if name in observed:
                duplicate_names.append(name)
            observed[name] = os.path.join(root, name)

    rows = []
    for name in sorted(set(expected) | set(observed)):
        path = observed.get(name, "")
        observed_size = os.path.getsize(path) if path else None
        expected_size = expected.get(name)
        status = "PASS" if expected_size is not None and observed_size == expected_size else "FAIL"
        rows.append({
            "file": name,
            "expected_bytes": "" if expected_size is None else str(expected_size),
            "observed_bytes": "" if observed_size is None else str(observed_size),
            "status": status,
        })

    output_path = os.path.join(output_dir, "raw_cel_manifest_v1.tsv")
    with open(output_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    gates = {
        "expected_cel_count": (len(expected), 114),
        "observed_cel_count": (len(observed), 114),
        "duplicate_cel_basenames": (len(duplicate_names), 0),
        "size_mismatches": (sum(row["status"] != "PASS" for row in rows), 0),
    }
    gate_path = os.path.join(output_dir, "raw_cel_gate_v1.tsv")
    with open(gate_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["gate", "observed", "expected", "status"])
        for gate, (observed_value, expected_value) in gates.items():
            status = "PASS" if observed_value == expected_value else "FAIL"
            writer.writerow([gate, observed_value, expected_value, status])
            print(gate, observed_value, expected_value, status, sep="\t")
    if any(observed_value != expected_value for observed_value, expected_value in gates.values()):
        raise SystemExit(2)


if __name__ == "__main__":
    main(*sys.argv[1:4])
