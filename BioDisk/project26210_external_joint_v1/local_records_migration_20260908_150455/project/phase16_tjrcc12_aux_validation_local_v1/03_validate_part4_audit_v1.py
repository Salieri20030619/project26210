#!/usr/bin/env python3
"""Validate the Part 4 reanalysis audit and mechanism registry."""

from __future__ import annotations

import csv
import hashlib
import os
from datetime import datetime
from pathlib import Path


HERE = Path(__file__).resolve().parent
AUDIT = HERE / "PART4_REANALYSIS_AUDIT_v1.tsv"
REGISTRY = HERE / "PART4_TYMP_MAINLINE_REANALYSIS_REGISTRY_v1.md"
OUTPUT = HERE / "PART4_REANALYSIS_AUDIT_VALIDATION_v1.tsv"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    with AUDIT.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fields = reader.fieldnames
    expected_fields = ["analysis_family", "old_dependency", "new_reference_impact", "action", "priority", "rationale"]
    allowed_actions = {
        "RETAIN_WITH_HASH_RECHECK", "RETAIN_AND_REMAP", "REPLACED_BY_NEW_MAPPING",
        "MUST_REDO", "MUST_REDO_MAINLINE", "MUST_REDO_FROM_REBUILD",
        "RETAIN_DOMAINS_REDO_ANNOTATION", "RETAIN_WITH_REINTERPRETATION",
        "REANALYZE_AS_SECONDARY", "REVIEW_THEN_SELECTIVE_RERUN", "RETAIN",
        "RETAIN_AS_AUXILIARY", "MUST_REBUILD_LAST",
    }
    registry = REGISTRY.read_text(encoding="utf-8")
    checks = [
        ("audit_schema", fields == expected_fields, f"fields={fields}"),
        ("audit_nonempty_unique", len(rows) == 20 and len({row['analysis_family'] for row in rows}) == len(rows), f"rows={len(rows)}"),
        ("audit_required_values", all(all(row[field].strip() for field in expected_fields) for row in rows), "no empty cells"),
        ("audit_action_enum", all(row["action"] in allowed_actions for row in rows), f"actions={sorted({row['action'] for row in rows})}"),
        ("audit_priority_enum", all(row["priority"] in {"P0", "P1", "P2", "P3"} for row in rows), f"priorities={sorted({row['priority'] for row in rows})}"),
        ("registry_customer_chain", all(term in registry for term in ["TYMP-high", "5'-MTA", "monocyte", "SPP1", "vascular", "immune exclusion", "therapy-resistance"]), "all main-line concepts present"),
        ("registry_f1_not_tymp", "F1-high/low is not a surrogate definition of TYMP-high/low" in registry, "identity boundary explicit"),
        ("registry_claim_boundaries", all(term in registry for term in ["transcriptomic context", "not lineage conversion proof", "Causal claims require"]), "measurement/lineage/causal boundaries explicit"),
        ("registry_execution_order", all(f"{index}." in registry for index in range(1, 7)), "six-step execution order present"),
    ]
    overall = all(passed for _, passed, _ in checks)
    checks.append(("OVERALL", overall, datetime.now().astimezone().isoformat(timespec="seconds")))
    part = Path(str(OUTPUT) + ".part")
    with part.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["check", "result", "detail"])
        for name, passed, detail in checks:
            writer.writerow([name, "PASS" if passed else "FAIL", detail])
        writer.writerow(["audit_sha256", digest(AUDIT), str(AUDIT)])
        writer.writerow(["registry_sha256", digest(REGISTRY), str(REGISTRY)])
    os.replace(part, OUTPUT)
    if not overall:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
