#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import pandas as pd


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--root", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args()
    if args.output.exists(): raise FileExistsError(args.output)
    part = args.output.with_name(args.output.name + ".part")
    if part.exists(): raise FileExistsError(part)
    part.mkdir(parents=True)
    checks = []
    def add(name, ok, detail): checks.append({"check":name,"status":"PASS" if bool(ok) else "FAIL","detail":str(detail)})
    components = [
        ("M11", "m11_validation_v2", 11),
        ("M12", "m12_validation_v1", 10),
        ("M13", "m13_validation_v1", 14),
    ]
    passed_sum = 0; total_sum = 0
    for module, directory, expected in components:
        status = json.loads((args.root / directory / "status.json").read_text(encoding="utf-8"))
        add(f"{module}_validation", status.get("status") == "PASS_INDEPENDENT_VALIDATION" and status.get("checks_passed") == expected and status.get("checks_total") == expected, status)
        passed_sum += int(status.get("checks_passed", 0)); total_sum += int(status.get("checks_total", 0))
    add("component_check_total", passed_sum == 35 and total_sum == 35, f"{passed_sum}/{total_sum}")
    m11 = pd.read_csv(args.root / "m11_analysis_v1" / "m11_final_mr_target_inference_status_v1.tsv", sep="\t")
    m12 = pd.read_csv(args.root / "m12_analysis_v1" / "m12_feature_registry_v1.tsv", sep="\t")
    m13 = pd.read_csv(args.root / "m13_analysis_v1" / "m13_feature_registry_v1.tsv", sep="\t")
    infer = set(m11.loc[m11.activity_row_present.astype(bool), "protein"])
    add("M11_M12_feature_continuity", infer == set(m12.loc[m12.evidence_family == "final_MR_inferable", "protein"]), infer)
    add("M12_M13_feature_continuity", infer == set(m13.loc[m13.evidence_family == "final_MR_inferable", "feature"]), infer)
    m11_status = json.loads((args.root / "m11_analysis_v1" / "status.json").read_text(encoding="utf-8"))
    m12_status = json.loads((args.root / "m12_analysis_v1" / "status.json").read_text(encoding="utf-8"))
    m13_status = json.loads((args.root / "m13_analysis_v1" / "status.json").read_text(encoding="utf-8"))
    add("common_scope", all(x.get("n_cells") == 11921 and x.get("n_patients") == 20 for x in (m11_status,m12_status,m13_status)), "11921 cells;20 patients")
    add("M13_scope_limitation", m13_status.get("TYMP_RRM2B_activity_available") is False and m13_status.get("TYMP_RRM2B_RNA_sensitivity") is True, "activity unavailable; RNA separate")
    analysis_dirs = [args.root / x for x in ("m11_analysis_v1","m12_analysis_v1","m13_expression_v1","m13_analysis_v1","m11_validation_v2","m12_validation_v1","m13_validation_v1")]
    add("no_figures", not any(any(d.rglob(ext)) for d in analysis_dirs for ext in ("*.png","*.pdf")), "tables only")
    frame = pd.DataFrame(checks)
    frame.to_csv(part / "terminal_checks_v1.tsv", sep="\t", index=False)
    failed = frame.loc[frame.status != "PASS", "check"].tolist()
    manifest=[]
    for directory in analysis_dirs:
        for path in sorted(directory.rglob("*")):
            if path.is_file(): manifest.append({"path":str(path),"bytes":path.stat().st_size,"sha256":sha256(path)})
    pd.DataFrame(manifest).to_csv(part / "VALIDATION_MANIFEST_v1.tsv", sep="\t", index=False)
    status={"status":"PASS_PHASE29_TERMINAL_VALIDATION" if not failed else "FAIL_PHASE29_TERMINAL_VALIDATION","terminal_checks_passed":int((frame.status=="PASS").sum()),"terminal_checks_total":len(frame),"component_checks_passed":passed_sum,"component_checks_total":total_sum,"failed_checks":failed}
    (part/"status.json").write_text(json.dumps(status,indent=2)+"\n",encoding="utf-8")
    os.replace(part,args.output)
    if failed: raise SystemExit("terminal failures: "+",".join(failed))


if __name__ == "__main__": main()

