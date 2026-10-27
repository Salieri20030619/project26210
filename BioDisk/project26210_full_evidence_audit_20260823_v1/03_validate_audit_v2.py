#!/usr/bin/env python3
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CHECKS = []


def add(cid, ok, detail):
    CHECKS.append({"check_id": cid, "status": "PASS" if ok else "FAIL", "detail": str(detail)})


def rows(path):
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()


required = [
    "PROJECT26210_FULL_EVIDENCE_AUDIT_20260823_v1.md",
    "CLIENT_WETLAB_EVIDENCE_REGISTER_v1.tsv",
    "MECHANISM_CHAIN_EVIDENCE_REGISTER_v1.tsv",
    "AUDIT_SCAN_COUNTS_v1.md",
    "inventory_local_v1/SCAN_METADATA_v1.json",
    "inventory_remote_v2/SCAN_METADATA_v1.json",
]
for rel in required:
    p = ROOT / rel
    add("exists_" + rel.replace("/", "_"), p.is_file() and p.stat().st_size > 0, p)

wet = rows(ROOT / "CLIENT_WETLAB_EVIDENCE_REGISTER_v1.tsv")
chain = rows(ROOT / "MECHANISM_CHAIN_EVIDENCE_REGISTER_v1.tsv")
add("wet_rows", len(wet) == 7, len(wet))
add("wet_unique", len({r["evidence_id"] for r in wet}) == len(wet), len(wet))
add("chain_rows", len(chain) == 13, len(chain))
add("chain_unique", len({r["link_id"] for r in chain}) == len(chain), len(chain))
add("client_44x_statement_only", any(r["evidence_id"] == "CW02" and r["current_status"] == "STATEMENT_ONLY" for r in wet), "CW02")
add("full_chain_not_established", any(r["link_id"] == "L13" and r["project_decision"] == "FULL_CAUSAL_CHAIN_NOT_ESTABLISHED" for r in chain), "L13")

report = (ROOT / required[0]).read_text(encoding="utf-8-sig")
for cid, token in {
    "report_tymp_mta": "TYMP\u21925\u2032-MTA",
    "report_mta_spp1": "5\u2032-MTA\u2192SPP1",
    "report_phase70": "Phase70",
    "report_44x": "44 \u500d MTA",
    "report_boundaries": "\u5f53\u524d\u4e0d\u80fd\u5199",
}.items():
    add(cid, token in report, token)
add("report_public_sources", report.count("https://") >= 5, report.count("https://"))

for scope in ("inventory_local_v1", "inventory_remote_v2"):
    meta = json.loads((ROOT / scope / "SCAN_METADATA_v1.json").read_text(encoding="utf-8-sig"))
    add(scope + "_nonempty", meta["counts"]["files"] > 0 and meta["counts"]["errors"] == 0, meta["counts"])
    manifest = rows(ROOT / scope / "OUTPUT_SHA256_v1.tsv")
    ok = all((ROOT / scope / r["file"]).is_file() and sha(ROOT / scope / r["file"]) == r["sha256"] for r in manifest)
    add(scope + "_hashes", ok, len(manifest))

out = ROOT / "validation_v2"
out.mkdir(exist_ok=True)
with (out / "AUDIT_VALIDATION_CHECKS_v2.tsv").open("w", encoding="utf-8", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["check_id", "status", "detail"], delimiter="\t")
    w.writeheader()
    w.writerows(CHECKS)
summary = {"checks": len(CHECKS), "passed": sum(x["status"] == "PASS" for x in CHECKS), "failed": sum(x["status"] == "FAIL" for x in CHECKS)}
(out / "AUDIT_VALIDATION_SUMMARY_v2.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
if summary["failed"] == 0:
    (out / "AUDIT_VALIDATION_PASS_v2").write_text("PASS\n", encoding="ascii")
print(json.dumps(summary))
raise SystemExit(0 if summary["failed"] == 0 else 1)
