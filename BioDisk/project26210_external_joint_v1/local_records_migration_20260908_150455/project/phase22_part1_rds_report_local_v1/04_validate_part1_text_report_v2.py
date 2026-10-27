from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
REPORT = ROOT / "PART1_EXTERNAL_SINGLE_CELL_REPORT_TEXT_DRAFT_v2.md"
TABLE_DIR = ROOT / "reporting_tables_v2"
OUT_DIR = ROOT / "validation_v2"


def read_tsv(name: str) -> list[dict[str, str]]:
    with (TABLE_DIR / name).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    OUT_DIR.mkdir(exist_ok=True)
    text = REPORT.read_text(encoding="utf-8")
    status = json.loads((TABLE_DIR / "reporting_tables_status_v2.json").read_text(encoding="utf-8"))
    assoc = read_tsv("association_reporting_v2.tsv")

    full = [r for r in assoc if r["reporting_tier"] == "FULL_STATISTICS_GE10"]
    directional = [r for r in assoc if r["reporting_tier"] == "DIRECTION_ONLY_LT10"]
    not_tested = [r for r in assoc if r["reporting_tier"] == "NOT_TESTED"]

    checks = {
        "report_exists_and_nonempty": REPORT.exists() and len(text) > 3000,
        "source_tables_pass": status["status"] == "PASS_REPORTING_TABLES_COMPLETE",
        "patient_threshold_is_10": status["independent_unit"] == "patient"
        and status["minimum_units_for_inferential_reporting"] == 10,
        "association_tiers_match": (len(full), len(directional), len(not_tested)) == (42, 24, 54),
        "full_rows_all_ge10": all(int(r["n_independent_units"]) >= 10 for r in full),
        "direction_rows_all_lt10": all(int(r["n_independent_units"]) < 10 for r in directional),
        "no_fdr_supported_rows_ge10": status["fdr_supported_rows_ge10"] == 0,
        "new_reporting_rule_stated": "有效患者数 **≥10**" in text and "有效患者数 **<10**" in text,
        "no_small_n_inferential_values": all(token not in text for token in ["0.0188", "0.0416", "0.02797"]),
        "no_defensive_phrases": all(
            phrase not in text
            for phrase in ["不能证明", "不能替代", "不能表述", "尚不足以", "证据边界", "局限性"]
        ),
        "no_figures_embedded": "![" not in text and "<img" not in text.lower(),
        "required_conclusions_present": all(
            phrase in text
            for phrase in ["SPP1 髓系状态与血管重塑", "F1 与 C1QC", "TYMP—SPP1 关系不支持"]
        ),
    }

    overall = "PASS" if all(checks.values()) else "FAIL"
    sha256 = hashlib.sha256(REPORT.read_bytes()).hexdigest()
    rows = [{"check": key, "result": "PASS" if value else "FAIL"} for key, value in checks.items()]
    with (OUT_DIR / "PART1_TEXT_REPORT_VALIDATION_v2.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["check", "result"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    result = {
        "schema_version": "part1_text_report_validation_v2",
        "status": overall,
        "report": str(REPORT),
        "report_sha256": sha256,
        "checks_passed": sum(checks.values()),
        "checks_total": len(checks),
        "figures_generated": False,
        "review_state": "AWAITING_USER_REVIEW",
    }
    (OUT_DIR / "PART1_TEXT_REPORT_VALIDATION_STATUS_v2.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(f"{overall}: {sum(checks.values())}/{len(checks)} checks; sha256={sha256}")
    if overall != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
