#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results_local_v1"
VALIDATION = ROOT / "validation_local_v1"
OUT = ROOT / "terminal_validation_local_v1"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    if OUT.exists():
        raise RuntimeError(f"Refusing to overwrite existing terminal output: {OUT}")
    OUT.mkdir(parents=False)

    checks: list[tuple[str, bool, str]] = []
    v = json.loads((VALIDATION / "status.json").read_text(encoding="utf-8"))
    checks.append(("independent_validation_pass", v["status"] == "PASS_INDEPENDENT_VALIDATION", v["status"]))
    checks.append(("independent_validation_16_of_16", v["checks_passed"] == 16 and v["checks_total"] == 16, f"{v['checks_passed']}/{v['checks_total']}"))

    s = json.loads((RESULTS / "status.json").read_text(encoding="utf-8"))
    checks.append(("analysis_cell_count", s["n_cells"] == 11950, str(s["n_cells"])))
    checks.append(("analysis_gene_count", s["n_genes"] == 12, str(s["n_genes"])))
    checks.append(("analysis_pair_count", s["n_gene_pairs"] == 66, str(s["n_gene_pairs"])))
    checks.append(("two_sided_no_prespecified_direction", s["two_sided_testing"] is True and s["direction_prespecified"] is False, "two-sided; no direction prespecified"))
    checks.append(("no_figures", s["figures_created"] is False and not any(ROOT.rglob("*.png")) and not any(ROOT.rglob("*.pdf")), "no PNG/PDF"))

    integrated = read_tsv(RESULTS / "m17_gene_integrated_evidence_v1.tsv")
    pairs = read_tsv(RESULTS / "m17_gene_pair_correlation_summary_v1.tsv")
    checks.append(("integrated_gene_rows", len(integrated) == 12, str(len(integrated))))
    checks.append(("pair_rows", len(pairs) == 66, str(len(pairs))))
    secondary = [r["gene"] for r in integrated if r["evidence_tier"] == "SECONDARY_ASSOCIATION"]
    strong = [r["gene"] for r in integrated if r["evidence_tier"] == "STRONG_ASSOCIATION"]
    supported_pairs = [f"{r['gene_a']}--{r['gene_b']}" for r in pairs if r["statistical_support"] == "True"]
    checks.append(("secondary_gene_sentinel", secondary == ["ODC1"], ",".join(secondary)))
    checks.append(("no_strong_gene", not strong, ",".join(strong) if strong else "none"))
    checks.append(("supported_pair_sentinels", set(supported_pairs) == {"TYMP--SLC29A1", "TYMP--SRM", "TK1--SMS", "SRM--PCYT2"}, ",".join(supported_pairs)))

    failed = [name for name, ok, _ in checks if not ok]
    check_path = OUT / "terminal_checks_v1.tsv.part"
    with check_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["check", "pass", "detail"])
        writer.writerows((name, str(ok), detail) for name, ok, detail in checks)
    check_path.replace(OUT / "terminal_checks_v1.tsv")

    manifest_rows = []
    for path in sorted(list(RESULTS.glob("*")) + list(VALIDATION.glob("*"))):
        if path.is_file():
            manifest_rows.append((str(path.relative_to(ROOT)), path.stat().st_size, sha256(path)))
    manifest_path = OUT / "validated_file_manifest_v1.tsv.part"
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["relative_path", "bytes", "sha256"])
        writer.writerows(manifest_rows)
    manifest_path.replace(OUT / "validated_file_manifest_v1.tsv")

    status = {
        "status": "PASS_PHASE31_M17_TERMINAL_VALIDATION" if not failed else "FAIL_PHASE31_M17_TERMINAL_VALIDATION",
        "timestamp_local": datetime.now().astimezone().isoformat(timespec="seconds"),
        "checks_passed": len(checks) - len(failed),
        "checks_total": len(checks),
        "independent_component_checks": "16/16",
        "failed_checks": failed,
        "n_manifest_files": len(manifest_rows),
        "figures_created": False,
    }
    status_path = OUT / "status.json.part"
    status_path.write_text(json.dumps(status, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    status_path.replace(OUT / "status.json")
    if failed:
        raise RuntimeError(f"Terminal validation failed: {failed}")


if __name__ == "__main__":
    main()
