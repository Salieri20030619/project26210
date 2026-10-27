#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path


ROOT = Path("/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/05_cnmf")
OUT = ROOT / "cnv_clone_program_function_v1"
VERIFY = ROOT / "cnv_clone_program_function_verification_v1"


def rows(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manual-pass", action="store_true")
    args = parser.parse_args()
    VERIFY.mkdir(parents=True, exist_ok=True)
    checks: list[tuple[str, bool, str]] = []

    def check(name: str, passed: bool, detail: str) -> None:
        checks.append((name, bool(passed), detail))

    required = [
        "RUN_COMPLETE", "run_metadata_v1.json", "file_manifest_v1.tsv",
        "pseudobulk_sample_qc_v1.tsv",
        "F1_high_clone_vs_other_paired_pseudobulk_DE_v1.tsv.gz",
        "rank7_program_pathway_ORA_v1.tsv", "F1_high_clone_DE_pathway_ORA_v1.tsv",
        "rank7_program_top100_genes_v1.tsv",
    ]
    for name in required:
        check(f"exists::{name}", (OUT / name).is_file(), str(OUT / name))
    if not all(p for _, p, _ in checks):
        raise SystemExit("Required output missing")

    metadata = json.loads((OUT / "run_metadata_v1.json").read_text())
    check("n_input_cells", metadata.get("n_input_cells") == 11950, str(metadata.get("n_input_cells")))
    check("n_evaluable_patients", metadata.get("n_evaluable_patients") == 14, str(metadata.get("n_evaluable_patients")))
    check("n_pseudobulk_samples", metadata.get("n_pseudobulk_samples") == 28, str(metadata.get("n_pseudobulk_samples")))
    check("n_rank7_programs", metadata.get("n_rank7_programs") == 7, str(metadata.get("n_rank7_programs")))

    samples = list(rows(OUT / "pseudobulk_sample_qc_v1.tsv"))
    patients = {r["curated_patient_id"] for r in samples}
    check("pseudobulk_28", len(samples) == 28, str(len(samples)))
    check("patients_14", len(patients) == 14, str(len(patients)))
    check("two_groups_per_patient", all(sum(r["curated_patient_id"] == p for r in samples) == 2 for p in patients), "paired")
    check("positive_cells_libraries", all(int(r["n_cells"]) > 0 and float(r["library_size"]) > 0 for r in samples), "all positive")

    de = list(rows(OUT / "F1_high_clone_vs_other_paired_pseudobulk_DE_v1.tsv.gz"))
    check("de_nonempty", len(de) > 1000, str(len(de)))
    check("de_unique_genes", len({r["gene"] for r in de}) == len(de), str(len(de)))
    check("de_p_ranges", all(0 <= float(r["edgeR_PValue"]) <= 1 and 0 <= float(r["edgeR_FDR"]) <= 1 for r in de), "[0,1]")

    ora = list(rows(OUT / "rank7_program_pathway_ORA_v1.tsv"))
    check("seven_programs_in_ora", {r["analysis"] for r in ora} == {f"F{i}" for i in range(1, 8)}, str(sorted({r["analysis"] for r in ora})))
    check("ora_p_ranges", all(0 <= float(r["PValue"]) <= 1 and 0 <= float(r["FDR"]) <= 1 for r in ora), "[0,1]")

    manifest = list(rows(OUT / "file_manifest_v1.tsv"))
    bad_manifest = []
    for record in manifest:
        path = OUT / record["file"]
        if not path.is_file() or path.stat().st_size != int(record["size_bytes"]) or sha256(path) != record["sha256"]:
            bad_manifest.append(record["file"])
    check("manifest_integrity", not bad_manifest, ",".join(bad_manifest) or "all matched")

    figures = [
        "rank7_program_pathway_ORA_v1", "F1_high_clone_pseudobulk_DE_volcano_v1",
        "F1_high_clone_top_DE_patient_heatmap_v1",
    ]
    for stem in figures:
        png, pdf = OUT / "figures" / f"{stem}.png", OUT / "figures" / f"{stem}.pdf"
        check(f"figure::{stem}", png.stat().st_size > 10000 and pdf.stat().st_size > 1000, f"png={png.stat().st_size};pdf={pdf.stat().st_size}")

    with (VERIFY / "verification_checks_v1.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["check", "status", "detail"])
        writer.writerows((n, "PASS" if p else "FAIL", d) for n, p, d in checks)
    if args.manual_pass:
        with (VERIFY / "manual_figure_review_v1.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["figure", "status", "review_scope"])
            writer.writerows((stem, "PASS", "nonblank, labels legible, geometry coherent") for stem in figures)
    summary = {"schema_version": "part1_cnv_cnmf_function_verification_v1", "n_checks": len(checks), "n_failed": sum(not p for _, p, _ in checks), "manual_review": bool(args.manual_pass)}
    (VERIFY / "verification_metadata_v1.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    if summary["n_failed"]:
        raise SystemExit(f"Verification failed: {summary['n_failed']}")
    (VERIFY / "RUN_COMPLETE").write_text("PASS\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
