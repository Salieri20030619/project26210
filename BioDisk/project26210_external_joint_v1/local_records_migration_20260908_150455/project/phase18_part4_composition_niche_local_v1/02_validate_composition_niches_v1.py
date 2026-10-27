#!/usr/bin/env python3
"""Independent structural/conservation validation for Phase 18."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd


EXPECTED_STATES = [
    "B cell", "C1QC-state Macrophage", "Dendritic cell", "Endothelial",
    "Epithelial-Normal", "Malignant", "Mast cell", "Monocyte", "NK cell",
    "Other Macrophage", "SMC/Pericyte", "SPP1-state Macrophage", "T cell",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
    part = Path(str(path) + ".part")
    frame.to_csv(part, sep="\t", index=False, lineterminator="\n")
    os.replace(part, path)


def atomic_json(path: Path, payload: dict) -> None:
    part = Path(str(path) + ".part")
    part.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    os.replace(part, path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-root", type=Path, required=True)
    parser.add_argument("--phase17-audit", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    root = args.analysis_root.resolve(strict=True)
    audit_path = args.phase17_audit.resolve(strict=True)
    output = args.output_root.resolve()
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite validation output: {output}")
    output.mkdir(parents=True)
    checks = []

    def check(name: str, passed: bool, detail: str) -> None:
        checks.append({"check": name, "status": "PASS" if passed else "FAIL", "detail": detail})

    status = json.loads((root / "PHASE18_COMPOSITION_NICHE_STATUS_v1.json").read_text(encoding="utf-8"))
    check("formal_status", status.get("status") == "PASS_ANALYSIS_COMPLETE", str(status.get("status")))
    parts = list(root.rglob("*.part"))
    check("no_part_files", len(parts) == 0, f"count={len(parts)}")

    manifest = pd.read_csv(root / "OUTPUT_MANIFEST.tsv", sep="\t")
    manifest_ok = True
    for row in manifest.itertuples(index=False):
        path = root / row.file
        manifest_ok &= path.is_file() and path.stat().st_size == int(row.size_bytes) and sha256(path) == row.sha256
    check("manifest_size_sha256", bool(manifest_ok), f"files={len(manifest)}")

    audit = pd.read_csv(audit_path, sep="\t")
    audit_ok = len(audit) == 14 and audit["sample"].nunique() == 14 and audit["status"].eq("completed").all() and audit["output_sha256_matches"].eq(True).all()
    check("phase17_14_of_14_gate", bool(audit_ok), f"rows={len(audit)}")

    tables = root / "tables"
    wide = pd.read_csv(tables / "spot_composition_wide.tsv.gz", sep="\t")
    provenance = pd.read_csv(tables / "input_provenance.tsv", sep="\t")
    state_summary = pd.read_csv(tables / "sample_state_summary.tsv", sep="\t")
    niche_summary = pd.read_csv(tables / "niche_summary.tsv", sep="\t")
    niche_comp = pd.read_csv(tables / "niche_state_composition.tsv", sep="\t")
    adjacency = pd.read_csv(tables / "niche_adjacency_descriptive.tsv", sep="\t")
    k_table = pd.read_csv(tables / "niche_k_selection.tsv", sep="\t")

    check("section_enumeration", wide["sample"].nunique() == 14 and provenance["sample"].nunique() == 14,
          f"wide={wide['sample'].nunique()}, provenance={provenance['sample'].nunique()}")
    keys_unique = not wide.duplicated(["sample", "spot_id"]).any()
    spot_counts = wide.groupby("sample").size().sort_index()
    provenance_counts = provenance.set_index("sample")["n_spots"].sort_index()
    check("spot_identity_conservation", keys_unique and spot_counts.equals(provenance_counts), f"spots={len(wide)}")

    fraction_cols = [f"fraction__{x}" for x in EXPECTED_STATES]
    exact_columns = all(x in wide.columns for x in fraction_cols) and len([x for x in wide.columns if x.startswith("fraction__")]) == 13
    fractions = wide[fraction_cols].to_numpy(float)
    fraction_ok = exact_columns and np.isfinite(fractions).all() and (fractions >= 0).all() and np.allclose(fractions.sum(axis=1), 1.0, atol=2e-7)
    check("exact_13_state_fractions", bool(fraction_ok), f"columns={len(fraction_cols)}")
    dominant = np.asarray(EXPECTED_STATES)[fractions.argmax(axis=1)]
    check("dominant_state_recomputed", np.array_equal(dominant, wide["dominant_state"].to_numpy()), "argmax of fractions")

    state_set_ok = len(state_summary) == 14 * 13 and set(state_summary["state"]) == set(EXPECTED_STATES)
    recomputed = wide.groupby("sample")[fraction_cols].mean().rename(columns=lambda x: x.removeprefix("fraction__")).stack().rename("recomputed").reset_index()
    recomputed.columns = ["sample", "state", "recomputed"]
    joined = state_summary.merge(recomputed, on=["sample", "state"], validate="one_to_one")
    summary_ok = state_set_ok and np.allclose(joined["mean_fraction"], joined["recomputed"], atol=2e-10)
    check("state_summary_recomputed", bool(summary_ok), f"rows={len(state_summary)}")

    niche_prefix = np.array([str(n).startswith(f"{s}-N") for s, n in zip(wide["sample"], wide["niche"])])
    assignment_counts = wide.groupby(["sample", "niche"]).size().rename("recomputed_n").reset_index()
    niche_join = niche_summary.merge(assignment_counts, on=["sample", "niche"], validate="one_to_one")
    niche_ok = niche_prefix.all() and len(niche_join) == len(niche_summary) and (niche_join["n_spots"] == niche_join["recomputed_n"]).all()
    check("niche_assignment_conservation", bool(niche_ok), f"niches={len(niche_summary)}")

    selected = k_table[k_table["selected"] == True].copy()
    selection_ok = len(selected) == 14 and selected["valid"].eq(True).all()
    for sample, group in k_table.groupby("sample"):
        valid = group[group["valid"] == True].sort_values(["silhouette", "k"], ascending=[False, True])
        expected_k = int(valid.iloc[0]["k"])
        actual_k = int(selected.loc[selected["sample"] == sample, "k"].iloc[0])
        selection_ok &= expected_k == actual_k and wide.loc[wide["sample"] == sample, "niche"].nunique() == actual_k
    check("independent_k_selection", bool(selection_ok), "per-section best valid silhouette, K=3..8")

    comp_ok = len(niche_comp) == len(niche_summary) * 13 and set(niche_comp["state"]) == set(EXPECTED_STATES)
    comp_sum = niche_comp.groupby(["sample", "niche"])["mean_neighborhood_fraction"].sum()
    comp_ok &= np.allclose(comp_sum.to_numpy(), 1.0, atol=2e-7)
    check("niche_composition_conservation", bool(comp_ok), f"rows={len(niche_comp)}")
    adjacency_ok = (adjacency["observed_edges"] >= 0).all() and (adjacency["expected_edges_independent_labels"] > 0).all() and np.isfinite(adjacency["observed_over_expected"]).all() and adjacency["interpretation"].eq("descriptive_not_inferential").all()
    check("adjacency_descriptive_boundary", bool(adjacency_ok), f"rows={len(adjacency)}")
    check("cohort_separation", set(wide.loc[wide["cohort"] == "GLMF", "sample"]) == {"GLMF1", "GLMF2"} and wide.loc[wide["cohort"] == "TJRCC", "sample"].nunique() == 12,
          "GLMF=2; TJRCC=12; no pooled clustering")

    checks_frame = pd.DataFrame(checks)
    passed = checks_frame["status"].eq("PASS").all()
    atomic_tsv(checks_frame, output / "PHASE18_INDEPENDENT_VALIDATION_STATUS_v1.tsv")
    validation = {
        "schema_version": "phase18_composition_niche_independent_validation_v1",
        "status": "PASS" if passed else "FAIL",
        "validated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "checks_passed": int(checks_frame["status"].eq("PASS").sum()),
        "checks_total": len(checks_frame), "n_sections": int(wide["sample"].nunique()),
        "n_spots": len(wide), "n_niches": len(niche_summary),
    }
    atomic_json(output / "PHASE18_INDEPENDENT_VALIDATION_SUMMARY_v1.json", validation)
    print(json.dumps(validation, indent=2))
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
