#!/usr/bin/env python3
"""Independent gate for the 13-state continuous-TYMP spatial reference."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd


EXPECTED = sorted([
    "B cell", "C1QC-state Macrophage", "Dendritic cell", "Endothelial",
    "Epithelial-Normal", "Malignant", "Mast cell", "Monocyte", "NK cell",
    "Other Macrophage", "SMC/Pericyte", "SPP1-state Macrophage", "T cell",
])


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--reference", type=Path, required=True)
    p.add_argument("--old-signature", type=Path, required=True)
    p.add_argument("--bridge", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError(f"Refusing to overwrite: {a.output}")
    work = a.output.with_name(a.output.name + ".part")
    if work.exists():
        raise FileExistsError(f"Stale work directory: {work}")
    work.mkdir(parents=True)
    sig_path = a.reference / "reference_tymp_continuous_v1_signatures_study_balanced_mean_counts.tsv.gz"
    assn_path = a.reference / "reference_tymp_continuous_v1_cell_assignments.tsv.gz"
    qc_path = a.reference / "reference_tymp_continuous_v1_state_qc.tsv"
    manifest_path = a.reference / "file_manifest_tymp_continuous_v1.tsv"
    needed = [sig_path, assn_path, qc_path, manifest_path]
    if not all(path.is_file() for path in needed):
        raise RuntimeError("Missing continuous-TYMP reference outputs")

    rows: list[dict[str, object]] = []
    def add(name: str, passed: bool, observed: object, expected: object) -> None:
        rows.append({"check": name, "pass": bool(passed), "observed": observed, "expected": expected})

    manifest = pd.read_csv(manifest_path, sep="\t")
    hashes_ok = all(
        (a.reference / r.file).is_file()
        and (a.reference / r.file).stat().st_size == int(r.size_bytes)
        and digest(a.reference / r.file) == r.sha256
        for r in manifest.itertuples(index=False)
    )
    add("manifest_hashes", hashes_ok, int(hashes_ok), 1)
    signature = pd.read_csv(sig_path, sep="\t", index_col=0)
    add("signature_shape", signature.shape == (12758, 13), f"{signature.shape[0]}x{signature.shape[1]}", "12758x13")
    add("state_names", sorted(signature.columns) == EXPECTED, signature.shape[1], 13)
    values = signature.to_numpy(dtype=float)
    add("signature_finite_nonnegative", np.isfinite(values).all() and (values >= 0).all(), float(values.min()), ">=0 finite")
    old = pd.read_csv(a.old_signature, sep="\t", index_col=0)
    nonmalignant = sorted(set(EXPECTED).difference({"Malignant"}))
    delta = float(np.abs(signature[nonmalignant].to_numpy() - old[nonmalignant].to_numpy()).max())
    add("nonmalignant_signatures_unchanged", delta <= 1e-12, delta, "<=1e-12")

    assignments = pd.read_csv(assn_path, sep="\t", index_col="exact_cell_id")
    add("assignment_ids_unique", not assignments.index.has_duplicates, int(assignments.index.duplicated().sum()), 0)
    counts = assignments["state_tymp_continuous_v1"].value_counts()
    qc = pd.read_csv(qc_path, sep="\t").set_index("state_tymp_continuous_v1")
    conserved = counts.sort_index().equals(qc["n_cells"].astype(int).sort_index())
    add("state_counts_conserved", conserved, int((counts.sort_index() - qc["n_cells"].astype(int).sort_index()).abs().sum()), 0)
    malignant = assignments[assignments["state_tymp_continuous_v1"].eq("Malignant")]
    bridge = pd.read_csv(a.bridge, sep="\t", index_col="exact_cell_id")
    add("malignant_cells_match_bridge", set(malignant.index) == set(bridge.index), len(malignant), len(bridge))
    add("malignant_coverage", len(malignant) == 11921 and malignant["patient"].nunique() == 20 and malignant["study_v1"].nunique() == 5,
        f"{len(malignant)} cells/{malignant['patient'].nunique()} patients/{malignant['study_v1'].nunique()} studies",
        "11921/20/5")
    bridge_match = np.allclose(
        malignant["TYMP_log1p_CP10K"].sort_index().to_numpy(),
        bridge["TYMP_log1p_CP10K"].sort_index().to_numpy(), equal_nan=True,
    )
    add("continuous_tymp_preserved", bridge_match, int(bridge_match), 1)
    f1_match = np.allclose(malignant["F1"].sort_index().to_numpy(), bridge["F1"].sort_index().to_numpy())
    add("continuous_f1_preserved", f1_match, int(f1_match), 1)
    add("single_malignant_state_only", not any("TYMP-high" in x or "TYMP-low" in x or "F1-high" in x or "F1-low" in x for x in signature.columns),
        "|".join(signature.columns), "one Malignant state")

    table = pd.DataFrame(rows)
    passed = bool(table["pass"].all())
    table.to_csv(work / "reference_tymp_continuous_v1_independent_checks.tsv", sep="\t", index=False)
    status = {
        "schema_version": "phase17_reference_tymp_continuous_independent_validation_v1",
        "status": "PASS" if passed else "FAIL", "completed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "n_checks": len(table), "n_pass": int(table["pass"].sum()), "n_fail": int((~table["pass"]).sum()),
        "reference_signature_sha256": digest(sig_path), "reference_manifest_sha256": digest(manifest_path),
        "design": "one identifiable Malignant Cell2location state plus continuous TYMP and F1 annotations",
        "rejected_design": "TYMP-high/low Cell2location subtypes; Spearman rho 0.9894368522597528",
    }
    part = work / "PHASE17_REFERENCE_TYMP_CONTINUOUS_INDEPENDENT_VALIDATION_v1.json.part"
    part.write_text(json.dumps(status, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(part, work / "PHASE17_REFERENCE_TYMP_CONTINUOUS_INDEPENDENT_VALIDATION_v1.json")
    os.replace(work, a.output)
    print(json.dumps(status, ensure_ascii=False, indent=2))
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
