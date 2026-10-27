#!/usr/bin/env python3
"""Audit gene overlap for every Phase 53 bridge relation."""

from __future__ import annotations

import importlib.util
import json
import os
from datetime import datetime
from pathlib import Path

import pandas as pd


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SOURCE_SCRIPT = HERE / "01_screen_reference_bridge_v1.py"
GENES = ROOT / "phase13_signature_lock_local_v1" / "registry_v1" / "frozen_signature_genes_v1.tsv"
REPORTABLE = HERE / "analysis_local_v2" / "bridge_screen_reportable_exact_v2.tsv"
FINAL = HERE / "signature_overlap_audit_local_v2"
PART = HERE / "signature_overlap_audit_local_v2.part"


def atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
    part = path.with_name(path.name + ".part")
    frame.to_csv(part, sep="\t", index=False, na_rep="NA", encoding="utf-8", lineterminator="\n")
    os.replace(part, path)


def main() -> int:
    if FINAL.exists() or PART.exists():
        raise FileExistsError("Refusing existing overlap-audit output")
    PART.mkdir(parents=True)
    spec = importlib.util.spec_from_file_location("bridge_v1", SOURCE_SCRIPT)
    bridge = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(bridge)
    genes = pd.read_csv(GENES, sep="\t")
    gene_sets = genes.groupby("signature_id")["gene"].apply(lambda x: set(x.astype(str))).to_dict()
    rows = []
    for relation, edge_group, x_sig, x_comp, y_sig, y_comp, expected in bridge.RELATIONS:
        if x_sig not in gene_sets or y_sig not in gene_sets:
            raise RuntimeError(f"Missing signature genes for {relation}")
        x = gene_sets[x_sig]
        y = gene_sets[y_sig]
        overlap = sorted(x & y)
        union = x | y
        rows.append({
            "relation": relation,
            "edge_group": edge_group,
            "x_signature": x_sig,
            "y_signature": y_sig,
            "x_gene_count": len(x),
            "y_gene_count": len(y),
            "overlap_count": len(overlap),
            "overlap_genes": ";".join(overlap) or "none",
            "jaccard": len(overlap) / len(union),
            "overlap_flag": "PARTIAL_OVERLAP_CAUTION" if overlap else "INDEPENDENT_GENE_SETS",
            "cross_compartment": x_comp != y_comp,
            "expected_direction": expected,
        })
    audit = pd.DataFrame(rows)
    reportable = pd.read_csv(REPORTABLE, sep="\t")
    audited = reportable.merge(
        audit[["relation", "overlap_count", "overlap_genes", "jaccard", "overlap_flag", "cross_compartment"]],
        on="relation",
        how="left",
        validate="many_to_one",
    )
    if len(audit) != len(bridge.RELATIONS) or audited["overlap_flag"].isna().any():
        raise RuntimeError("Signature overlap validation failed")
    atomic_tsv(audit, PART / "bridge_relation_signature_overlap_v1.tsv")
    atomic_tsv(audited, PART / "bridge_reportable_overlap_audited_v1.tsv")
    status = {
        "status": "PASS",
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "relations": len(audit),
        "relations_with_overlap": int((audit["overlap_count"] > 0).sum()),
        "reportable_rows": len(audited),
        "reportable_rows_with_overlap": int((audited["overlap_count"] > 0).sum()),
    }
    (PART / "SIGNATURE_OVERLAP_AUDIT_STATUS_v1.json").write_text(
        json.dumps(status, indent=2) + "\n", encoding="utf-8"
    )
    (PART / "SIGNATURE_OVERLAP_AUDIT_PASS_v1").write_text(
        f"status=PASS\nrelations={len(audit)}\nrelations_with_overlap={status['relations_with_overlap']}\n",
        encoding="utf-8",
    )
    os.replace(PART, FINAL)
    print(f"PASS relations={len(audit)} overlap={status['relations_with_overlap']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
