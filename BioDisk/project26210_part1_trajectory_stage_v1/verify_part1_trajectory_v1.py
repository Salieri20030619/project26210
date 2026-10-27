#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
OUT = PROJECT / "results/part1_rebuild_v1/09_metabolism_trajectory/trajectory_v1"
GATE = PROJECT / "check_memory_gate.sh"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""): h.update(block)
    return h.hexdigest()


def main() -> None:
    checks = []
    def check(name: str, passed: bool, detail: str) -> None: checks.append({"check": name, "status": "PASS" if passed else "FAIL", "detail": detail})
    required_markers = [OUT / "input/INPUT_COMPLETE", OUT / "R_METHODS_COMPLETE", OUT / "PYTHON_METHODS_COMPLETE"]
    check("stage_markers", all(x.exists() for x in required_markers), ";".join(str(x) for x in required_markers))
    cells = pd.read_csv(OUT / "input/cells_v1.tsv", sep="\t", header=None, dtype=str).iloc[:,0]
    check("strict_cell_registry", len(cells)==11950 and cells.nunique()==11950, f"n={len(cells)};unique={cells.nunique()}")
    r = pd.read_csv(OUT / "source_data/trajectory_r_cell_scores_v1.tsv.gz", sep="\t", dtype={"exact_cell_id":str})
    p = pd.read_csv(OUT / "source_data/trajectory_all_methods_cell_scores_v1.tsv.gz", sep="\t", index_col=0)
    check("R_cell_order", r["exact_cell_id"].tolist()==cells.tolist(), f"n={len(r)}")
    check("Python_cell_set", set(p.index.astype(str))==set(cells), f"n={len(p)}")
    for col in ("Monocle3_pseudotime_v1", "Slingshot_pseudotime_v1"):
        finite = int(np.isfinite(pd.to_numeric(r[col],errors="coerce")).sum()); check(col+"_finite", finite>=int(.95*len(cells)), f"finite={finite}")
    for col in ("CytoTRACE2_Score", "scTour_pseudotime_patient_v1", "CellRank_max_fate_probability_v1"):
        finite = int(np.isfinite(pd.to_numeric(p[col],errors="coerce")).sum()); check(col+"_finite", finite>=int(.95*len(cells)), f"finite={finite}")
    graph = pd.read_csv(OUT / "tables/monocle3_graph_test_hvg_v1.tsv.gz", sep="\t")
    check("graph_test_rows", len(graph)==2500 and "q_value" in graph, f"rows={len(graph)}")
    summary = pd.read_csv(OUT / "tables/trajectory_all_methods_correlation_summary_v1.tsv", sep="\t")
    check("patient_summary", len(summary)==48 and summary["n_patients"].between(1,26).all(), f"rows={len(summary)}")
    h5 = ad.read_h5ad(OUT / "objects/part1_malignant_trajectory_all_methods_v1.h5ad", backed="r")
    check("trajectory_h5ad_shape", h5.shape==(11950,2500), f"shape={h5.shape}"); h5.file.close()
    figures = list((OUT / "figures").glob("*.png")) + list((OUT / "figures").glob("*.pdf"))
    check("figure_files", len(figures)>=10 and all(x.stat().st_size>1000 for x in figures), f"n={len(figures)}")
    table = pd.DataFrame(checks); table.to_csv(OUT / "qc/trajectory_verification_checks_v1.tsv",sep="\t",index=False)
    failures = table.loc[table.status!="PASS","check"].tolist()
    metadata = {"schema_version":"part1_trajectory_verification_v1","completed_at":datetime.now().astimezone().isoformat(timespec="seconds"),"n_checks":len(checks),"n_failed":len(failures),"failed_checks":failures,"manual_figure_review_pending":True}
    (OUT / "qc/trajectory_verification_metadata_v1.json").write_text(json.dumps(metadata,indent=2),encoding="utf-8")
    if failures: raise SystemExit("Verification failed: "+",".join(failures))
    files = [x for x in OUT.rglob("*") if x.is_file() and x.name not in {"RUN_COMPLETE","file_manifest_v1.tsv"} and "logs" not in x.parts]
    manifest = pd.DataFrame([{"file":str(x.relative_to(OUT)),"size_bytes":x.stat().st_size,"sha256":sha256(x)} for x in sorted(files)])
    manifest.to_csv(OUT / "file_manifest_v1.tsv",sep="\t",index=False)
    (OUT / "RUN_COMPLETE").write_text(datetime.now().astimezone().isoformat(timespec="seconds")+"\n",encoding="utf-8")
    print(json.dumps(metadata,indent=2))


if __name__ == "__main__": main()
