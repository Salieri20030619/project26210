#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import subprocess
from datetime import datetime
from pathlib import Path


ROOT = Path("/mnt/BioDisk/wangmingyuan/project26210_envs/part1_trajectory_v1")
R = ROOT / "bin/Rscript"
PY = ROOT / "bin/python"
SOURCE_R = Path("/mnt/BioDisk/ludingyu/project26210/envs/part2_6_2_7_20260726_v1/R_library")
SOURCE_PY = Path("/mnt/BioDisk/ludingyu/project26210/envs/part2_6_2_7_20260726_v1/python/bin/python")


def execute(name: str, command: list[str], checks: list[dict[str, str]]) -> str:
    result = subprocess.run(command, text=True, capture_output=True, timeout=300)
    detail = (result.stdout + "\n" + result.stderr).strip()[-4000:]
    checks.append({"check": name, "status": "PASS" if result.returncode == 0 else "FAIL", "detail": detail.replace("\t", " ").replace("\n", " | ")})
    return detail


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    checks: list[dict[str, str]] = []
    for name, path in (("R_wrapper_executable", R), ("Python_wrapper_executable", PY), ("source_R_library", SOURCE_R), ("source_Python", SOURCE_PY)):
        passed = path.exists() and (path.is_dir() or path.stat().st_size > 0)
        checks.append({"check": name, "status": "PASS" if passed else "FAIL", "detail": str(path)})

    r_code = r'''
suppressPackageStartupMessages({library(monocle3);library(slingshot);library(SingleCellExperiment);library(Matrix)})
set.seed(26210)
x<-Matrix(matrix(rpois(500*80,2),500,80),sparse=TRUE)
rownames(x)<-paste0("G",seq_len(nrow(x)));colnames(x)<-paste0("C",seq_len(ncol(x)))
cds<-new_cell_data_set(x,cell_metadata=data.frame(row.names=colnames(x)),gene_metadata=data.frame(row.names=rownames(x),gene_short_name=rownames(x)))
cds<-preprocess_cds(cds,num_dim=10,method="PCA",norm_method="log")
cds<-reduce_dimension(cds,reduction_method="UMAP",preprocess_method="PCA",umap.metric="cosine",umap.n_neighbors=10,cores=1,verbose=FALSE)
cds<-cluster_cells(cds,reduction_method="UMAP",k=10,verbose=FALSE)
cds<-learn_graph(cds,use_partition=TRUE,close_loop=FALSE,learn_graph_control=list(ncenter=20),verbose=FALSE)
stopifnot(ncol(cds)==80,nrow(principal_graph_aux(cds)[["UMAP"]]$dp_mst)>0)
n<-120;t<-seq(0,1,length.out=n);xy<-cbind(t+rnorm(n,0,.03),sin(pi*t)+rnorm(n,0,.03));cl<-factor(rep(c("Base","Middle","Late"),each=40),levels=c("Base","Middle","Late"))
sce<-SingleCellExperiment(assays=list(counts=matrix(rpois(30*n,2),30,n)));reducedDims(sce)$PCA<-xy;colData(sce)$cluster<-cl
sce<-slingshot(sce,clusterLabels="cluster",reducedDim="PCA",start.clus="Base")
pt<-slingPseudotime(sce);stopifnot(nrow(pt)==n,sum(is.finite(pt))>=100)
cat("monocle3=",as.character(packageVersion("monocle3")),";slingshot=",as.character(packageVersion("slingshot")),";MONOCLE3_SLINGSHOT_MINIMAL_PASS\n",sep="")
'''
    execute("R_packages_and_minimal_trajectory", [str(R), "-e", r_code], checks)

    py_code = r'''
from importlib.metadata import version
from pathlib import Path
import numpy as np, anndata as ad, scanpy as sc, cellrank as cr
import sctour, cytotrace2_py, scvelo, torch
np.random.seed(26210)
a=ad.AnnData(np.random.poisson(2,(120,60)).astype("float32"))
sc.pp.normalize_total(a);sc.pp.log1p(a);sc.pp.pca(a,n_comps=15);sc.pp.neighbors(a,n_neighbors=10,n_pcs=15)
a.obs["pt"]=np.linspace(0,1,a.n_obs)
k=cr.kernels.PseudotimeKernel(a,time_key="pt").compute_transition_matrix(n_jobs=2,show_progress_bar=False)
c=cr.kernels.ConnectivityKernel(a).compute_transition_matrix();combo=0.8*k+0.2*c
assert combo.transition_matrix.shape==(120,120) and combo.transition_matrix.nnz>0
models=list((Path(cytotrace2_py.__file__).resolve().parent/"resources/models").glob("*.pt"));assert len(models)>=10
names=["scanpy","anndata","numpy","pandas","scipy","torch","sctour","cellrank","scvelo","cytotrace2-py"]
print(";".join(f"{x}={version(x)}" for x in names)+f";cytotrace_models={len(models)};CELLRANK_MINIMAL_PASS")
'''
    execute("Python_packages_models_and_CellRank_minimal", [str(PY), "-W", "ignore", "-c", py_code], checks)

    ROOT.mkdir(parents=True, exist_ok=True)
    with (ROOT / "environment_verification_v1.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["check", "status", "detail"], delimiter="\t")
        writer.writeheader(); writer.writerows(checks)
    failures = [row["check"] for row in checks if row["status"] != "PASS"]
    metadata = {
        "schema_version": "part1_trajectory_environment_v1",
        "completed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "R_executable": "/mnt/BioDisk/ludingyu/miniconda3/envs/r453/bin/Rscript",
        "R_library": str(SOURCE_R),
        "Python_executable": str(SOURCE_PY),
        "dependency_changes": False,
        "environment_strategy": "versioned Part1 wrappers over previously installed project-private trajectory runtimes",
        "shared_environment_modified": False,
        "n_checks": len(checks),
        "n_failed": len(failures),
        "failed_checks": failures,
    }
    (ROOT / "run_metadata_v1.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")
    files = [R, PY, ROOT / "verify_part1_trajectory_environment_v1.py", ROOT / "memory_gate_v1.tsv", ROOT / "environment_verification_v1.tsv", ROOT / "run_metadata_v1.json"]
    with (ROOT / "file_manifest_v1.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t"); writer.writerow(["file", "size_bytes", "sha256"])
        for path in files:
            writer.writerow([str(path.relative_to(ROOT)), path.stat().st_size, sha256(path)])
    if failures:
        raise SystemExit("Trajectory environment verification failed: " + ",".join(failures))
    (ROOT / "RUN_COMPLETE").write_text("PASS\n", encoding="utf-8")
    print(json.dumps(metadata, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
