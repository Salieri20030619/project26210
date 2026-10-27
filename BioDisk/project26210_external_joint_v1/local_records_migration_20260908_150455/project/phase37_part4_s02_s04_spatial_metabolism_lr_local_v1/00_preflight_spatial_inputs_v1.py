#!/usr/bin/env python3
from __future__ import annotations
import argparse, hashlib, json, sys
from datetime import datetime
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

FOUR = ["Pyrimidine metabolism", "Cysteine and methionine metabolism",
        "Arginine and proline metabolism", "Glycerophospholipid metabolism"]
LR_GENES = ["MST1", "MST1R", "ICAM5", "ITGAL", "ITGB2", "ITGAM", "CD209", "PRTG", "IGDCC4"]

def sha(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(8*1024*1024),b""): h.update(b)
    return h.hexdigest()
def read_gmt(path):
    out={}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        z=line.rstrip().split("\t")
        if len(z)>=3: out[z[0]]=list(dict.fromkeys(x for x in z[2:] if x))
    return out
def dense(x):
    if hasattr(x,"to_memory"): x=x.to_memory()
    if sparse.issparse(x): x=x.toarray()
    return np.asarray(x)
def atomic_tsv(df,p):
    q=Path(str(p)+".part"); df.to_csv(q,sep="\t",index=False); q.replace(p)
def atomic_json(x,p):
    q=Path(str(p)+".part"); q.write_text(json.dumps(x,indent=2),encoding="utf-8"); q.replace(p)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--phase18-provenance",type=Path,required=True); ap.add_argument("--phase18-spots",type=Path,required=True); ap.add_argument("--phase19-spots",type=Path,required=True); ap.add_argument("--custom-gmt",type=Path,required=True); ap.add_argument("--kegg-gmt",type=Path,required=True); ap.add_argument("--pathology-dir",type=Path,required=True); ap.add_argument("--output",type=Path,required=True); a=ap.parse_args()
    if a.output.exists(): raise FileExistsError(a.output)
    a.output.mkdir(parents=True)
    prov=pd.read_csv(a.phase18_provenance,sep="\t"); p18=pd.read_csv(a.phase18_spots,sep="\t"); p19=pd.read_csv(a.phase19_spots,sep="\t")
    custom=read_gmt(a.custom_gmt); kegg=read_gmt(a.kegg_gmt)
    if not all(x in custom and x in kegg for x in FOUR): raise RuntimeError("Missing one or more frozen pathway sets")
    rows=[]; cover=[]
    for r in prov.itertuples(index=False):
        hp=Path(r.input_h5ad); actual=sha(hp)
        h=ad.read_h5ad(hp,backed="r"); genes=set(map(str,h.var_names)); wanted=sorted(set(LR_GENES+sum([custom[x]+kegg[x] for x in FOUR],[])) & genes)
        idx=[h.var_names.get_loc(g) for g in wanted]; m=dense(h.X[:,idx])
        spatial_key="spatial" if "spatial" in h.obsm_keys() else None
        ids=set(map(str,h.obs_names)); x18=p18[p18["sample"].eq(r.sample)]; x19=p19[p19["sample"].eq(r.sample)]
        rows.append({"cohort":r.cohort,"sample":r.sample,"h5ad":str(hp),"recorded_sha256":r.input_h5ad_sha256,"recomputed_sha256":actual,"hash_match":actual==r.input_h5ad_sha256,"n_spots":h.n_obs,"n_genes":h.n_vars,"x_storage":type(h.X).__name__,"layers":";".join(h.layers.keys()),"spatial_key":spatial_key or "","phase18_spot_match":ids==set(x18.spot_id.astype(str)),"phase19_spot_match":ids==set(x19.spot_id.astype(str)),"candidate_values_finite":bool(np.isfinite(m).all()),"candidate_values_nonnegative":bool((m>=0).all()),"candidate_values_integer":bool(np.allclose(m,np.rint(m))),"n_lr_genes":sum(g in genes for g in LR_GENES),"all_lr_genes":all(g in genes for g in LR_GENES)})
        for source,sets in [("CUSTOM_PRIMARY",custom),("KEGG_SENSITIVITY",kegg)]:
            for pathway in FOUR:
                present=[g for g in sets[pathway] if g in genes]
                cover.append({"cohort":r.cohort,"sample":r.sample,"source":source,"pathway":pathway,"n_set_genes":len(sets[pathway]),"n_present":len(present),"coverage_fraction":len(present)/len(sets[pathway]),"present_genes":";".join(present)})
        h.file.close()
    rows=pd.DataFrame(rows); cover=pd.DataFrame(cover)
    path_rows=[]
    for sample in ["GLMF1","GLMF2"]:
        p=a.pathology_dir/f"{sample}_pathology_polygons.json"
        x=json.loads(p.read_text(encoding="utf-8")) if p.exists() else {}
        labels=[z.get("label") for z in x.get("polygons",[])]
        path_rows.append({"sample":sample,"file":str(p),"exists":p.exists(),"review_status":x.get("review_status",""),"annotation_source":x.get("annotation_source",""),"n_polygons":len(labels),"has_tumor":("tumor" in labels),"has_nat":("NAT" in labels),"formal_pathology_ready":x.get("review_status")=="pathologist_reviewed" and "tumor" in labels and "NAT" in labels})
    path_rows=pd.DataFrame(path_rows)
    atomic_tsv(rows,a.output/"spatial_h5ad_preflight_v1.tsv"); atomic_tsv(cover,a.output/"pathway_gene_coverage_v1.tsv"); atomic_tsv(path_rows,a.output/"boundary_annotation_preflight_v1.tsv")
    ok=(len(rows)==14 and rows.hash_match.all() and rows.phase18_spot_match.all() and rows.phase19_spot_match.all() and rows.candidate_values_finite.all() and rows.candidate_values_nonnegative.all() and rows.candidate_values_integer.all() and rows.all_lr_genes.all() and (cover.n_present>=10).all())
    payload={"status":"PASS_S02_S04_SPATIAL_PREFLIGHT" if ok else "FAIL_S02_S04_SPATIAL_PREFLIGHT","completed_at":datetime.now().astimezone().isoformat(timespec="seconds"),"n_sections":len(rows),"n_spots":int(rows.n_spots.sum()),"n_pathway_coverage_rows":len(cover),"all_hashes_match":bool(rows.hash_match.all()),"all_spot_ids_match":bool(rows.phase18_spot_match.all() and rows.phase19_spot_match.all()),"all_candidate_counts_raw_safe":bool(rows.candidate_values_finite.all() and rows.candidate_values_nonnegative.all() and rows.candidate_values_integer.all()),"formal_pathology_boundary_sections":int(path_rows.formal_pathology_ready.sum()),"draft_boundary_sections":int((path_rows.has_tumor & path_rows.has_nat).sum()),"figures_created":False}
    atomic_json(payload,a.output/"status.json"); print(json.dumps(payload,indent=2)); return 0 if ok else 1
if __name__=="__main__": sys.exit(main())
