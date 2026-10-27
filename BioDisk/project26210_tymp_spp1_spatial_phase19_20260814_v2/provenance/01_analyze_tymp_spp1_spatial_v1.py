#!/usr/bin/env python3
from __future__ import annotations

import argparse, gzip, hashlib, json, math, os, platform
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy
from scipy.sparse import issparse
from scipy.spatial import cKDTree
from scipy.stats import binomtest, rankdata, spearmanr, wilcoxon

SEED = 20260814
BOOTSTRAPS = 9999
EXPECTED_EPOCH = 10000
EXPECTED_STATES = ["B cell", "C1QC-state Macrophage", "Dendritic cell", "Endothelial", "Epithelial-Normal", "Malignant", "Mast cell", "Monocyte", "NK cell", "Other Macrophage", "SMC/Pericyte", "SPP1-state Macrophage", "T cell"]
SIGNATURES = ["tymp_expression", "spp1_state_no_spp1", "spp1_state_full", "c1qc_state", "polyamine_mta_context", "f1_top100_weighted"]
PRIMARY = ["neighbor_TYMP_to_SPP1_vs_C1QC_expression", "neighbor_TYMP_to_SPP1_component_polarization", "neighbor_TYMP_to_SPP1_state_no_SPP1"]

def now(): return datetime.now().astimezone().isoformat(timespec="seconds")
def sha256(path):
    d=hashlib.sha256()
    with Path(path).open("rb") as h:
        for b in iter(lambda:h.read(8*1024*1024),b""): d.update(b)
    return d.hexdigest()
def atomic_json(path,payload):
    part=Path(str(path)+".part"); part.write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8"); os.replace(part,path)
def atomic_tsv(frame,path):
    part=Path(str(path)+".part")
    if str(path).endswith(".gz"):
        with gzip.open(part,"wt",encoding="utf-8",newline="") as h: frame.to_csv(h,sep="\t",index=False,lineterminator="\n")
    else: frame.to_csv(part,sep="\t",index=False,lineterminator="\n")
    os.replace(part,path)
def zscore(x):
    x=np.asarray(x,float); s=np.nanstd(x); return np.zeros(len(x)) if not np.isfinite(s) or s<=0 else (x-np.nanmean(x))/s
def residualize(y,covs):
    x=np.column_stack([np.ones(len(y)),*[zscore(c) for c in covs]]); target=zscore(y); return zscore(target-x@np.linalg.lstsq(x,target,rcond=None)[0])
def safe_rho(x,y):
    x,y=np.asarray(x,float),np.asarray(y,float); ok=np.isfinite(x)&np.isfinite(y)
    return math.nan if ok.sum()<5 or np.ptp(x[ok])==0 or np.ptp(y[ok])==0 else float(spearmanr(x[ok],y[ok]).statistic)
def bh(p):
    p=np.asarray(p,float); order=np.argsort(p); q=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]; out=np.empty(len(p)); out[order]=np.minimum(q,1); return out
def clean_state(x):
    x=str(x); prefix="q05cell_abundance_w_sf_"; return x[len(prefix):] if x.startswith(prefix) else x

def discover(root):
    rows=[]
    for cohort,sub in [("GLMF","glmf12_cell2location_tymp_continuous_v1"),("TJRCC","tjrcc12_cell2location_tymp_continuous_v1")]:
        for status_path in sorted((root/sub).glob("*/status.json")):
            s=json.loads(status_path.read_text(encoding="utf-8")); sample=s["sample"]; h=Path(s["output_h5ad"])
            if s["status"]!="completed" or int(s["epoch"])!=EXPECTED_EPOCH or h.stat().st_size!=int(s["output_size_bytes"]): raise RuntimeError(f"Phase17 gate failed {sample}")
            actual_hash=sha256(h)
            if actual_hash!=s["output_sha256"]: raise RuntimeError(f"Phase17 H5AD SHA256 mismatch {sample}")
            rows.append({"cohort":cohort,"sample":sample,"h5ad":h,"status":s,"status_path":status_path,"h5ad_sha256":actual_hash})
    if len(rows)!=14 or len({r["sample"] for r in rows})!=14: raise RuntimeError("Expected 14 unique mappings")
    return rows

def block_labels(coords,grid=4):
    bins=[]
    for j in range(2): bins.append(pd.qcut(pd.Series(coords[:,j]).rank(method="first"),q=grid,labels=False).to_numpy())
    return bins[0]*grid+bins[1]
def block_rho(x,y,labels,mask):
    f=pd.DataFrame({"block":labels[mask],"x":np.asarray(x)[mask],"y":np.asarray(y)[mask]}).groupby("block",as_index=False).agg(x=("x","median"),y=("y","median"),n=("x","size"))
    return safe_rho(f.x,f.y),len(f),int(mask.sum())

def load_scores(rec,genes):
    obj=ad.read_h5ad(rec["h5ad"])
    try:
        wanted=genes[genes.signature_id.isin(SIGNATURES)].gene.drop_duplicates().tolist(); present=[g for g in wanted if g in obj.var_names]
        idx=[obj.var_names.get_loc(g) for g in present]; matrix=obj.X[:,idx]
        if issparse(matrix): matrix=matrix.toarray()
        matrix=np.asarray(matrix,float); totals=np.asarray(obj.X.sum(axis=1)).ravel().astype(float)
        if np.any(totals<=0): raise RuntimeError("nonpositive library")
        lognorm=np.log1p(matrix*(10000/totals[:,None])); pos={g:i for i,g in enumerate(present)}
        abundance=obj.obsm["q05_cell_abundance_w_sf"]; av=np.asarray(abundance,float)
        states=[clean_state(x) for x in getattr(abundance,"columns",[])]
        if len(states)!=av.shape[1]: states=list(map(str,obj.uns["mod"]["factor_names"]))
        if states!=EXPECTED_STATES or av.shape!=(obj.n_obs,13) or np.any(~np.isfinite(av)) or np.any(av<0): raise RuntimeError(f"abundance gate {rec['sample']}")
        coords=np.asarray(obj.obsm["spatial"],float)[:,:2]; spots=obj.obs_names.astype(str).tolist()
    finally: del obj
    frame=pd.DataFrame({"cohort":rec["cohort"],"sample":rec["sample"],"spot_id":spots,"x":coords[:,0],"y":coords[:,1],"total_counts":totals})
    qc=[]
    for sid in SIGNATURES:
        d=genes[genes.signature_id==sid].sort_values("gene_order"); used=[g for g in d.gene if g in pos]; minimum=int(d.minimum_present_genes.iloc[0]) if "minimum_present_genes" in d else (70 if sid=="f1_top100_weighted" else 1)
        if len(used)<minimum: raise RuntimeError(f"signature coverage {rec['sample']} {sid} {len(used)}/{minimum}")
        columns=[pos[g] for g in used]
        if sid=="f1_top100_weighted":
            w=d.set_index("gene").loc[used,"raw_weight"].astype(float).to_numpy(); primary=lognorm[:,columns]@w/w.sum()
        elif sid=="tymp_expression": primary=lognorm[:,columns[0]]
        else: primary=np.column_stack([zscore(lognorm[:,j]) for j in columns]).mean(axis=1)
        frame[sid]=primary; frame["sensitivity__"+sid]=lognorm[:,columns].mean(axis=1)
        qc.append({"cohort":rec["cohort"],"sample":rec["sample"],"signature_id":sid,"defined_genes":len(d),"present_genes":len(used),"missing_genes":";".join([g for g in d.gene if g not in pos])})
    for j,state in enumerate(states): frame[state]=av[:,j]
    return frame,pd.DataFrame(qc)

def main():
    p=argparse.ArgumentParser(); p.add_argument("--phase17-root",type=Path,required=True); p.add_argument("--phase18-root",type=Path,required=True); p.add_argument("--phase18-validation",type=Path,required=True); p.add_argument("--signature-genes",type=Path,required=True); p.add_argument("--signature-registry",type=Path,required=True); p.add_argument("--analysis-registry",type=Path,required=True); p.add_argument("--output",type=Path,required=True); a=p.parse_args()
    if a.output.exists(): raise FileExistsError(a.output)
    a.output.mkdir(parents=True); tables=a.output/"tables"; tables.mkdir(); atomic_json(a.output/"RUNNING.json",{"status":"RUNNING","started_at":now()})
    phase18_status=json.loads((a.phase18_root/"PHASE18_COMPOSITION_NICHE_STATUS_v1.json").read_text());
    if phase18_status.get("status")!="PASS_ANALYSIS_COMPLETE": raise RuntimeError("Phase18 status gate")
    phase18_validation=json.loads((a.phase18_validation/"PHASE18_INDEPENDENT_VALIDATION_SUMMARY_v1.json").read_text())
    if phase18_validation.get("status")!="PASS": raise RuntimeError("Phase18 independent validation gate")
    wide18=pd.read_csv(a.phase18_root/"tables/spot_composition_wide.tsv.gz",sep="\t",usecols=["sample","spot_id","niche"])
    registry=pd.read_csv(a.signature_registry,sep="\t"); genes=pd.read_csv(a.signature_genes,sep="\t").merge(registry[["signature_id","minimum_present_genes"]],on="signature_id",how="left")
    records=discover(a.phase17_root); all_frames=[]; qcs=[]; stats=[]; niche_rows=[]; provenance=[]
    for i,rec in enumerate(records):
        print(f"[{now()}] {rec['cohort']} {rec['sample']}",flush=True); frame,qc=load_scores(rec,genes); qcs.append(qc)
        join=wide18[wide18["sample"]==rec["sample"]][["spot_id","niche"]]
        frame=frame.merge(join,on="spot_id",how="left",validate="one_to_one")
        if frame.niche.isna().any(): raise RuntimeError(f"Phase18 identity join {rec['sample']}")
        macro=frame[["C1QC-state Macrophage","SPP1-state Macrophage","Other Macrophage"]].sum(axis=1).to_numpy(); malignant=frame["Malignant"].to_numpy(); depth=np.log1p(frame.total_counts.to_numpy())
        tymp=residualize(frame.tymp_expression,[depth,np.log1p(malignant),frame.f1_top100_weighted]); tymp_sens=residualize(frame.tymp_expression,[depth,np.log1p(malignant)])
        expr_contrast=residualize(frame.spp1_state_no_spp1-frame.c1qc_state,[depth,np.log1p(macro)]); spp1_resid=residualize(frame.spp1_state_no_spp1,[depth,np.log1p(macro)])
        component=frame["SPP1-state Macrophage"].to_numpy()/np.maximum(macro,1e-12); mta=residualize(frame.polyamine_mta_context,[depth,np.log1p(malignant),frame.f1_top100_weighted])
        neighbors=cKDTree(frame[["x","y"]]).query(frame[["x","y"]],k=7)[1][:,1:]; ntymp=tymp[neighbors].mean(axis=1); ntymp_sens=tymp_sens[neighbors].mean(axis=1); nmta=mta[neighbors].mean(axis=1); labels=block_labels(frame[["x","y"]].to_numpy())
        variables={"tymp_exposure":tymp,"neighbor_tymp_exposure":ntymp,"neighbor_tymp_no_f1_sensitivity":ntymp_sens,"neighbor_mta_context":nmta,"spp1_vs_c1qc_expression":expr_contrast,"spp1_state_no_spp1_residual":spp1_resid,"spp1_component_polarization":component}
        for name,val in variables.items(): frame[name]=val
        relations=[("neighbor_TYMP_to_SPP1_vs_C1QC_expression",ntymp,expr_contrast,"primary"),("neighbor_TYMP_to_SPP1_component_polarization",ntymp,component,"primary"),("neighbor_TYMP_to_SPP1_state_no_SPP1",ntymp,spp1_resid,"primary"),("same_spot_TYMP_to_SPP1_vs_C1QC_expression",tymp,expr_contrast,"sensitivity"),("neighbor_TYMP_no_F1_to_SPP1_vs_C1QC_expression",ntymp_sens,expr_contrast,"sensitivity"),("neighbor_MTA_context_to_SPP1_vs_C1QC_expression",nmta,expr_contrast,"context_proxy")]
        for gate,q in [("q25",.25),("q50",.5)]:
            mask=(malignant>=np.quantile(malignant,q))&(macro>=np.quantile(macro,q))
            for relation,x,y,role in relations:
                rho,nblocks,nspots=block_rho(x,y,labels,mask); stats.append({"cohort":rec["cohort"],"sample":rec["sample"],"relation":relation,"role":role,"gate":gate,"grid":"4x4","n_spots":nspots,"n_blocks":nblocks,"spearman_rho":rho})
        for niche,g in frame.groupby("niche"):
            niche_rows.append({"cohort":rec["cohort"],"sample":rec["sample"],"niche":niche,"n_spots":len(g),"median_neighbor_tymp":float(g.neighbor_tymp_exposure.median()),"median_spp1_vs_c1qc_expression":float(g.spp1_vs_c1qc_expression.median()),"median_spp1_component_polarization":float(g.spp1_component_polarization.median()),"median_neighbor_mta_context":float(g.neighbor_mta_context.median())})
        provenance.append({"cohort":rec["cohort"],"sample":rec["sample"],"h5ad":str(rec["h5ad"]),"size_bytes":rec["h5ad"].stat().st_size,"sha256_recomputed":rec["h5ad_sha256"],"sha256_recorded":rec["status"]["output_sha256"],"n_spots":len(frame)})
        all_frames.append(frame[["cohort","sample","spot_id","x","y","niche",*variables.keys()]])
    spot=pd.concat(all_frames,ignore_index=True); section=pd.DataFrame(stats); primary=section[(section.cohort=="TJRCC")&(section.role=="primary")&(section.gate=="q25")].copy()
    if len(primary)!=36 or primary.spearman_rho.isna().any(): raise RuntimeError("primary conservation")
    rng=np.random.default_rng(SEED); agg=[]; loo=[]
    for relation in PRIMARY:
        v=primary.loc[primary.relation==relation,"spearman_rho"].to_numpy(); boot=np.median(v[rng.integers(0,len(v),size=(BOOTSTRAPS,len(v)))],axis=1); nonzero=v[v!=0]
        try: wp=float(wilcoxon(v,alternative="two-sided").pvalue)
        except ValueError: wp=1.0
        sp=float(binomtest(int((nonzero>0).sum()),len(nonzero),.5).pvalue) if len(nonzero) else 1.0
        lv=[]
        for sample in sorted(primary["sample"].unique()):
            med=float(np.median(primary[(primary.relation==relation)&(primary["sample"]!=sample)].spearman_rho)); lv.append(med); loo.append({"relation":relation,"left_out_sample":sample,"median_rho":med})
        agg.append({"relation":relation,"sections":len(v),"median_rho":float(np.median(v)),"q25_rho":float(np.quantile(v,.25)),"q75_rho":float(np.quantile(v,.75)),"positive_sections":int((v>0).sum()),"negative_sections":int((v<0).sum()),"wilcoxon_p":wp,"sign_test_p":sp,"bootstrap_ci_low":float(np.quantile(boot,.025)),"bootstrap_ci_high":float(np.quantile(boot,.975)),"loo_positive_fraction":float((np.array(lv)>0).mean())})
    aggregate=pd.DataFrame(agg); aggregate["wilcoxon_fdr_bh"]=bh(aggregate.wilcoxon_p)
    aggregate["decision"]=np.where((aggregate.median_rho>=.1)&(aggregate.positive_sections>=8)&(aggregate.wilcoxon_fdr_bh<.05)&(aggregate.bootstrap_ci_low>0)&(aggregate.loo_positive_fraction==1),"PRIMARY_SPATIAL_SUPPORT",np.where(aggregate.positive_sections>=8,"DIRECTIONAL_NOT_FDR",np.where((aggregate.positive_sections>=4)&(aggregate.negative_sections>=4),"HETEROGENEOUS","NOT_SUPPORTED")))
    outputs={"spot_scores_v1.tsv.gz":spot,"section_statistics_v1.tsv":section,"primary_tjrcc_section_statistics_v1.tsv":primary,"primary_across_section_inference_v1.tsv":aggregate,"primary_leave_one_section_out_v1.tsv":pd.DataFrame(loo),"niche_tymp_spp1_summary_v1.tsv":pd.DataFrame(niche_rows),"signature_gene_coverage_v1.tsv":pd.concat(qcs,ignore_index=True),"input_provenance_v1.tsv":pd.DataFrame(provenance)}
    for name,f in outputs.items(): atomic_tsv(f,tables/name)
    atomic_json(a.output/"software_versions_v1.json",{"python":platform.python_version(),"anndata":ad.__version__,"numpy":np.__version__,"pandas":pd.__version__,"scipy":scipy.__version__,"script_sha256":sha256(Path(__file__))})
    manifest=[]
    for path in sorted(p for p in a.output.rglob("*") if p.is_file() and p.name!="RUNNING.json"):
        manifest.append({"file":path.relative_to(a.output).as_posix(),"size_bytes":path.stat().st_size,"sha256":sha256(path)})
    atomic_tsv(pd.DataFrame(manifest),a.output/"OUTPUT_MANIFEST_v1.tsv")
    status={"schema_version":"phase19_tymp_spp1_spatial_v2","status":"PASS_ANALYSIS_COMPLETE","completed_at":now(),"sections":14,"tjrcc_sections":12,"glmf_sections":2,"spots":len(spot),"primary_rows":len(primary),"primary_relations":PRIMARY,"section_as_biological_unit":True,"spot_p_values":False,"mta_interpretation":"transcriptomic context proxy only","analysis_registry":str(a.analysis_registry),"analysis_registry_sha256":sha256(a.analysis_registry),"signature_registry_sha256":sha256(a.signature_registry),"phase18_spot_table_sha256":sha256(a.phase18_root/"tables/spot_composition_wide.tsv.gz"),"phase18_validation_sha256":sha256(a.phase18_validation/"PHASE18_INDEPENDENT_VALIDATION_STATUS_v1.tsv"),"version_note":"v1 stopped before scoring because two old-script signature IDs were not in the Phase13 frozen registry; v2 uses only frozen signatures plus validated abundances"}
    atomic_json(a.output/"PHASE19_TYMP_SPP1_SPATIAL_STATUS_v1.json",status); (a.output/"RUNNING.json").unlink(); print(json.dumps(status,indent=2))

if __name__=="__main__": main()
