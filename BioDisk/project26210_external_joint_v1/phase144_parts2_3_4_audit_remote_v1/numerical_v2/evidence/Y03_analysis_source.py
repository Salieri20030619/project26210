#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import os
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon


CANDIDATES = ["JUN", "JUNB", "CEBPD"]
ENDPOINTS = ["Monocle3_continuous", "DPT_continuous", "SPP1_minus_Monocyte", "SPP1_minus_C1QC"]


def rho(x: pd.Series, y: pd.Series) -> float:
    a,b=x.to_numpy(float),y.to_numpy(float); keep=np.isfinite(a)&np.isfinite(b)
    if keep.sum()<20 or np.ptp(a[keep])==0 or np.ptp(b[keep])==0:return math.nan
    return float(spearmanr(a[keep],b[keep]).statistic)


def bh(series: pd.Series) -> np.ndarray:
    p=series.to_numpy(float);out=np.full(len(p),np.nan);valid=np.isfinite(p)
    if not valid.any():return out
    v=p[valid];o=np.argsort(v);q=np.minimum.accumulate((v[o]*len(v)/np.arange(1,len(v)+1))[::-1])[::-1];r=np.empty(len(v));r[o]=np.minimum(q,1);out[valid]=r;return out


def signed_rank(values: np.ndarray) -> float:
    values=values[np.isfinite(values)]
    if not len(values) or np.all(np.abs(values)<1e-15):return 1.0
    return float(wilcoxon(values,alternative="two-sided",zero_method="wilcox",correction=False).pvalue)


def atomic_tsv(frame: pd.DataFrame,path:Path,compression:str|None=None)->None:
    part=Path(str(path)+".part");frame.to_csv(part,sep="\t",index=False,compression=compression,lineterminator="\n");os.replace(part,path)


def main() -> None:
    parser=argparse.ArgumentParser();parser.add_argument("--auc",type=Path,required=True);parser.add_argument("--registry",type=Path,required=True);parser.add_argument("--trajectory",type=Path,required=True);parser.add_argument("--output",type=Path,required=True);args=parser.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    args.output.mkdir(parents=True)
    auc=pd.read_csv(args.auc);registry=pd.read_csv(args.registry,sep="\t");traj=pd.read_csv(args.trajectory,sep="\t",dtype={"cell_id":str})
    cell_col="Cell" if "Cell" in auc.columns else auc.columns[0]
    auc=auc.rename(columns={cell_col:"cell_id"});auc["cell_id"]=auc["cell_id"].astype(str)
    if len(auc)!=50761 or auc["cell_id"].duplicated().any() or set(auc["cell_id"])!=set(traj["cell_id"]):raise RuntimeError("AUCell and trajectory identity mismatch")
    reg_rows=[];candidate_cols={}
    for tf in CANDIDATES:
        choices=[c for c in auc.columns if c==f"{tf}(+)"]
        reg=registry.loc[registry["regulon"].eq(f"{tf}(+)")]
        available=len(choices)==1 and len(reg)==1
        candidate_cols[tf]=choices[0] if available else None
        reg_rows.append({"tf":tf,"expected_regulon":f"{tf}(+)","available":available,"auc_column":choices[0] if choices else "","n_targets":int(reg.iloc[0]["n_targets"]) if len(reg)==1 else 0})
    candidate_registry=pd.DataFrame(reg_rows)
    keep=["cell_id"]+[col for col in candidate_cols.values() if col]
    cells=traj.merge(auc[keep],on="cell_id",how="left",validate="one_to_one")
    for tf,col in candidate_cols.items():
        if col:cells[f"{tf}_AUCell_v1"]=cells[col].to_numpy(float)

    effects=[]
    for (study,patient),frame in cells.groupby(["study_v1","patient"],sort=True):
        base=len(frame)>=50 and frame["state_v2"].nunique()>=2
        for tf,col in candidate_cols.items():
            if not col:
                for endpoint in ENDPOINTS:effects.append({"study":study,"patient":patient,"tf":tf,"endpoint":endpoint,"eligible":False,"n_cells":len(frame),"patient_effect":math.nan})
                continue
            score=f"{tf}_AUCell_v1"
            effects.append({"study":study,"patient":patient,"tf":tf,"endpoint":"Monocle3_continuous","eligible":base,"n_cells":len(frame),"patient_effect":rho(frame["monocle3_pseudotime_v1"],frame[score]) if base else math.nan})
            effects.append({"study":study,"patient":patient,"tf":tf,"endpoint":"DPT_continuous","eligible":base,"n_cells":len(frame),"patient_effect":rho(frame["dpt_pseudotime"],frame[score]) if base else math.nan})
            spp1=frame.loc[frame["state_v2"].eq("SPP1-state Macrophage"),score];mono=frame.loc[frame["state_v2"].eq("Monocyte"),score];c1qc=frame.loc[frame["state_v2"].eq("C1QC-state Macrophage"),score]
            elig_sm=len(spp1)>=10 and len(mono)>=10;elig_sc=len(spp1)>=10 and len(c1qc)>=10
            effects.append({"study":study,"patient":patient,"tf":tf,"endpoint":"SPP1_minus_Monocyte","eligible":elig_sm,"n_cells":len(frame),"patient_effect":float(spp1.median()-mono.median()) if elig_sm else math.nan})
            effects.append({"study":study,"patient":patient,"tf":tf,"endpoint":"SPP1_minus_C1QC","eligible":elig_sc,"n_cells":len(frame),"patient_effect":float(spp1.median()-c1qc.median()) if elig_sc else math.nan})
    effects=pd.DataFrame(effects)
    summaries=[]
    for (endpoint,tf),frame in effects.groupby(["endpoint","tf"],sort=True):
        values=frame.loc[frame["eligible"]&frame["patient_effect"].notna(),"patient_effect"].to_numpy(float);n=len(values);p=signed_rank(values) if n>=3 else math.nan
        summaries.append({"endpoint":endpoint,"tf":tf,"n_eligible_patients":n,"n_positive":int((values>0).sum()),"n_negative":int((values<0).sum()),"n_zero":int((values==0).sum()),"median_patient_effect":float(np.median(values)) if n else math.nan,"wilcoxon_p_two_sided":p})
    summary=pd.DataFrame(summaries);summary["BH_FDR_within_endpoint3TF"]=summary.groupby("endpoint",sort=False)["wilcoxon_p_two_sided"].transform(bh);summary["reporting_tier"]=np.select([summary["n_eligible_patients"]>=10,summary["n_eligible_patients"].between(3,9)],["FDR_AT_LEAST_10","NOMINAL_P_3_TO_9"],default="DIRECTION_ONLY_LT3");summary["statistical_support"]=np.where(summary["n_eligible_patients"]>=10,summary["BH_FDR_within_endpoint3TF"]<.05,np.where(summary["n_eligible_patients"]>=3,summary["wilcoxon_p_two_sided"]<.05,False));summary["direction"]=np.where(summary["median_patient_effect"]>0,"positive",np.where(summary["median_patient_effect"]<0,"negative","zero"))
    integrated=[]
    for tf,frame in summary.groupby("tf",sort=False):
        supported=frame.loc[frame["statistical_support"]];pos=int(supported["direction"].eq("positive").sum());neg=int(supported["direction"].eq("negative").sum());n=len(supported)
        tier={4:"STRONG_ASSOCIATION",3:"SECONDARY_ASSOCIATION",2:"EXPLORATORY_ASSOCIATION",1:"EXPLORATORY_ASSOCIATION",0:"NO_SUPPORTED_ASSOCIATION"}[n]
        direction="positive" if pos and not neg else ("negative" if neg and not pos else ("mixed" if pos and neg else "unsupported"))
        integrated.append({"tf":tf,"n_supported_endpoints_of_4":n,"n_positive_supported":pos,"n_negative_supported":neg,"association_direction":direction,"evidence_tier":tier})
    integrated=pd.DataFrame(integrated)
    state_source=cells.melt(id_vars=["study_v1","patient","state_v2","cell_id"],value_vars=[f"{tf}_AUCell_v1" for tf,col in candidate_cols.items() if col],var_name="tf_score",value_name="AUCell").groupby(["study_v1","patient","state_v2","tf_score"],observed=True).agg(n_cells=("cell_id","size"),median_AUCell=("AUCell","median")).reset_index()
    study_rows=[]
    for (study,endpoint,tf),frame in effects.groupby(["study","endpoint","tf"],sort=True):
        values=frame.loc[frame["eligible"]&frame["patient_effect"].notna(),"patient_effect"].to_numpy(float);n=len(values);study_rows.append({"study":study,"endpoint":endpoint,"tf":tf,"n_eligible_patients":n,"median_patient_effect":float(np.median(values)) if n else math.nan,"positive_fraction":float((values>0).mean()) if n else math.nan,"wilcoxon_p_two_sided":signed_rank(values) if n>=3 else math.nan})
    study_summary=pd.DataFrame(study_rows)
    atomic_tsv(candidate_registry,args.output/"y03_candidate_regulon_registry_v1.tsv");atomic_tsv(cells[["cell_id","study_v1","patient","state_v2","dpt_pseudotime","monocle3_pseudotime_v1",*[f"{tf}_AUCell_v1" for tf,col in candidate_cols.items() if col]]],args.output/"y03_candidate_regulon_cell_scores_v1.tsv.gz",compression="gzip");atomic_tsv(effects,args.output/"y03_candidate_patient_effects_v1.tsv");atomic_tsv(summary,args.output/"y03_candidate_endpoint_summary_v1.tsv");atomic_tsv(integrated,args.output/"y03_candidate_integrated_evidence_v1.tsv");atomic_tsv(state_source,args.output/"y03_candidate_state_heatmap_source_v1.tsv");atomic_tsv(study_summary,args.output/"y03_candidate_study_stratified_summary_v1.tsv")
    status={"status":"PASS_Y03_CANDIDATE_ANALYSIS_COMPLETE_PRESENTATION_PENDING","completed_at":datetime.now().astimezone().isoformat(timespec="seconds"),"n_cells":len(cells),"n_candidates":len(CANDIDATES),"n_available_candidates":int(candidate_registry["available"].sum()),"direction_prespecified":False,"two_sided_testing":True,"figures_created":False,"interpretation":"expression-derived regulon enrichment; not direct TF binding, protein activity, perturbation evidence or causality"};(args.output/"status.json").write_text(json.dumps(status,indent=2)+"\n",encoding="utf-8")


if __name__=="__main__":main()
