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
from scipy.stats import wilcoxon


ENDPOINTS = ["malignant_to_macrophage", "macrophage_to_malignant", "malignant_to_endothelial", "endothelial_to_malignant"]


def bh(values: pd.Series) -> np.ndarray:
    p = values.to_numpy(float); order = np.argsort(p); ranked = p[order] * len(p) / np.arange(1, len(p)+1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]; out = np.empty(len(p)); out[order] = np.minimum(ranked, 1); return out


def signed_rank(values: np.ndarray) -> float:
    values = values[np.isfinite(values)]
    if len(values) < 3 or np.all(values == 0): return 1.0
    return float(wilcoxon(values, alternative="two-sided", zero_method="wilcox", correction=False, method="auto").pvalue)


def atomic_tsv(frame: pd.DataFrame, path: Path, gzip_output: bool = False) -> None:
    part = Path(str(path) + ".part")
    frame.to_csv(part, sep="\t", index=False, compression="gzip" if gzip_output else None, lineterminator="\n")
    os.replace(part, path)


def main() -> None:
    parser = argparse.ArgumentParser(); parser.add_argument("--batches", nargs="+", type=Path, required=True); parser.add_argument("--output", type=Path, required=True); args = parser.parse_args()
    if args.output.exists(): raise FileExistsError(args.output)
    args.output.mkdir(parents=True)
    interactions=[];effects=[];lr=[];qc=[];statuses=[]
    for batch in args.batches:
        status=json.loads((batch/"status.json").read_text());statuses.append(status)
        if status["status"]!="PASS_C01_CELLCHAT_ANALYSIS_COMPLETE_PRESENTATION_PENDING": raise RuntimeError(f"Batch failed: {batch}")
        interactions.append(pd.read_csv(batch/"c01_cellchat_patient_interactions_v1.tsv.gz",sep="\t"))
        effects.append(pd.read_csv(batch/"c01_cellchat_patient_endpoint_effects_v1.tsv",sep="\t"))
        lr.append(pd.read_csv(batch/"c01_cellchat_lr_pair_differences_v1.tsv.gz",sep="\t"))
        qc.append(pd.read_csv(batch/"c01_cellchat_patient_run_qc_v1.tsv",sep="\t"))
    interactions=pd.concat(interactions,ignore_index=True);effects=pd.concat(effects,ignore_index=True);lr=pd.concat(lr,ignore_index=True);qc=pd.concat(qc,ignore_index=True)
    patients=sorted(effects.patient.unique())
    if len(patients)!=14 or effects.duplicated(["patient","endpoint"]).any() or set(effects.endpoint)!=set(ENDPOINTS): raise RuntimeError("Final CellChat endpoint grid mismatch")
    summary=[]
    for endpoint,frame in effects.groupby("endpoint",sort=True):
        values=frame.loc[frame.eligible.astype(bool),"high_minus_low"].to_numpy(float)
        summary.append({"endpoint":endpoint,"n_eligible_patients":len(values),"n_positive":int((values>0).sum()),"n_negative":int((values<0).sum()),"n_zero":int((values==0).sum()),"median_high_minus_low":float(np.median(values)),"wilcoxon_p_two_sided":signed_rank(values)})
    summary=pd.DataFrame(summary);summary["BH_FDR_4_endpoints"]=bh(summary.wilcoxon_p_two_sided);summary["statistical_support"]=summary.BH_FDR_4_endpoints.lt(0.05);summary["direction"]=np.where(summary.median_high_minus_low>0,"positive",np.where(summary.median_high_minus_low<0,"negative","zero"))
    pair_rows=[]
    for (family,pair,pathway),frame in lr.groupby(["direction_family","pair_id","pathway_name"],sort=True):
        by_patient=frame.groupby("patient").high_minus_low.sum().reindex(patients,fill_value=0).to_numpy(float)
        pair_rows.append({"direction_family":family,"pair_id":pair,"pathway_name":pathway,"n_patients":len(patients),"n_positive":int((by_patient>0).sum()),"n_negative":int((by_patient<0).sum()),"n_zero":int((by_patient==0).sum()),"median_high_minus_low":float(np.median(by_patient)),"wilcoxon_p_two_sided":signed_rank(by_patient)})
    pair_summary=pd.DataFrame(pair_rows)
    pair_summary["BH_FDR_within_direction_family"]=pair_summary.groupby("direction_family",sort=False).wilcoxon_p_two_sided.transform(bh)
    atomic_tsv(interactions,args.output/"c01_cellchat_patient_interactions_v1.tsv.gz",True);atomic_tsv(effects,args.output/"c01_cellchat_patient_endpoint_effects_v1.tsv");atomic_tsv(summary,args.output/"c01_cellchat_endpoint_summary_v1.tsv");atomic_tsv(lr,args.output/"c01_cellchat_lr_pair_differences_v1.tsv.gz",True);atomic_tsv(pair_summary,args.output/"c01_cellchat_lr_pair_summary_v1.tsv");atomic_tsv(qc,args.output/"c01_cellchat_patient_run_qc_v1.tsv")
    payload={"status":"PASS_C01_CELLCHAT_ANALYSIS_COMPLETE_PRESENTATION_PENDING","completed_at":datetime.now().astimezone().isoformat(timespec="seconds"),"n_batches":len(args.batches),"n_patients":len(patients),"n_interaction_rows":len(interactions),"two_sided_testing":True,"direction_prespecified":False,"figures_created":False,"interpretation":"patient-stratified CellChat RNA-based communication probability; not measured protein binding or causality"}
    part=args.output/"status.json.part";part.write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8");os.replace(part,args.output/"status.json")


if __name__=="__main__":main()
