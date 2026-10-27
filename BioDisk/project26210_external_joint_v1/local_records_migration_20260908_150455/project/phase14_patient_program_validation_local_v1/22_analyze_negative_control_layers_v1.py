#!/usr/bin/env python3
"""Analyze negative-control expanded layers for key malignant relations."""

from __future__ import annotations

import json
import os
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
PROGRAMS = HERE / "posthoc_malignancy_evidence_tier_programs_v5" / "cross_sectional_patient_evidence_tier_programs_v1.tsv"
EXPANDED_RUN = HERE / "runs" / "20260812_expanded_malignancy_v2"
OUT = HERE / "posthoc_negative_control_layers_v1"
DATASETS = ["GSE178481", "GSE202374", "GSE207493", "GSE264586", "GSE294109"]
GROUPS = ["expanded", "expanded_exclude_low_annotation", "high_plus_infercnv_only_supportive"]
METHODS = [("primary_gene_z_mean", "mean_primary"), ("sensitivity_mean_log_normalized", "mean_sensitivity")]
RELATIONS = [
    ("H2_F1_to_C1QC", "f1_top100_weighted", "c1qc_state", "myeloid"),
    ("H2_F1_to_APOE_C1QC", "f1_top100_weighted", "apoe_tam_expression_proxy", "myeloid"),
    ("H3_F1_to_hypoxia", "f1_top100_weighted", "hypoxia_angiogenesis", "malignant"),
    ("H1_TYMP_to_SPP1", "tymp_expression", "spp1_state_no_spp1", "myeloid"),
]

def atomic_tsv(x, path):
    part = path.with_name(path.name + ".part"); x.to_csv(part, sep="\t", index=False, na_rep="NA", lineterminator="\n"); os.replace(part, path)

def get_nonmal(dataset):
    x = pd.read_csv(EXPANDED_RUN / dataset / "tables" / "sample_program_summary_v1.tsv", sep="\t")
    x = x[x.analysis_compartment != "malignant_provisional"].copy()
    if dataset == "GSE178481": x = x[~x.sample_id.astype(str).str.contains("Normal", case=False, regex=False)]
    if dataset == "GSE264586": x = x[x.timepoint == "Baseline"]
    rows=[]; keys=["dataset","patient_key","patient_id","signature_id","analysis_compartment"]
    for k,g in x.groupby(keys,dropna=False):
        w=g.n_cells.to_numpy(float); rows.append(dict(zip(keys,k),n_cells=int(w.sum()),mean_primary=np.average(g.mean_primary,weights=w),mean_sensitivity=np.average(g.mean_sensitivity,weights=w)))
    return pd.DataFrame(rows)

def main():
    if OUT.exists(): raise FileExistsError(OUT)
    OUT.mkdir(); prog=pd.read_csv(PROGRAMS,sep="\t"); rows=[]
    for dataset in DATASETS:
        non=get_nonmal(dataset)
        for method,value in METHODS:
            n=non[non.n_cells>=20].copy(); n["key"]=n.signature_id+"__"+n.analysis_compartment
            nw=n.pivot(index=["dataset","patient_key","patient_id"],columns="key",values=value).reset_index()
            for group in GROUPS:
                m=prog[(prog.dataset==dataset)&(prog.evidence_tier==group)&(prog.n_cells>=20)]
                mw=m.pivot(index=["dataset","patient_key","patient_id"],columns="signature_id",values=value).reset_index()
                wide=mw.merge(nw,on=["dataset","patient_key","patient_id"],how="outer")
                for relation,xsig,ysig,ycomp in RELATIONS:
                    ycol=ysig if ycomp=="malignant" else f"{ysig}__{ycomp}"; cols=[c for c in [xsig,ycol] if c in wide]
                    sub=wide[["patient_key",*cols]].dropna() if len(cols)==2 else pd.DataFrame(); count=len(sub)
                    if count>=5 and sub[xsig].nunique()>1 and sub[ycol].nunique()>1:
                        t=stats.spearmanr(sub[xsig],sub[ycol]); rho=float(t.statistic); p=float(t.pvalue)
                        xx=sub[xsig].to_numpy(float); yy=sub[ycol].to_numpy(float); loo=[float(stats.spearmanr(np.delete(xx,i),np.delete(yy,i)).statistic) for i in range(count)]
                        stable=float(np.mean(np.sign(loo)==np.sign(rho))) if rho else np.nan
                    else: rho=p=stable=np.nan
                    rows.append({"dataset":dataset,"method":method,"control_layer":group,"relation":relation,"n_patients":count,"spearman_rho":rho,"p_value":p,"loo_sign_stability_fraction":stable})
    out=pd.DataFrame(rows); atomic_tsv(out,OUT/"negative_control_layer_key_relations_v1.tsv")
    key=out[(out.dataset=="GSE207493")&(out.relation=="H2_F1_to_C1QC")]
    status={"schema_version":"phase14_negative_control_layers_v1","status":"PASS_POSTHOC_NEGATIVE_CONTROL_LAYERS","completed_at":datetime.now().astimezone().isoformat(timespec="seconds"),"rows":len(out),"gse207493_f1_c1qc_testable":int(key.spearman_rho.notna().sum()),"gse207493_f1_c1qc_positive":int((key.spearman_rho>0).sum()),"formal_strict_v6_replaced":False}
    part=OUT/"NEGATIVE_CONTROL_LAYER_STATUS_v1.json.part"; part.write_text(json.dumps(status,indent=2)+"\n",encoding="utf-8"); os.replace(part,OUT/"NEGATIVE_CONTROL_LAYER_STATUS_v1.json"); print(json.dumps(status,indent=2))

if __name__=="__main__": main()
