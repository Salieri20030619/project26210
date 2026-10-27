#!/usr/bin/env python3
"""Independent structural and claim-level validation for Phase 15 spatial analysis."""

from __future__ import annotations

import hashlib
import json
import os
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
POLAR = HERE / "spatial_polarization_v1"
EXPR = HERE / "server_expression_results_v1"
PHASE12 = ROOT / "phase12_glmf_spatial_analysis_local_v1" / "server_results_v1"
SUMMARY = HERE / "SPATIAL_HYPOTHESIS_VALIDATION_SUMMARY_v1.md"
OUT = HERE / "independent_validation_v3"

def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""): h.update(b)
    return h.hexdigest()

def atomic_tsv(x,path):
    part=path.with_name(path.name+".part"); x.to_csv(part,sep="\t",index=False,lineterminator="\n"); os.replace(part,path)

def main():
    if OUT.exists(): raise FileExistsError(OUT)
    OUT.mkdir(); checks=[]
    def check(name,ok,detail): checks.append({"check":name,"result":"PASS" if ok else "FAIL","detail":detail})
    required=[
        POLAR/"SPATIAL_POLARIZATION_STATUS_v1.json",
        POLAR/"tables"/"spot_within_lineage_polarization_v1.tsv.gz",
        POLAR/"tables"/"spatial_polarization_correlations_v1.tsv",
        POLAR/"tables"/"spatial_block_stability_summary_v1.tsv",
        EXPR/"SPATIAL_EXPRESSION_VALIDATION_STATUS_v1.json",
        EXPR/"spatial_program_scores_v1.tsv.gz",
        EXPR/"spatial_program_correlations_v1.tsv",
        EXPR/"spatial_program_block_stability_summary_v1.tsv",
        EXPR/"signature_gene_overlap_v1.tsv",
        SUMMARY,
    ]
    check("required_outputs_nonempty",all(p.is_file() and p.stat().st_size>0 for p in required),f"n={len(required)}")
    pstatus=json.loads((POLAR/"SPATIAL_POLARIZATION_STATUS_v1.json").read_text()); estatus=json.loads((EXPR/"SPATIAL_EXPRESSION_VALIDATION_STATUS_v1.json").read_text())
    check("status_and_no_spot_inference",pstatus["status"]=="PASS_DESCRIPTIVE_SPATIAL_VALIDATION" and estatus["status"]=="PASS_DESCRIPTIVE_SPATIAL_EXPRESSION_VALIDATION" and pstatus["spots"]==estatus["spots"]==7079 and pstatus["biological_units"]==estatus["biological_units"]==2 and not pstatus["spot_as_replicate_inference"] and not estatus["spot_as_replicate_inference"],"spots=7079 units=2 spot_inference=false")
    pv=pd.read_csv(POLAR/"tables"/"spot_within_lineage_polarization_v1.tsv.gz",sep="\t"); ev=pd.read_csv(EXPR/"spatial_program_scores_v1.tsv.gz",sep="\t")
    id_ok=len(pv)==len(ev)==7079 and not pv[["sample","spot_id"]].duplicated().any() and not ev[["sample","spot_id"]].duplicated().any() and set(pv["sample"])==set(ev["sample"])=={"GLMF1","GLMF2"}
    check("spot_identity_conservation",id_ok,"GLMF1=3199 GLMF2=3880")
    overlap=pd.read_csv(EXPR/"signature_gene_overlap_v1.tsv",sep="\t"); minima={"f1_top100_weighted":70,"c1qc_state":2,"spp1_state_no_spp1":3,"spp1_state_full":4,"hypoxia_angiogenesis":4,"tcell_dysfunction":5,"cytotoxicity":3}
    overlap_ok=all((overlap.loc[overlap.signature_id==s,"present_genes"]>=m).all() for s,m in minima.items()) and set(overlap["sample"])=={"GLMF1","GLMF2"}
    check("frozen_signature_gene_coverage",overlap_ok,"F1=74/100; C1QC=3/3; hypoxia=5/6; Tdysfunction=7/8")
    corr=pd.read_csv(EXPR/"spatial_program_correlations_v1.tsv",sep="\t"); pcorr=pd.read_csv(POLAR/"tables"/"spatial_polarization_correlations_v1.tsv",sep="\t")
    check("no_inferential_p_value_columns",not any("p_value" in c.lower() for c in list(corr.columns)+list(pcorr.columns)),"descriptive rho only")
    controls=corr[(corr.relation.isin(["C1QC_mapping_positive_control","SPP1_mapping_positive_control"]))&(corr.outcome=="lineage_depth_residual")]
    check("mapping_positive_controls",len(controls)==48 and (controls.spearman_rho_descriptive>0).all(),f"rows={len(controls)} rho_min={controls.spearman_rho_descriptive.min():.3f}")
    f1c=corr[(corr.relation=="F1polar_to_C1QC_expression")&(corr.outcome=="lineage_depth_residual")]
    hetero=len(f1c)==24 and (f1c[f1c["sample"]=="GLMF1"].spearman_rho_descriptive>0).all() and (f1c[f1c["sample"]=="GLMF2"].spearman_rho_descriptive<0).all()
    check("f1_c1qc_section_heterogeneity",hetero,"GLMF1 all positive; GLMF2 all negative")
    hyp=corr[(corr.relation=="F1polar_to_hypoxia_expression")&(corr.outcome=="lineage_depth_residual")]
    td=corr[(corr.relation=="SPP1polar_to_Tdysfunction_expression")&(corr.outcome=="lineage_depth_residual")]
    check("hypoxia_not_supported_both_sections",len(hyp)==24 and (hyp.spearman_rho_descriptive<0).all(),f"rows={len(hyp)}")
    check("t_dysfunction_coarse_test_not_supported",len(td)==24 and (td.spearman_rho_descriptive<0).all(),f"rows={len(td)}")
    blocks=pd.read_csv(EXPR/"spatial_program_block_stability_summary_v1.tsv",sep="\t"); key=blocks[(blocks.outcome=="lineage_depth_residual")&blocks.relation.isin(["F1polar_to_C1QC_expression","F1polar_to_hypoxia_expression","SPP1polar_to_Tdysfunction_expression"])]
    block_ok=len(key)==24 and ((key.positive_fraction==1)|(key.negative_fraction==1)).all()
    check("leave_one_block_direction_stability",block_ok,f"rows={len(key)}")
    text=" ".join(SUMMARY.read_text(encoding="utf-8").split()); terms=["two sections","not independently reproduced","not supported by the two spatial sections","failed coarse spatial validation"]
    check("summary_claim_boundaries",all(t in text for t in terms),";".join(terms))
    passed=all(x["result"]=="PASS" for x in checks); now=datetime.now().astimezone().isoformat(timespec="seconds"); checks.append({"check":"OVERALL","result":"PASS" if passed else "FAIL","detail":now}); atomic_tsv(pd.DataFrame(checks),OUT/"SPATIAL_HYPOTHESIS_VALIDATION_STATUS_v1.tsv")
    manifest=pd.DataFrame([{"relative_path":str(p.relative_to(ROOT)).replace("\\","/"),"bytes":p.stat().st_size,"sha256":sha256(p)} for p in required]); atomic_tsv(manifest,OUT/"SPATIAL_HYPOTHESIS_VALIDATION_FILE_MANIFEST_v1.tsv")
    status={"schema_version":"phase15_spatial_hypothesis_independent_validation_v1","status":"PASS" if passed else "FAIL","completed_at":now,"checks_passed":sum(x["result"]=="PASS" for x in checks),"checks_total":len(checks),"biological_units":2,"spot_as_replicate_inference":False}; part=OUT/"SPATIAL_HYPOTHESIS_VALIDATION_SUMMARY_v1.json.part"; part.write_text(json.dumps(status,indent=2)+"\n"); os.replace(part,OUT/"SPATIAL_HYPOTHESIS_VALIDATION_SUMMARY_v1.json"); print(json.dumps(status,indent=2));
    if not passed: raise SystemExit(1)

if __name__=="__main__": main()
