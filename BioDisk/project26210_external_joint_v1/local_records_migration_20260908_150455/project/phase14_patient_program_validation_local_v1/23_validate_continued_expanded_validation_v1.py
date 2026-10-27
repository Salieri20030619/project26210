#!/usr/bin/env python3
"""Independent reconciliation of continued expanded-malignancy validation."""

from __future__ import annotations

import hashlib
import json
import os
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
TIER = HERE / "posthoc_malignancy_evidence_tier_programs_v5"
THRESH = HERE / "posthoc_tier_threshold_validation_v1"
METHOD = HERE / "posthoc_method_specific_supportive_validation_v1"
NEG = HERE / "posthoc_negative_control_layers_v1"
STRICT = HERE / "patient_analysis_final_cnv_v6" / "within_dataset_patient_associations_v6.tsv"
SUMMARY = HERE / "CONTINUED_EXPANDED_MALIGNANCY_VALIDATION_SUMMARY_v1.md"
OUT = HERE / "posthoc_continued_expanded_validation_check_v1"

def sha256(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""): h.update(b)
    return h.hexdigest()

def atomic_tsv(x,path):
    part=path.with_name(path.name+".part"); x.to_csv(part,sep="\t",index=False,lineterminator="\n"); os.replace(part,path)

def main():
    if OUT.exists(): raise FileExistsError(OUT)
    OUT.mkdir(); rows=[]
    def check(name,ok,detail): rows.append({"check":name,"result":"PASS" if ok else "FAIL","detail":detail})
    required=[TIER/"MALIGNANCY_EVIDENCE_TIER_PROGRAM_STATUS_v1.tsv",TIER/"cell_evidence_method_and_annotation_audit_v1.tsv",TIER/"cross_sectional_patient_evidence_tier_programs_v1.tsv",THRESH/"TIER_THRESHOLD_VALIDATION_STATUS_v1.json",THRESH/"tier_threshold_malignant_relations_v1.tsv",METHOD/"METHOD_SPECIFIC_SUPPORTIVE_VALIDATION_STATUS_v1.json",METHOD/"method_specific_supportive_key_relations_v1.tsv",NEG/"NEGATIVE_CONTROL_LAYER_STATUS_v1.json",NEG/"negative_control_layer_key_relations_v1.tsv",SUMMARY]
    check("required_outputs_nonempty",all(p.is_file() and p.stat().st_size>0 for p in required),f"n={len(required)}")
    ts=pd.read_csv(TIER/"MALIGNANCY_EVIDENCE_TIER_PROGRAM_STATUS_v1.tsv",sep="\t").iloc[0]
    check("tier_program_status",ts.status=="PASS_POSTHOC_TIER_PROGRAMS" and ts.tumor_like_cells==40302 and ts.high_confidence_cells==16645 and ts.supportive_cells==17775 and bool(ts.rebuild_priority),f"tumor={ts.tumor_like_cells} high={ts.high_confidence_cells} support={ts.supportive_cells}")
    audit=pd.read_csv(TIER/"cell_evidence_method_and_annotation_audit_v1.tsv",sep="\t"); sup=audit[audit.evidence_tier=="supportive"]
    copy=int(sup.loc[sup.evidence_detail.str.contains("copykat_only"),"n_cells"].sum()); inf=int(sup.loc[sup.evidence_detail.str.contains("infercnv_only"),"n_cells"].sum())
    check("supportive_method_composition",copy==15119 and inf==2656 and copy+inf==17775,f"copykat_only={copy} infercnv_only={inf}")
    th=json.loads((THRESH/"TIER_THRESHOLD_VALIDATION_STATUS_v1.json").read_text(encoding="utf-8")); x=pd.read_csv(THRESH/"tier_threshold_malignant_relations_v1.tsv",sep="\t")
    check("threshold_grid",th["status"]=="PASS_POSTHOC_TIER_THRESHOLD_VALIDATION" and th["thresholds"]==[10,20,30,50] and len(x)==1080 and th["expanded_threshold20_tested"]==54 and th["expanded_threshold20_uniform9_fdr"]==0,"thresholds=10,20,30,50 rows=1080 expanded20_tested=54 fdr=0")
    key=x[(x.dataset=="GSE207493")&(x.relation=="H2_F1_to_C1QC")&(x.evidence_tier.isin(["expanded","supportive"]))]
    key_ok=len(key)==16 and key.spearman_rho.notna().all() and (key.spearman_rho>0).all() and (key.loo_sign_stability_fraction==1).all()
    check("gse207493_f1_c1qc_tier_threshold_stability",key_ok,f"rows={len(key)} rho_range={key.spearman_rho.min():.3f}-{key.spearman_rho.max():.3f}")
    ms=json.loads((METHOD/"METHOD_SPECIFIC_SUPPORTIVE_VALIDATION_STATUS_v1.json").read_text()); m=pd.read_csv(METHOD/"method_specific_supportive_key_relations_v1.tsv",sep="\t"); mk=m[(m.dataset=="GSE207493")&(m.relation=="H2_F1_to_C1QC")&m.spearman_rho.notna()]
    check("method_specific_direction",ms["status"]=="PASS_POSTHOC_METHOD_SPECIFIC_SUPPORTIVE_VALIDATION" and len(mk)==4 and (mk.spearman_rho>0).all(),f"testable={len(mk)} positive={(mk.spearman_rho>0).sum()}")
    ns=json.loads((NEG/"NEGATIVE_CONTROL_LAYER_STATUS_v1.json").read_text()); n=pd.read_csv(NEG/"negative_control_layer_key_relations_v1.tsv",sep="\t"); nk=n[(n.dataset=="GSE207493")&(n.relation=="H2_F1_to_C1QC")]
    check("negative_control_layers",ns["status"]=="PASS_POSTHOC_NEGATIVE_CONTROL_LAYERS" and len(nk)==6 and nk.spearman_rho.notna().all() and (nk.spearman_rho>0).all() and (nk.loo_sign_stability_fraction==1).all(),f"rows={len(nk)} positive={(nk.spearman_rho>0).sum()}")
    strict=pd.read_csv(STRICT,sep="\t"); check("formal_strict_v6_unchanged",int((strict.p_fdr_bh_within_dataset_method<.05).sum())==0,"strict_v6_fdr=0")
    text=" ".join(SUMMARY.read_text(encoding="utf-8").split()); terms=["rebuild-first","15,119","10, 20, 30, and 50","formal strict v6 result remains unchanged"]
    check("summary_scope_and_claims",all(t in text for t in terms),";".join(terms))
    passed=all(r["result"]=="PASS" for r in rows); now=datetime.now().astimezone().isoformat(timespec="seconds"); rows.append({"check":"OVERALL","result":"PASS" if passed else "FAIL","detail":now}); atomic_tsv(pd.DataFrame(rows),OUT/"CONTINUED_EXPANDED_VALIDATION_STATUS_v1.tsv")
    manifest=pd.DataFrame([{"relative_path":str(p.relative_to(HERE.parent)).replace("\\","/"),"bytes":p.stat().st_size,"sha256":sha256(p)} for p in required]); atomic_tsv(manifest,OUT/"CONTINUED_EXPANDED_VALIDATION_FILE_MANIFEST_v1.tsv")
    status={"schema_version":"phase14_continued_expanded_validation_check_v1","status":"PASS" if passed else "FAIL","completed_at":now,"checks_passed":sum(r["result"]=="PASS" for r in rows),"checks_total":len(rows),"formal_strict_v6_replaced":False}; part=OUT/"CONTINUED_EXPANDED_VALIDATION_SUMMARY_v1.json.part"; part.write_text(json.dumps(status,indent=2)+"\n"); os.replace(part,OUT/"CONTINUED_EXPANDED_VALIDATION_SUMMARY_v1.json"); print(json.dumps(status,indent=2));
    if not passed: raise SystemExit(1)

if __name__=="__main__": main()
