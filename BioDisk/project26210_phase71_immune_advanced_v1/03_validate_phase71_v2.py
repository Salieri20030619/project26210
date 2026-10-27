#!/usr/bin/env python3
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT=Path("/mnt/BioDisk/wangmingyuan/project26210_phase71_immune_advanced_v1")
RES=ROOT/"results_v2"
CHECKS=[]


def add(cid,ok,detail): CHECKS.append({"check_id":cid,"status":"PASS" if ok else "FAIL","detail":str(detail)})
def read(name): return pd.read_csv(RES/name,sep="\t")
def bh(p):
    p=np.asarray(p,float); n=len(p); order=np.argsort(p); ranked=p[order]; q=np.minimum.accumulate((ranked*n/np.arange(1,n+1))[::-1])[::-1]; q=np.minimum(q,1); out=np.empty(n); out[order]=q; return out
def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""): h.update(b)
    return h.hexdigest()

summary=read("phase71_relation_summary_all_v1.tsv")
cross=read("phase71_cross_cohort_decisions_v1.tsv")
bridge=read("phase71_primary_bridge_decision_v1.tsv")
coverage=read("phase71_program_coverage_all_v1.tsv")
status=read("phase71_dataset_status_all_v1.tsv")
add("dataset_status",len(status)==3 and status.status.eq("PASS_PHASE71_DATASET_V1").all(),status[["dataset","status"]].to_dict("records"))
add("summary_dimensions",len(summary)==39 and len(summary.query("tier=='PRIMARY'"))==18,len(summary))
add("cross_dimensions",len(cross)==13 and cross.relation_id.nunique()==13,len(cross))
add("coverage",len(coverage)==24 and coverage.eligible.all(),len(coverage))
q_ok=True
decision_ok=True
for (dataset,tier),g in summary.groupby(["dataset","tier"],sort=False):
    q_ok &= np.allclose(g.q_value.to_numpy(float),bh(g.p_value.to_numpy(float)),rtol=1e-10,atol=1e-12)
for r in summary.itertuples(index=False):
    if r.expected_direction=="two_sided": expected="ASSOCIATED_TWO_SIDED" if r.q_value<0.05 else "NOT_SUPPORTED"
    else:
        supported=r.q_value<0.05 and ((r.expected_direction=="positive" and r.median_adjusted_rho>0 and r.positive_fraction>=2/3) or (r.expected_direction=="negative" and r.median_adjusted_rho<0 and r.positive_fraction<=1/3))
        expected="COHORT_SUPPORTED" if supported else ("DIRECTION_CONFLICT" if r.q_value<0.05 else "NOT_SUPPORTED")
    decision_ok &= expected==r.decision
add("bh_recalculation",q_ok,q_ok)
add("dataset_decisions_recalculated",decision_ok,decision_ok)
cross_ok=True
for r in cross.itertuples(index=False):
    if r.expected_direction=="two_sided":
        n=sum(getattr(r,f"{d}_decision")=="ASSOCIATED_TWO_SIDED" for d in ("GSE264586","GSE294109")); expected="TWO_EXTERNAL_ASSOCIATED" if n==2 else ("ONE_EXTERNAL_ASSOCIATED" if n==1 else "NOT_SUPPORTED")
    else:
        n=sum(getattr(r,f"{d}_decision")=="COHORT_SUPPORTED" for d in ("GSE264586","GSE294109")); expected="TWO_EXTERNAL_SUPPORTED" if n==2 else ("ONE_EXTERNAL_SUPPORTED" if n==1 else "NOT_SUPPORTED")
    cross_ok &= expected==r.external_decision
add("cross_decisions_recalculated",cross_ok,cross_ok)
required=cross.loc[cross.relation_id.isin(["P01","P02","P03"])]
expected_bridge="TRANSCRIPTIONAL_BRIDGE_SUPPORTED" if required.external_decision.eq("TWO_EXTERNAL_SUPPORTED").all() else "TRANSCRIPTIONAL_BRIDGE_NOT_FULLY_SUPPORTED"
add("bridge_decision",len(bridge)==1 and bridge.iloc[0].decision==expected_bridge,expected_bridge)
manifest=read("FILE_MANIFEST_SHA256_v1.tsv")
hash_ok=all((RES/r.file).is_file() and sha(RES/r.file)==r.sha256 for r in manifest.itertuples(index=False))
add("output_hashes",hash_ok,len(manifest))
out=ROOT/"independent_validation_v2"; out.mkdir(exist_ok=True)
pd.DataFrame(CHECKS).to_csv(out/"validation_checks_v2.tsv",sep="\t",index=False,lineterminator="\n")
payload={"checks":len(CHECKS),"passed":sum(x["status"]=="PASS" for x in CHECKS),"failed":sum(x["status"]=="FAIL" for x in CHECKS)}
(out/"validation_summary_v2.json").write_text(json.dumps(payload,indent=2)+"\n")
if payload["failed"]==0: (out/"PHASE71_INDEPENDENT_VALIDATION_PASS_v2").write_text("PASS\n")
print(json.dumps(payload))
raise SystemExit(0 if payload["failed"]==0 else 1)
