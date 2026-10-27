#!/usr/bin/env python3
from __future__ import annotations
import argparse,hashlib,json,os
from datetime import datetime
from pathlib import Path
import numpy as np,pandas as pd

STUDIES={"GSE156632","GSE159115","GSE178481","GSE207493","GSE242299"}
STATES={"Monocyte","Other Macrophage","C1QC-state Macrophage","SPP1-state Macrophage"}
PRIMARY={"monocyte_program":-1,"spp1_state_no_spp1":1,"spp1_vs_c1qc":1}
H5SHA="e4779c46af13a8643e8dec584eeff4322258181f5edab3e55dcd8c20bc28a9e7"
ASHA="f8d0d01dcb76b3a19055accae42be3a2e14c6ba3a9da85db30c9592f2e635e6d"
def sha(p):
 d=hashlib.sha256()
 with Path(p).open("rb") as h:
  for b in iter(lambda:h.read(8*1024*1024),b""):d.update(b)
 return d.hexdigest()
def bh(p):
 p=np.asarray(p,float);out=np.full(len(p),np.nan);m=np.isfinite(p)
 if not m.any():return out
 v=p[m];o=np.argsort(v);q=np.minimum.accumulate((v[o]*len(v)/np.arange(1,len(v)+1))[::-1])[::-1];r=np.empty(len(v));r[o]=np.minimum(q,1);out[m]=r;return out
def main():
 p=argparse.ArgumentParser();p.add_argument("--analysis",type=Path,required=True);p.add_argument("--assignments",type=Path,required=True);p.add_argument("--output",type=Path,required=True);a=p.parse_args()
 if a.output.exists():raise FileExistsError(a.output)
 a.output.mkdir(parents=True);checks=[]
 def ck(name,value,detail):checks.append({"check":name,"status":"PASS" if bool(value) else "FAIL","detail":detail})
 r=a.analysis;status=json.loads((r/"PHASE21_MONOCYTE_SPP1_TRAJECTORY_STATUS_v1.json").read_text());ck("formal_status",status.get("status")=="PASS_ANALYSIS_COMPLETE",str(status.get("status")));ck("input_hash_gate",status.get("input_h5ad_sha256")==H5SHA and status.get("assignments_sha256")==ASHA and sha(a.assignments)==ASHA,"frozen H5AD status plus independently rehashed assignments");parts=list(r.rglob("*.part"));ck("no_part_files",not parts,f"count={len(parts)}")
 manifest=pd.read_csv(r/"OUTPUT_MANIFEST_v1.tsv",sep="\t");ok=True
 for row in manifest.itertuples(index=False):
  f=r/row.file;ok &= f.is_file() and f.stat().st_size==int(row.size_bytes) and sha(f)==row.sha256
 ck("manifest",ok,f"files={len(manifest)}");t=r/"tables";cell=pd.read_csv(t/"myeloid_trajectory_cell_scores_v1.tsv.gz",sep="\t");patient=pd.read_csv(t/"patient_pseudotime_associations_v1.tsv",sep="\t");agg=pd.read_csv(t/"study_primary_inference_v1.tsv",sep="\t");loo=pd.read_csv(t/"leave_one_patient_out_v1.tsv",sep="\t");synth=pd.read_csv(t/"cross_study_direction_synthesis_v1.tsv",sep="\t");state=pd.read_csv(t/"state_trajectory_summary_v1.tsv",sep="\t");graph=pd.read_csv(t/"graph_state_connectivity_v1.tsv",sep="\t");qc=pd.read_csv(t/"signature_gene_coverage_v1.tsv",sep="\t");roots=pd.read_csv(t/"trajectory_root_qc_v1.tsv",sep="\t")
 assignments=pd.read_csv(a.assignments,sep="\t");expected=assignments[assignments.state_v2.isin(STATES)][["exact_cell_id","patient","study_v1","state_v2"]].rename(columns={"exact_cell_id":"cell_id"}).sort_values("cell_id").reset_index(drop=True);got=cell[["cell_id","patient","study_v1","state_v2"]].sort_values("cell_id").reset_index(drop=True);ck("cell_identity",len(cell)==50761 and cell.cell_id.nunique()==50761 and expected.equals(got),f"cells={len(cell)}")
 numeric=cell.drop(columns=["cell_id","patient","study_v1","state_v2"]).to_numpy(float);ck("finite_scores",np.isfinite(numeric).all() and ((cell.dpt_pseudotime>=0)&(cell.dpt_pseudotime<=1)).all(),"finite scores and DPT [0,1]");ck("study_separation",set(cell.study_v1)==STUDIES and cell.groupby("cell_id").study_v1.nunique().max()==1,"five disjoint study-labelled graphs");ck("root_qc",len(roots)==5 and set(roots.study)==STUDIES and roots.root_state.eq("Monocyte").all() and (roots.root_component_fraction>=.95).all(),f"rows={len(roots)}")
 ck("state_summary",len(state)==20 and set(state.state)==STATES and set(state.study)==STUDIES,f"rows={len(state)}");ck("graph_matrix",len(graph)==80 and set(graph.source)==STATES and set(graph.target)==STATES and np.isfinite(graph.source_normalized_fraction).all(),f"rows={len(graph)}");ck("signature_coverage",len(qc)==35 and (qc.present>=qc.minimum_present).all() and set(qc.study)==STUDIES,f"rows={len(qc)}")
 endpoints=set(PRIMARY)|{"spp1_state_full","vascular_remodeling","myeloid_immunosuppression","polyamine_mta_context"};ck("patient_units",set(patient.endpoint)==endpoints and not patient.duplicated(["study","patient","endpoint"]).any() and patient.loc[patient.eligible,"n_cells"].ge(50).all() and patient.loc[patient.eligible,"n_states"].ge(2).all(),f"rows={len(patient)}")
 expected_agg={(s,e) for s in STUDIES for e in PRIMARY};ck("primary_grid",len(agg)==15 and set(zip(agg.study,agg.endpoint))==expected_agg and agg.expected_direction.eq(agg.endpoint.map(PRIMARY)).all(),f"rows={len(agg)}")
 med=patient[patient.eligible&patient.spearman_rho.notna()&patient.endpoint.isin(PRIMARY)].groupby(["study","endpoint"]).spearman_rho.median();joined=agg.set_index(["study","endpoint"]).join(med.rename("recomputed"));ck("aggregate_medians",np.allclose(joined.median_rho,joined.recomputed,equal_nan=True),"patient-level medians")
 fdr_ok=True
 for study_name,f in agg.groupby("study",sort=False):fdr_ok &= np.allclose(f.wilcoxon_fdr_bh_within_study,bh(f.wilcoxon_p),equal_nan=True)
 ck("within_study_fdr",fdr_ok,"three primary endpoints per study");rec=[]
 for row in loo.itertuples(index=False):
  f=patient[(patient.study==row.study)&(patient.endpoint==row.endpoint)&patient.eligible&patient.spearman_rho.notna()&patient.patient.ne(row.left_out_patient)];rec.append(float(f.spearman_rho.median()))
 ck("loo",len(loo)==sum(agg.eligible_patients.fillna(0).astype(int)) and np.allclose(loo.median_rho,rec),f"rows={len(loo)}");ck("cross_study_boundary",len(synth)==3 and set(synth.endpoint)==set(PRIMARY) and synth.interpretation.eq("DESCRIPTIVE_DIRECTION_SYNTHESIS_ONLY").all(),"no pooled cell inference");ck("claim_boundary",status.get("study_separate_graphs") is True and status.get("patient_as_biological_unit") is True and status.get("cross_study_cell_pooling") is False and status.get("lineage_conversion_claim") is False,"state continuity only")
 out=pd.DataFrame(checks);q=a.output/"PHASE21_INDEPENDENT_VALIDATION_STATUS_v1.tsv.part";out.to_csv(q,sep="\t",index=False);os.replace(q,a.output/"PHASE21_INDEPENDENT_VALIDATION_STATUS_v1.tsv");passed=out.status.eq("PASS").all();summary={"schema_version":"phase21_independent_validation_v1","status":"PASS" if passed else "FAIL","validated_at":datetime.now().astimezone().isoformat(timespec="seconds"),"checks_passed":int(out.status.eq("PASS").sum()),"checks_total":len(out),"studies":5,"cells":len(cell),"primary_rows":len(agg)};(a.output/"PHASE21_INDEPENDENT_VALIDATION_SUMMARY_v1.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8");print(json.dumps(summary,indent=2));raise SystemExit(0 if passed else 1)
if __name__=="__main__":main()
