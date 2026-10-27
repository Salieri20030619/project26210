#!/usr/bin/env python3
import argparse,hashlib,json,os
from datetime import datetime
from pathlib import Path
import numpy as np,pandas as pd

PRIMARY=["neighbor_TYMP_to_SPP1_vs_C1QC_expression","neighbor_TYMP_to_SPP1_component_polarization","neighbor_TYMP_to_SPP1_state_no_SPP1"]
def sha256(p):
 d=hashlib.sha256()
 with Path(p).open("rb") as h:
  for b in iter(lambda:h.read(8*1024*1024),b""): d.update(b)
 return d.hexdigest()
def write(frame,path):
 part=Path(str(path)+".part"); frame.to_csv(part,sep="\t",index=False,lineterminator="\n"); os.replace(part,path)
def main():
 p=argparse.ArgumentParser();p.add_argument("--analysis",type=Path,required=True);p.add_argument("--phase18-validation",type=Path,required=True);p.add_argument("--output",type=Path,required=True);a=p.parse_args()
 if a.output.exists(): raise FileExistsError(a.output)
 a.output.mkdir(parents=True); checks=[]
 def ck(n,x,d): checks.append({"check":n,"status":"PASS" if bool(x) else "FAIL","detail":d})
 root=a.analysis; status=json.loads((root/"PHASE19_TYMP_SPP1_SPATIAL_STATUS_v1.json").read_text()); ck("formal_status",status.get("status")=="PASS_ANALYSIS_COMPLETE",str(status.get("status"))); ck("phase18_validation_gate",json.loads((a.phase18_validation/"PHASE18_INDEPENDENT_VALIDATION_SUMMARY_v1.json").read_text()).get("status")=="PASS","authoritative v2")
 parts=list(root.rglob("*.part"));ck("no_part_files",not parts,f"count={len(parts)}")
 manifest=pd.read_csv(root/"OUTPUT_MANIFEST_v1.tsv",sep="\t"); ok=True
 for r in manifest.itertuples(index=False):
  f=root/r.file;ok &= f.is_file() and f.stat().st_size==int(r.size_bytes) and sha256(f)==r.sha256
 ck("manifest_size_sha256",ok,f"files={len(manifest)}")
 t=root/"tables";spot=pd.read_csv(t/"spot_scores_v1.tsv.gz",sep="\t");section=pd.read_csv(t/"section_statistics_v1.tsv",sep="\t");primary=pd.read_csv(t/"primary_tjrcc_section_statistics_v1.tsv",sep="\t");agg=pd.read_csv(t/"primary_across_section_inference_v1.tsv",sep="\t");loo=pd.read_csv(t/"primary_leave_one_section_out_v1.tsv",sep="\t");niche=pd.read_csv(t/"niche_tymp_spp1_summary_v1.tsv",sep="\t");qc=pd.read_csv(t/"signature_gene_coverage_v1.tsv",sep="\t")
 ck("section_spot_conservation",len(spot)==41688 and spot["sample"].nunique()==14 and not spot.duplicated(["sample","spot_id"]).any(),f"spots={len(spot)} sections={spot['sample'].nunique()}")
 vals=spot[["tymp_exposure","neighbor_tymp_exposure","spp1_vs_c1qc_expression","spp1_state_no_spp1_residual","spp1_component_polarization"]].to_numpy(float);ck("spot_scores_finite",np.isfinite(vals).all() and ((spot.spp1_component_polarization>=0)&(spot.spp1_component_polarization<=1)).all(),"finite and polarization in [0,1]")
 ck("niche_identity_conservation",len(niche)==42 and niche["sample"].nunique()==14 and set(niche.niche)==set(spot.niche),f"niches={len(niche)}")
 ck("signature_coverage",len(qc)==14*6 and (qc.present_genes>=1).all(),f"rows={len(qc)}")
 expected=section[(section.cohort=="TJRCC")&(section.role=="primary")&(section.gate=="q25")].sort_values(["relation","sample"]).reset_index(drop=True); got=primary.sort_values(["relation","sample"]).reset_index(drop=True);ck("primary_endpoint_conservation",len(got)==36 and got["sample"].nunique()==12 and got.relation.nunique()==3 and expected.equals(got),f"rows={len(got)}")
 ck("no_spot_p_values",not any("p_value" in c for c in spot.columns) and not any("p_value" in c for c in section.columns),"section rhos only")
 med=primary.groupby("relation").spearman_rho.median();joined=agg.set_index("relation");ck("aggregate_medians_recomputed",set(agg.relation)==set(PRIMARY) and np.allclose(joined.loc[med.index,"median_rho"],med),"3 relations")
 rec=[]
 for relation in PRIMARY:
  sub=primary[primary.relation==relation]
  for sample in sorted(sub["sample"]): rec.append((relation,sample,float(np.median(sub[sub["sample"]!=sample].spearman_rho))))
 ref=pd.DataFrame(rec,columns=["relation","left_out_sample","recomputed"]);m=loo.merge(ref,on=["relation","left_out_sample"]);ck("leave_one_section_out_recomputed",len(m)==36 and np.allclose(m.median_rho,m.recomputed),f"rows={len(m)}")
 allowed={"PRIMARY_SPATIAL_SUPPORT","DIRECTIONAL_NOT_FDR","HETEROGENEOUS","NOT_SUPPORTED"};ck("decision_rule_domain",set(agg.decision)<=allowed and ((agg.decision=="PRIMARY_SPATIAL_SUPPORT")<=((agg.median_rho>=.1)&(agg.positive_sections>=8)&(agg.wilcoxon_fdr_bh<.05)&(agg.bootstrap_ci_low>0)&(agg.loo_positive_fraction==1))).all(),"frozen rule")
 ck("cohort_separation",set(spot.loc[spot.cohort=="GLMF","sample"])=={"GLMF1","GLMF2"} and spot.loc[spot.cohort=="TJRCC","sample"].nunique()==12,"GLMF descriptive; TJRCC inference")
 frame=pd.DataFrame(checks);write(frame,a.output/"PHASE19_INDEPENDENT_VALIDATION_STATUS_v1.tsv");passed=frame.status.eq("PASS").all();summary={"schema_version":"phase19_independent_validation_v1","status":"PASS" if passed else "FAIL","validated_at":datetime.now().astimezone().isoformat(timespec="seconds"),"checks_passed":int(frame.status.eq("PASS").sum()),"checks_total":len(frame),"sections":14,"spots":len(spot),"primary_rows":len(primary)};(a.output/"PHASE19_INDEPENDENT_VALIDATION_SUMMARY_v1.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8");print(json.dumps(summary,indent=2));raise SystemExit(0 if passed else 1)
if __name__=="__main__":main()
