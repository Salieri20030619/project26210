#!/usr/bin/env python3
import argparse,hashlib,json,os
from datetime import datetime
from pathlib import Path
import numpy as np,pandas as pd
PRIMARY=["SPP1_polarization_to_vascular_compartment","SPP1_vascular_niche_to_T_cell_balance","SPP1_vascular_niche_to_cytotoxicity","SPP1_vascular_niche_to_hypoxia_angiogenesis","SPP1_vascular_niche_to_invasion_EMT"]
def sha(p):
 d=hashlib.sha256()
 with Path(p).open("rb") as h:
  for b in iter(lambda:h.read(8*1024*1024),b""):d.update(b)
 return d.hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument("--analysis",type=Path,required=True);p.add_argument("--output",type=Path,required=True);a=p.parse_args();
 if a.output.exists():raise FileExistsError(a.output)
 a.output.mkdir(parents=True);c=[]
 def ck(n,x,d):c.append({"check":n,"status":"PASS" if bool(x) else "FAIL","detail":d})
 r=a.analysis;s=json.loads((r/"PHASE20_SPP1_VASCULAR_STATUS_v1.json").read_text());ck("formal_status",s.get("status")=="PASS_ANALYSIS_COMPLETE",str(s.get("status")));parts=list(r.rglob("*.part"));ck("no_part_files",not parts,f"count={len(parts)}");m=pd.read_csv(r/"OUTPUT_MANIFEST_v1.tsv",sep="\t");ok=True
 for x in m.itertuples(index=False):
  f=r/x.file;ok &= f.is_file() and f.stat().st_size==int(x.size_bytes) and sha(f)==x.sha256
 ck("manifest",ok,f"files={len(m)}");t=r/"tables";spot=pd.read_csv(t/"spot_vascular_immune_scores_v1.tsv.gz",sep="\t");sec=pd.read_csv(t/"section_statistics_v1.tsv",sep="\t");pri=pd.read_csv(t/"primary_tjrcc_section_statistics_v1.tsv",sep="\t");agg=pd.read_csv(t/"primary_across_section_inference_v1.tsv",sep="\t");loo=pd.read_csv(t/"primary_leave_one_section_out_v1.tsv",sep="\t");niche=pd.read_csv(t/"niche_vascular_immune_summary_v1.tsv",sep="\t");qc=pd.read_csv(t/"signature_gene_coverage_v1.tsv",sep="\t")
 ck("identity",len(spot)==41688 and spot.sample.nunique()==14 and not spot.duplicated(["sample","spot_id"]).any(),f"spots={len(spot)}");v=spot.drop(columns=["cohort","sample","spot_id","niche"]).to_numpy(float);ck("finite_ranges",np.isfinite(v).all() and ((spot.spp1_polarization>=0)&(spot.spp1_polarization<=1)).all() and ((spot.spp1_vascular_niche>=0)&(spot.spp1_vascular_niche<=1)).all(),"finite and bounded");ck("niches",len(niche)==42 and set(niche.niche)==set(spot.niche),f"n={len(niche)}");ck("signature_coverage",len(qc)==42 and (qc.present>=qc.defined).all(),f"rows={len(qc)}")
 exp=sec[(sec.cohort=="TJRCC")&(sec.role=="primary")].sort_values(["relation","sample"]).reset_index(drop=True);got=pri.sort_values(["relation","sample"]).reset_index(drop=True);ck("primary_conservation",len(got)==60 and exp.equals(got),f"rows={len(got)}");ck("no_spot_p",not any("p_value" in x for x in spot.columns) and not any("p_value" in x for x in sec.columns),"section stats only");med=pri.groupby("relation").spearman_rho.median();j=agg.set_index("relation");ck("aggregate",set(agg.relation)==set(PRIMARY) and np.allclose(j.loc[med.index,"median_rho"],med),"5 relations");rec=[]
 for relation in PRIMARY:
  sub=pri[pri.relation==relation]
  for sample in sorted(sub["sample"]):rec.append((relation,sample,float(np.median(sub[sub["sample"]!=sample].spearman_rho))))
 z=pd.DataFrame(rec,columns=["relation","left_out_sample","recomputed"]);q=loo.merge(z,on=["relation","left_out_sample"]);ck("loo",len(q)==60 and np.allclose(q.median_rho,q.recomputed),f"rows={len(q)}");ck("context_secondary",sec[sec.relation.str.contains("TYMP|MTA")].role.eq("secondary_context").all(),"TYMP/MTA not primary");ck("cohort_separation",spot.loc[spot.cohort=="TJRCC","sample"].nunique()==12 and set(spot.loc[spot.cohort=="GLMF","sample"])=={"GLMF1","GLMF2"},"12+2")
 f=pd.DataFrame(c);part=a.output/"PHASE20_INDEPENDENT_VALIDATION_STATUS_v1.tsv.part";f.to_csv(part,sep="\t",index=False);os.replace(part,a.output/"PHASE20_INDEPENDENT_VALIDATION_STATUS_v1.tsv");passed=f.status.eq("PASS").all();summary={"schema_version":"phase20_independent_validation_v1","status":"PASS" if passed else "FAIL","validated_at":datetime.now().astimezone().isoformat(timespec="seconds"),"checks_passed":int(f.status.eq("PASS").sum()),"checks_total":len(f),"sections":14,"spots":len(spot),"primary_rows":len(pri)};(a.output/"PHASE20_INDEPENDENT_VALIDATION_SUMMARY_v1.json").write_text(json.dumps(summary,indent=2)+"\n");print(json.dumps(summary,indent=2));raise SystemExit(0 if passed else 1)
if __name__=="__main__":main()
