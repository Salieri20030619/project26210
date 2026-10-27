#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, math, sys
from datetime import datetime
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

def inside(points,poly):
 x=points[:,0]; y=points[:,1]; p=np.asarray(poly,float); ans=np.zeros(len(points),bool); j=len(p)-1
 for i in range(len(p)):
  xi,yi=p[i]; xj,yj=p[j]; hit=((yi>y)!=(yj>y))&(x<(xj-xi)*(y-yi)/(yj-yi+1e-300)+xi); ans^=hit; j=i
 return ans
def keyseg(a,b):
 a=tuple(map(float,a)); b=tuple(map(float,b)); return tuple(sorted([a,b]))
def shared_segments(tumor,nats):
 def segs(poly): return {keyseg(poly[i],poly[(i+1)%len(poly)]) for i in range(len(poly))}
 ts=segs(tumor); ns=set().union(*(segs(p) for p in nats)); return list(ts&ns)
def distance_segments(points,segs):
 out=np.full(len(points),np.inf)
 for aa,bb in segs:
  a=np.asarray(aa,float); b=np.asarray(bb,float); ab=b-a; t=np.clip(((points-a)@ab)/(ab@ab),0,1); d=np.sqrt(((points-(a+t[:,None]*ab))**2).sum(axis=1)); out=np.minimum(out,d)
 return out
def atomic_tsv(x,p,gz=False):
 q=Path(str(p)+".part"); x.to_csv(q,sep="\t",index=False,compression="gzip" if gz else None); q.replace(p)
def atomic_json(x,p):
 q=Path(str(p)+".part"); q.write_text(json.dumps(x,indent=2),encoding="utf-8"); q.replace(p)

def main():
 ap=argparse.ArgumentParser(); ap.add_argument("--phase18-spots",type=Path,required=True); ap.add_argument("--phase19-spots",type=Path,required=True); ap.add_argument("--s02-spots",type=Path,required=True); ap.add_argument("--s04-spots",type=Path,required=True); ap.add_argument("--pathology-dir",type=Path,required=True); ap.add_argument("--output",type=Path,required=True); a=ap.parse_args()
 if a.output.exists(): raise FileExistsError(a.output)
 a.output.mkdir(parents=True); tables=a.output/"tables"; tables.mkdir()
 p18=pd.read_csv(a.phase18_spots,sep="\t"); p19=pd.read_csv(a.phase19_spots,sep="\t"); s02=pd.read_csv(a.s02_spots,sep="\t"); s04=pd.read_csv(a.s04_spots,sep="\t")
 spot_rows=[]; profiles=[]; summaries=[]; qc=[]
 for sample in ["GLMF1","GLMF2"]:
  meta=json.loads((a.pathology_dir/f"{sample}_pathology_polygons.json").read_text(encoding="utf-8")); tumor=[z["vertices"] for z in meta["polygons"] if z["label"]=="tumor"][0]; nats=[z["vertices"] for z in meta["polygons"] if z["label"]=="NAT"]; segs=shared_segments(tumor,nats)
  if not segs: raise RuntimeError(f"no shared draft interface for {sample}")
  base=p18[p18["sample"]==sample][["cohort","sample","spot_id","x","y"]].merge(p19[p19["sample"]==sample][["spot_id","tymp_exposure"]],on="spot_id",validate="one_to_one"); scale=float(meta["spatial_scale_factor"]); pts=base[["x","y"]].to_numpy(float)*scale; it=inside(pts,tumor); inn=np.column_stack([inside(pts,p) for p in nats]).any(axis=1); unique=it^inn; region=np.where(it&~inn,"tumor",np.where(inn&~it,"NAT","excluded")); dist=distance_segments(pts,segs); um_per_hires=(float(meta["spot_diameter_um"])/float(meta["spot_diameter_fullres_px"]))/scale; signed=np.where(region=="tumor",dist,-dist)*um_per_hires
  base["draft_region"]=region; base["signed_distance_um"]=signed; base=base[unique].copy(); base=base[np.abs(base.signed_distance_um)<=1000].copy()
  w02=s02[s02["sample"]==sample].pivot(index="spot_id",columns=["source","pathway"],values="spatial_flux_score"); w02.columns=["S02__"+x+"__"+y.replace(" ","_") for x,y in w02.columns]; w04=s04[s04["sample"]==sample].pivot(index="spot_id",columns="pair",values="spatial_lr_score"); w04.columns=["S04__"+x for x in w04.columns]; base=base.set_index("spot_id").join(w02).join(w04).reset_index(); variables=["tymp_exposure"]+list(w02.columns)+list(w04.columns); base["distance_bin_lower_um"]=(np.floor(base.signed_distance_um/100)*100).astype(int)
  for var in variables:
   z=base.groupby("distance_bin_lower_um",as_index=False).agg(n_spots=(var,"count"),median_value=(var,"median")); z["sample"]=sample; z["variable"]=var; profiles.append(z)
   if len(z)>=5 and np.ptp(z.median_value)>0: rho=float(spearmanr(z.distance_bin_lower_um,z.median_value).statistic)
   else: rho=math.nan
   summaries.append({"sample":sample,"variable":var,"n_distance_bins":len(z),"n_spots":int(z.n_spots.sum()),"spearman_rho_across_bin_medians":rho,"direction_toward_tumor":("increasing" if rho>0 else "decreasing" if rho<0 else "unevaluable")})
  spot_rows.append(base); qc.append({"sample":sample,"review_status":meta.get("review_status"),"annotation_source":meta.get("annotation_source"),"n_shared_interface_segments":len(segs),"um_per_hires_pixel":um_per_hires,"n_tumor_spots":int((base.draft_region=="tumor").sum()),"n_nat_spots":int((base.draft_region=="NAT").sum()),"n_spots_within_1000um":len(base),"formal_pathology_ready":False})
 atomic_tsv(pd.concat(spot_rows,ignore_index=True),tables/"s05_draft_signed_distance_spot_scores_v1.tsv.gz",True); atomic_tsv(pd.concat(profiles,ignore_index=True),tables/"s05_draft_boundary_100um_profiles_v1.tsv"); atomic_tsv(pd.DataFrame(summaries),tables/"s05_draft_boundary_gradient_summary_v1.tsv"); atomic_tsv(pd.DataFrame(qc),tables/"s05_draft_boundary_qc_v1.tsv")
 ok=len(qc)==2 and all(x["n_tumor_spots"]>0 and x["n_nat_spots"]>0 and x["n_shared_interface_segments"]>0 for x in qc)
 payload={"status":"PASS_S05_AI_DRAFT_BOUNDARY_SENSITIVITY_COMPLETE_FORMAL_PATHOLOGY_PENDING" if ok else "FAIL_S05_DRAFT_BOUNDARY_SENSITIVITY","completed_at":datetime.now().astimezone().isoformat(timespec="seconds"),"n_sections":2,"formal_pathology_boundary_sections":0,"draft_boundary_sections":2,"inference":"descriptive per-section gradients only; no P values","direction_prespecified":False,"figures_created":False}; atomic_json(payload,a.output/"status.json"); print(json.dumps(payload,indent=2)); return 0 if ok else 1
if __name__=="__main__": sys.exit(main())
