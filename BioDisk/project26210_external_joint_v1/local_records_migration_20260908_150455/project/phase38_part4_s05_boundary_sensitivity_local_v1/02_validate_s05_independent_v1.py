#!/usr/bin/env python3
import hashlib,json,sys
from datetime import datetime
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import spearmanr
ROOT=Path(__file__).resolve().parent; SRC=ROOT/"server_results_sync_v2"/"analysis_v1"; OUT=ROOT/"independent_validation_local_v1"; checks=[]
def ck(n,o,d):checks.append({"check":n,"passed":bool(o),"detail":str(d)})
def close(a,b):return np.allclose(np.asarray(a,float),np.asarray(b,float),rtol=1e-8,atol=1e-10,equal_nan=True)
def sha(p):
 h=hashlib.sha256()
 with open(p,"rb") as f:
  for x in iter(lambda:f.read(1048576),b""):h.update(x)
 return h.hexdigest()
def main():
 st=json.loads((SRC/"status.json").read_text()); t=SRC/"tables"; spot=pd.read_csv(t/"s05_draft_signed_distance_spot_scores_v1.tsv.gz",sep="\t"); prof=pd.read_csv(t/"s05_draft_boundary_100um_profiles_v1.tsv",sep="\t"); sm=pd.read_csv(t/"s05_draft_boundary_gradient_summary_v1.tsv",sep="\t"); qc=pd.read_csv(t/"s05_draft_boundary_qc_v1.tsv",sep="\t")
 ck("analysis status PASS draft sensitivity",st["status"].startswith("PASS_S05_AI_DRAFT"),st["status"]); ck("formal pathology scope is zero",st["formal_pathology_boundary_sections"]==0 and not qc.formal_pathology_ready.astype(bool).any(),"0"); ck("two draft sections",st["draft_boundary_sections"]==2 and set(qc["sample"])=={"GLMF1","GLMF2"},list(qc["sample"])); ck("review status explicitly not pathologist",set(qc.review_status)=={"not_pathologist_reviewed"},set(qc.review_status)); ck("both regions present per section",(qc.n_tumor_spots>0).all() and (qc.n_nat_spots>0).all(),qc[["sample","n_tumor_spots","n_nat_spots"]].to_dict("records")); ck("shared physical interfaces present",(qc.n_shared_interface_segments>0).all(),qc.n_shared_interface_segments.tolist()); ck("micron calibration plausible",qc.um_per_hires_pixel.between(3.9,4.0).all(),qc.um_per_hires_pixel.tolist()); ck("spot counts conserved",len(spot)==qc.n_spots_within_1000um.sum()==1932,len(spot)); ck("signed distance range",np.isfinite(spot.signed_distance_um).all() and (spot.signed_distance_um.abs()<=1000).all(),spot.signed_distance_um.abs().max()); ck("distance sign matches region",((spot.draft_region=="tumor")==spot.signed_distance_um.ge(0)).all(),"positive=tumor"); ck("no duplicated section-spot",not spot.duplicated(["sample","spot_id"]).any(),spot.shape); ck("profile bin width and range",(prof.distance_bin_lower_um%100==0).all() and prof.distance_bin_lower_um.between(-1000,900).all(),prof.distance_bin_lower_um.min());
 merged=[]
 for (sample,var),f in prof.groupby(["sample","variable"],sort=True):
  z=spot[spot["sample"].eq(sample)].groupby("distance_bin_lower_um")[var].agg(["count","median"]).reset_index(); f=f.sort_values("distance_bin_lower_um"); merged.append(len(f)==len(z) and np.array_equal(f.n_spots,z["count"]) and close(f.median_value,z["median"]))
 ck("100-um profiles independently reproduced",all(merged),f"{sum(merged)}/{len(merged)}")
 rhos=[]
 for _,r in sm.iterrows():
  z=prof[(prof["sample"]==r["sample"])&(prof.variable==r.variable)].sort_values("distance_bin_lower_um"); rho=float(spearmanr(z.distance_bin_lower_um,z.median_value).statistic) if len(z)>=5 and np.ptp(z.median_value)>0 else np.nan; rhos.append(rho)
 ck("gradient rhos independently reproduced",close(sm.spearman_rho_across_bin_medians,rhos),f"{len(rhos)} variables"); ck("no inferential P values",not any("p_" in c.lower() or c.lower().endswith("_p") for c in sm.columns),list(sm.columns)); ck("no prespecified direction and no figures",not st["direction_prespecified"] and not st["figures_created"],"contract"); figs=[p for p in ROOT.rglob("*") if p.suffix.lower() in {".png",".jpg",".jpeg",".pdf",".svg",".tif",".tiff"}]; ck("filesystem contains no figures",not figs,len(figs))
 OUT.mkdir(parents=False,exist_ok=False); df=pd.DataFrame(checks); df.to_csv(OUT/"phase38_s05_independent_validation_checks_v1.tsv",sep="\t",index=False); n=int(df.passed.sum()); files={str(p.relative_to(ROOT)):sha(p) for p in SRC.rglob("*") if p.is_file()}; payload={"status":"PASS_PHASE38_S05_DRAFT_SENSITIVITY_INDEPENDENT_VALIDATION" if n==len(df) else "FAIL_PHASE38_S05_DRAFT_SENSITIVITY_INDEPENDENT_VALIDATION","completed_at":datetime.now().astimezone().isoformat(timespec="seconds"),"checks_passed":n,"checks_total":len(df),"formal_pathology_complete":False,"validated_files":files,"figures_created":False}; (OUT/"status.json").write_text(json.dumps(payload,indent=2),encoding="utf-8"); print(json.dumps(payload,indent=2)); return 0 if n==len(df) else 1
if __name__=="__main__":sys.exit(main())
