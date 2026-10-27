#!/usr/bin/env python3
import hashlib,json,sys
from datetime import datetime
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parent; SRC=ROOT/"server_results_sync_v2"/"analysis_v1"; V=ROOT/"independent_validation_local_v1"/"status.json"; OUT=ROOT/"terminal_validation_local_v1"; c=[]
def ck(n,o,d):c.append({"check":n,"passed":bool(o),"detail":str(d)})
def sha(p):
 h=hashlib.sha256()
 with open(p,"rb") as f:
  for x in iter(lambda:f.read(1048576),b""):h.update(x)
 return h.hexdigest()
def main():
 v=json.loads(V.read_text()); a=json.loads((SRC/"status.json").read_text()); qc=pd.read_csv(SRC/"tables"/"s05_draft_boundary_qc_v1.tsv",sep="\t")
 ck("analysis PASS draft sensitivity",a["status"].startswith("PASS_S05_AI_DRAFT"),a["status"]); ck("independent PASS",v["status"].startswith("PASS"),v["status"]); ck("independent 17/17",v["checks_passed"]==v["checks_total"]==17,"17/17"); ck("two draft sections",len(qc)==2 and a["draft_boundary_sections"]==2,"2"); ck("formal pathology remains pending",a["formal_pathology_boundary_sections"]==0 and v["formal_pathology_complete"] is False,"0 formal"); ck("both region scopes populated",(qc.n_tumor_spots>0).all() and (qc.n_nat_spots>0).all(),qc.n_spots_within_1000um.sum()); current={str(p.relative_to(ROOT)):sha(p) for p in SRC.rglob("*") if p.is_file()}; ck("validated hashes unchanged",current==v["validated_files"],len(current)); ck("registry present",(ROOT/"PHASE38_S05_ANALYSIS_REGISTRY_v1.md").exists(),"registry"); figs=[p for p in ROOT.rglob("*") if p.suffix.lower() in {".png",".jpg",".jpeg",".pdf",".svg",".tif",".tiff"}]; ck("no figures",not figs,len(figs))
 OUT.mkdir(parents=False,exist_ok=False); d=pd.DataFrame(c); d.to_csv(OUT/"phase38_terminal_checks_v1.tsv",sep="\t",index=False); n=int(d.passed.sum()); payload={"status":"PASS_PHASE38_S05_DRAFT_SENSITIVITY_TERMINAL_VALIDATION_FORMAL_PATHOLOGY_PENDING" if n==len(d) else "FAIL_PHASE38_S05_TERMINAL_VALIDATION","completed_at":datetime.now().astimezone().isoformat(timespec="seconds"),"checks_passed":n,"checks_total":len(d),"component_checks":"17/17","formal_pathology_complete":False,"figures_created":False}; (OUT/"status.json").write_text(json.dumps(payload,indent=2),encoding="utf-8"); print(json.dumps(payload,indent=2)); return 0 if n==len(d) else 1
if __name__=="__main__":sys.exit(main())
