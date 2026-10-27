#!/usr/bin/env python3
import hashlib,json,sys
from datetime import datetime
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parent; SYNC=ROOT/"server_results_sync_v2"; V=ROOT/"independent_validation_local_v2"/"status.json"; OUT=ROOT/"terminal_validation_local_v1"; c=[]
def ck(n,o,d):c.append({"check":n,"passed":bool(o),"detail":str(d)})
def sha(p):
 h=hashlib.sha256()
 with open(p,"rb") as f:
  for x in iter(lambda:f.read(1048576),b""):h.update(x)
 return h.hexdigest()
def main():
 v=json.loads(V.read_text()); pre=json.loads((SYNC/"preflight_final_v3"/"status.json").read_text()); a=json.loads((SYNC/"analysis_v1"/"status.json").read_text()); f=json.loads((SYNC/"s04_final_v3"/"status.json").read_text())
 ck("independent PASS",v["status"].startswith("PASS"),v["status"]); ck("independent 23/23",v["checks_passed"]==v["checks_total"]==23,"23/23"); ck("preflight 11/11",pre["checks_passed"]==pre["checks_total"]==11,"11/11"); ck("S02 analysis PASS",a["status"].startswith("PASS"),a["status"]); ck("S04 final PASS",f["status"].startswith("PASS"),f["status"]); ck("S02 8 inference rows",len(pd.read_csv(SYNC/"analysis_v1"/"tables"/"s02_tjrcc12_inference_v1.tsv",sep="\t"))==8,"8"); ck("S04 6 pair rows",len(pd.read_csv(SYNC/"s04_reanalysis_v2"/"tables"/"s04_tjrcc12_inference_block_mean_v2.tsv",sep="\t"))==6,"6"); current={str(p.relative_to(ROOT)):sha(p) for p in SYNC.rglob("*") if p.is_file()}; ck("validated hashes unchanged",current==v["validated_files"],len(current)); ck("failed attempts preserved",(ROOT/"independent_validation_local_v1"/"status.json").exists() and (SYNC/"s04_reanalysis_v2"/"status.json").exists(),"preserved"); ck("three registries present",all((ROOT/x).exists() for x in ["PHASE37_S02_S04_ANALYSIS_REGISTRY_v1.md","PHASE37_S02_S04_ANALYSIS_REGISTRY_v2.md","PHASE37_S04_ANALYSIS_REGISTRY_v3.md"]),"v1-v3"); figs=[p for p in ROOT.rglob("*") if p.suffix.lower() in {".png",".jpg",".jpeg",".pdf",".svg",".tif",".tiff"}]; ck("no figures",not figs,len(figs))
 OUT.mkdir(parents=False,exist_ok=False); d=pd.DataFrame(c); d.to_csv(OUT/"phase37_terminal_checks_v1.tsv",sep="\t",index=False); n=int(d.passed.sum()); payload={"status":"PASS_PHASE37_S02_S04_TERMINAL_VALIDATION" if n==len(d) else "FAIL_PHASE37_S02_S04_TERMINAL_VALIDATION","completed_at":datetime.now().astimezone().isoformat(timespec="seconds"),"checks_passed":n,"checks_total":len(d),"component_checks":"23/23","figures_created":False}; (OUT/"status.json").write_text(json.dumps(payload,indent=2),encoding="utf-8"); print(json.dumps(payload,indent=2)); return 0 if n==len(d) else 1
if __name__=="__main__":sys.exit(main())
