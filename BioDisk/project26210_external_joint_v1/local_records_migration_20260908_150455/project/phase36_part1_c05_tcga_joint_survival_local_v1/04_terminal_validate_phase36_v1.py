#!/usr/bin/env python3
import hashlib, json, sys
from datetime import datetime
from pathlib import Path
import pandas as pd

ROOT=Path(__file__).resolve().parent; SRC=ROOT/"server_results_sync_v1"/"analysis_v1"; OUT=ROOT/"terminal_validation_local_v1"; V=ROOT/"independent_validation_local_v2"/"status.json"
checks=[]
def ck(name,ok,detail): checks.append({"check":name,"passed":bool(ok),"detail":str(detail)})
def sha(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(1048576),b""): h.update(b)
    return h.hexdigest()
def main():
    v=json.loads(V.read_text(encoding="utf-8")); a=json.loads((SRC/"status.json").read_text(encoding="utf-8"))
    ck("analysis status PASS",a["status"].startswith("PASS"),a["status"])
    ck("independent validation PASS",v["status"]=="PASS_PHASE36_C05_INDEPENDENT_VALIDATION",v["status"])
    ck("independent checks complete",v["checks_passed"]==v["checks_total"]==21,"21/21")
    ck("median-only frozen cutoffs",a.get("median_cutoffs_only") is True,"median only")
    ck("no prespecified direction",a.get("direction_prespecified") is False,"two-sided associations")
    cohort=pd.read_csv(SRC/"c05_tcga_analysis_cohort_v1.tsv",sep="\t"); lr=pd.read_csv(SRC/"c05_joint_group_logrank_summary_v1.tsv",sep="\t"); cox=pd.read_csv(SRC/"c05_joint_group_cox_models_v1.tsv",sep="\t"); inter=pd.read_csv(SRC/"c05_continuous_interaction_cox_v1.tsv",sep="\t")
    ck("cohort identity",len(cohort)==510 and cohort.patient.nunique()==510,cohort.shape)
    ck("log-rank grid",len(lr)==12 and lr.receptor_axis.nunique()==6 and lr.endpoint.nunique()==2,lr.shape)
    ck("group Cox grid populated",len(cox)==108 and cox.receptor_axis.nunique()==6 and cox.endpoint.nunique()==2,cox.shape)
    ck("interaction Cox grid populated",len(inter)==96 and inter.receptor_axis.nunique()==6 and inter.endpoint.nunique()==2,inter.shape)
    current={str(p.relative_to(ROOT)):sha(p) for p in sorted(SRC.iterdir()) if p.is_file()}
    ck("validated file hashes unchanged",current==v["validated_files"],f"{len(current)} files")
    figs=[p for p in ROOT.rglob("*") if p.suffix.lower() in {".png",".jpg",".jpeg",".pdf",".svg",".tif",".tiff"}]
    ck("no figures generated",not figs,len(figs))
    ck("failed v1 validator preserved",(ROOT/"independent_validation_local_v1"/"status.json").exists(),"rounding-contract failure retained")
    ck("analysis registry present",(ROOT/"PHASE36_C05_ANALYSIS_REGISTRY_v1.md").exists(),"registry")
    OUT.mkdir(parents=False,exist_ok=False); df=pd.DataFrame(checks); df.to_csv(OUT/"phase36_terminal_checks_v1.tsv",sep="\t",index=False); pd.DataFrame([{"path":k,"sha256":z} for k,z in current.items()]).to_csv(OUT/"phase36_terminal_file_manifest_v1.tsv",sep="\t",index=False)
    n=int(df.passed.sum()); payload={"status":"PASS_PHASE36_C05_TERMINAL_VALIDATION" if n==len(df) else "FAIL_PHASE36_C05_TERMINAL_VALIDATION","completed_at":datetime.now().astimezone().isoformat(timespec="seconds"),"checks_passed":n,"checks_total":len(df),"component_checks":"21/21","figures_created":False}; (OUT/"status.json").write_text(json.dumps(payload,indent=2),encoding="utf-8"); print(json.dumps(payload,indent=2)); return 0 if n==len(df) else 1
if __name__=="__main__": sys.exit(main())
