#!/usr/bin/env python3
import hashlib, json, sys
from datetime import datetime
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent
SYNC = ROOT / "server_results_sync_v1"
OUT = ROOT / "terminal_validation_local_v1"
V = ROOT / "independent_validation_local_v2" / "status.json"
checks = []
def ck(name, ok, detail): checks.append({"check": name, "passed": bool(ok), "detail": str(detail)})
def sha(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(1048576),b""): h.update(b)
    return h.hexdigest()

def main():
    v=json.loads(V.read_text(encoding="utf-8"))
    ck("independent validation PASS",v["status"]=="PASS_PHASE35_C01_C04_INDEPENDENT_VALIDATION",v["status"])
    ck("independent checks complete",v["checks_passed"]==v["checks_total"]==28,"28/28")
    statuses=[json.loads(p.read_text(encoding="utf-8")) for p in sorted(SYNC.glob("*/status.json"))]
    ck("three component status files",len(statuses)==3,len(statuses))
    ck("all components PASS",all(x["status"].startswith("PASS") for x in statuses),[x["status"] for x in statuses])
    ck("patient-level inference frozen",v.get("patient_is_inferential_unit") is True,"patient")
    ck("two-sided no-direction contract",v.get("two_sided_testing") is True and v.get("direction_prespecified") is False,"two-sided; direction not prespecified")
    s1=pd.read_csv(SYNC/"c01_cellchat_integrated_v1"/"c01_cellchat_endpoint_summary_v1.tsv",sep="\t")
    s3=pd.read_csv(SYNC/"c03_nichenet_v1"/"c03_nichenet_ligand_activity_summary_v1.tsv",sep="\t")
    s4=pd.read_csv(SYNC/"c02_c04_analysis_v1"/"c04_receptor_tymp_f1_summary_v1.tsv",sep="\t")
    ck("C01 endpoint scope",len(s1)==4 and (s1.n_eligible_patients==14).all(),s1.shape)
    ck("C03 ligand-endpoint scope",len(s3)==6 and set(s3.ligand)=={"MST1","ICAM5"},s3.shape)
    ck("C04 outcome-feature scope",len(s4)==16 and s4.outcome.nunique()==2 and s4.receptor_feature.nunique()==8,s4.shape)
    current={str(p.relative_to(ROOT)):sha(p) for p in sorted(SYNC.rglob("*")) if p.is_file()}
    ck("validated file hashes unchanged",current==v["validated_files"],f"{len(current)} files")
    figs=[p for p in ROOT.rglob("*") if p.suffix.lower() in {".png",".jpg",".jpeg",".pdf",".svg",".tif",".tiff"}]
    ck("no figures generated",not figs,len(figs))
    ck("failed v1 validator preserved",(ROOT/"independent_validation_local_v1"/"status.json").exists(),"validator-contract failure retained")
    ck("analysis registry present",(ROOT/"PHASE35_C01_C04_ANALYSIS_REGISTRY_v1.md").exists(),"registry")
    OUT.mkdir(parents=False,exist_ok=False)
    df=pd.DataFrame(checks); df.to_csv(OUT/"phase35_terminal_checks_v1.tsv",sep="\t",index=False)
    pd.DataFrame([{"path":k,"sha256":z} for k,z in current.items()]).to_csv(OUT/"phase35_terminal_file_manifest_v1.tsv",sep="\t",index=False)
    n=int(df.passed.sum()); payload={"status":"PASS_PHASE35_C01_C04_TERMINAL_VALIDATION" if n==len(df) else "FAIL_PHASE35_C01_C04_TERMINAL_VALIDATION","completed_at":datetime.now().astimezone().isoformat(timespec="seconds"),"checks_passed":n,"checks_total":len(df),"component_checks":"28/28","figures_created":False}
    (OUT/"status.json").write_text(json.dumps(payload,indent=2),encoding="utf-8"); print(json.dumps(payload,indent=2)); return 0 if n==len(df) else 1
if __name__=="__main__": sys.exit(main())
