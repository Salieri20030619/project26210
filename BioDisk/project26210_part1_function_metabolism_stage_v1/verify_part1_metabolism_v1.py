#!/usr/bin/env python3
import argparse,csv,gzip,hashlib,json
from datetime import datetime
from pathlib import Path
parser=argparse.ArgumentParser(); parser.add_argument("--record-recovery",action="store_true"); args=parser.parse_args()
ROOT=Path("/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/09_metabolism_trajectory/metabolism_v1")
V=ROOT.parent/"metabolism_verification_v1"; V.mkdir(parents=True,exist_ok=True)
checks=[]
def ck(n,p,d): checks.append((n,bool(p),str(d)))
def read(p):
    op=gzip.open if p.suffix==".gz" else open
    with op(p,"rt",encoding="utf-8",newline="") as h:return list(csv.DictReader(h,delimiter="\t"))
req=["RUN_COMPLETE","AUC_COMPLETE","SCFEA_INPUT_COMPLETE","run_metadata_v1.json","file_manifest_v1.tsv","source_data/KEGG_AUCell_all_pathways_v1.tsv.gz","source_data/scFEA_M168_cell_flux.tsv.gz","tables/four_pathway_AUCell_paired_statistics_v1.tsv","tables/scFEA_four_axis_paired_statistics_v1.tsv","tables/scFEA_scope_limitations_v1.tsv"]
for f in req: ck("exists::"+f,(ROOT/f).is_file(),f)
if not all(x[1] for x in checks): raise SystemExit("missing metabolism outputs")
meta=json.loads((ROOT/"run_metadata_v1.json").read_text()); ck("n_cells",meta.get("n_cells")==11950,meta.get("n_cells")); ck("n_modules",meta.get("n_modules")==168,meta.get("n_modules")); ck("n_axes",meta.get("n_axes")==4,meta.get("n_axes"))
for f in ["tables/four_pathway_AUCell_paired_statistics_v1.tsv","tables/scFEA_four_axis_paired_statistics_v1.tsv"]:
    r=read(ROOT/f); ck("eight_contrast_rows::"+f,len(r)==8,len(r)); ck("p_ranges::"+f,all(0<=float(x["paired_wilcox_PValue"])<=1 and 0<=float(x["FDR"])<=1 for x in r),"[0,1]")
for f in ["figures/four_pathway_AUCell_paired_heatmap_v1.png","figures/scFEA_four_axis_paired_heatmap_v1.png","figures/F_scFEA_training_loss.png"]: ck("figure::"+f,(ROOT/f).stat().st_size>10000,(ROOT/f).stat().st_size)
bad=[]
for r in read(ROOT/"file_manifest_v1.tsv"):
    p=ROOT/r["file"]
    h=hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else ""
    if not p.is_file() or p.stat().st_size!=int(r["size_bytes"]) or h!=r["sha256"]:bad.append(r["file"])
ck("manifest_integrity",not bad,",".join(bad) or "all matched")
with (V/"verification_checks_v1.tsv").open("w",encoding="utf-8",newline="") as h:
    w=csv.writer(h,delimiter="\t");w.writerow(["check","status","detail"]);w.writerows((n,"PASS" if p else "FAIL",d) for n,p,d in checks)
summary={"schema_version":"part1_metabolism_verification_v1","n_checks":len(checks),"n_failed":sum(not p for _,p,_ in checks)};(V/"verification_metadata_v1.json").write_text(json.dumps(summary,indent=2))
if summary["n_failed"]:raise SystemExit(summary)
(V/"RUN_COMPLETE").write_text("PASS\n");print(json.dumps(summary,indent=2))
if args.record_recovery:
    now=datetime.now().astimezone().isoformat(timespec="seconds")
    with (ROOT/"background_status_v1.tsv").open("a",encoding="utf-8") as h:
        h.write(f"{now}\tverification_recovery\tPASS\tdynamic files excluded from static manifest\n")
        h.write(f"{now}\tcontroller_recovery\tPASS\tall metabolism stages complete\n")
