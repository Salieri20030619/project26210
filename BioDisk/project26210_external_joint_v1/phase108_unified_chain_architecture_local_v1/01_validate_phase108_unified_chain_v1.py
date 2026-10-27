from __future__ import annotations
import csv, hashlib, json
from pathlib import Path
ROOT=Path(r"C:\Users\31937\project26210_phase108_unified_chain_architecture_v1")
def tsv(name):
    with (ROOT/name).open("r",encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
    return h.hexdigest()
checks=[]
def check(name,passed,detail=""):checks.append({"check":name,"passed":bool(passed),"detail":str(detail)})
required=["PHASE108_UNIFIED_CHAIN_RULE_v1.md","UNIFIED_CHAIN_NODE_MATRIX_v1.tsv","CLIENT_REPORT_ARCHITECTURE_v4.md","UNIFIED_FIGURE_PLAN_v1.tsv"]
for n in required:check(f"required_{n}",(ROOT/n).is_file() and (ROOT/n).stat().st_size>0)
nodes=tsv("UNIFIED_CHAIN_NODE_MATRIX_v1.tsv");check("nodes_11",len(nodes)==11,len(nodes));check("node_ids_unique",len({r['node_id'] for r in nodes})==11)
check("tumour_origin",any(r['node_id']=='U01' and r['current_evidence_state']=='SUPPORTED_WITH_HETEROGENEITY' for r in nodes))
check("direct_mta_not_supported",any(r['node_id']=='U03' and r['current_evidence_state']=='CURRENT_DATA_DO_NOT_SUPPORT_OR_NOT_EVALUABLE' for r in nodes))
check("mta_spp1_pending",any(r['node_id']=='U06' and r['current_evidence_state']=='NOT_CLOSED_PENDING_CLIENT' for r in nodes))
check("spp1_vascular_specific",any(r['node_id']=='U08' and r['current_evidence_state']=='DATASET_SPECIFIC_SPATIAL_SUPPORT' for r in nodes))
check("clinical_unsupported",any(r['node_id']=='U10' and r['current_evidence_state']=='CURRENT_DATA_DO_NOT_SUPPORT' for r in nodes))
figs=tsv("UNIFIED_FIGURE_PLAN_v1.tsv");check("figures_6",len(figs)==6,len(figs));check("figure_ids_unique",len({r['figure_id'] for r in figs})==6)
check("client_slots_no_data",any(r['figure_id']=='F06' and 'without project data' in r['client_validation_interface'] for r in figs))
rule=(ROOT/"PHASE108_UNIFIED_CHAIN_RULE_v1.md").read_text(encoding="utf-8");arch=(ROOT/"CLIENT_REPORT_ARCHITECTURE_v4.md").read_text(encoding="utf-8")
check("no_part_divisions","original Part 1–4 divisions are retired" in rule)
check("one_chain","one continuous biological evidence chain" in rule)
check("no_statistical_mixing","does not authorize statistical mixing" in rule)
check("negative_edges_retained","negative, contradictory and non-estimable edges remain" in rule)
check("client_reply_no_parts","最终报告不再保留原来一到四部分的分割" in arch)
check("architecture_seven_sections",all(f"### {i}." in arch for i in range(1,8)))
check("conclusion_not_full_validation","should not say that the complete chain is validated" in arch)
check("v4_supersession","supersedes the client-facing organisation" in arch)
failed=[c for c in checks if not c['passed']];payload={"status":"PASS_PHASE108_UNIFIED_CHAIN_VALIDATION" if not failed else "FAIL_PHASE108_UNIFIED_CHAIN_VALIDATION","checks_passed":len(checks)-len(failed),"checks_total":len(checks),"failed_checks":failed}
(ROOT/"PHASE108_UNIFIED_CHAIN_VALIDATION_v1.json").write_text(json.dumps(payload,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
rows=[]
for p in sorted(x for x in ROOT.iterdir() if x.is_file() and x.name!="FILE_MANIFEST_SHA256_v1.tsv"):rows.append({"file":p.name,"size_bytes":p.stat().st_size,"sha256":sha(p)})
with (ROOT/"FILE_MANIFEST_SHA256_v1.tsv").open("w",encoding="utf-8",newline="") as f:
    w=csv.DictWriter(f,fieldnames=["file","size_bytes","sha256"],delimiter="\t",lineterminator="\n");w.writeheader();w.writerows(rows)
print(json.dumps(payload,ensure_ascii=False,indent=2));raise SystemExit(0 if not failed else 1)
