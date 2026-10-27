from __future__ import annotations
import csv,hashlib,json,re
from pathlib import Path
S=Path(r"C:\Users\31937\project26210_phase108_state_update_v1");B=S/'preupdate_backup'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for x in iter(lambda:f.read(1024*1024),b''):h.update(x)
 return h.hexdigest()
def read(p):return p.read_text(encoding='utf-8-sig')
def tsv(p):
 with p.open('r',encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
checks=[]
def check(n,p,d=''):checks.append({'check':n,'passed':bool(p),'detail':str(d)})
req=['PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv','CODEX_HANDOFF.md']
for n in req:check('required_'+n,(S/n).is_file() and (S/n).stat().st_size>0)
markers={'PROJECT_STATUS.md':'## Phase 107-108 evidence escalation and unified client architecture (2026-08-31)','DECISIONS.md':'## 2026-08-31 - Phase 107-108 unified-chain decision','CODEX_HANDOFF.md':'## Phase 107-108 current unified-chain state (2026-08-31)'}
for n,m in markers.items():
 old=read(B/n);new=read(S/n);check('append_'+n,new.startswith(old));check('marker_'+n,new.count(m)==1)
old=tsv(B/'RESULTS_MANIFEST.tsv');new=tsv(S/'RESULTS_MANIFEST.tsv');check('append18',len(new)==len(old)+18,f'{len(old)}->{len(new)}');check('old_preserved',new[:len(old)]==old)
ids=[r['result_id'] for r in new];check('unique_ids',len(ids)==len(set(ids)))
groups={'phase107_evidence_escalation_audit':8,'phase108_unified_chain_architecture':7,'phase108_project_state_update':3}
selected=[]
for stage,count in groups.items():
 rows=[r for r in new if r['stage']==stage];check('count_'+stage,len(rows)==count,len(rows));check('pass_'+stage,all(r['validation_status']=='PASS' for r in rows));selected+=rows
for r in selected:
 p=Path(r['absolute_path'])
 if r['stage']=='phase108_project_state_update':p=S/Path(r['absolute_path']).name
 check('path_'+r['result_id'],p.is_file(),p);m=re.search(r'SHA256 ([0-9A-Fa-f]{64})',r['validation_evidence'])
 if m and p.is_file():check('hash_'+r['result_id'],sha(p).lower()==m.group(1).lower())
v107=json.loads(Path(r"D:\project26210_external_joint_v1\phase107_evidence_escalation_audit_local_v1\PHASE107_EVIDENCE_ESCALATION_VALIDATION_v1.json").read_text(encoding='utf-8'));v108=json.loads(Path(r"D:\project26210_external_joint_v1\phase108_unified_chain_architecture_local_v1\PHASE108_UNIFIED_CHAIN_VALIDATION_v1.json").read_text(encoding='utf-8'))
check('v107',v107.get('status')=='PASS_PHASE107_EVIDENCE_ESCALATION_VALIDATION' and v107.get('checks_passed')==26);check('v108',v108.get('status')=='PASS_PHASE108_UNIFIED_CHAIN_VALIDATION' and v108.get('checks_passed')==22)
status=read(S/'PROJECT_STATUS.md');dec=read(S/'DECISIONS.md');hand=read(S/'CODEX_HANDOFF.md')
check('no_parts','retires the original Part1-4 divisions' in status);check('one_chain','one continuous biological chain' in status);check('no_pool','Do not pool incompatible datasets' in dec);check('negative_edges','Retain negative, contradictory and non-estimable edges' in dec);check('handoff_chain','Narrative chain:' in hand);check('no_full_claim','Do not claim a fully validated chain' in hand);check('no_job','No Phase107/108 background job remains' in hand)
failed=[c for c in checks if not c['passed']];payload={'status':'PASS_PHASE108_COMPLETION_GATE' if not failed else 'FAIL_PHASE108_COMPLETION_GATE','checks_passed':len(checks)-len(failed),'checks_total':len(checks),'failed_checks':failed,'key_file_sha256':{n:sha(S/n) for n in req}}
(S/'STATE_UPDATE_VALIDATION_v1.json').write_text(json.dumps(payload,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');print(json.dumps(payload,ensure_ascii=False,indent=2));raise SystemExit(0 if not failed else 1)
