#!/usr/bin/env python3
import csv,hashlib,json
from pathlib import Path
ROOT=Path('/mnt/BioDisk/wangmingyuan/project26210_phase71_immune_advanced_v1');C=[]
def add(i,o,d):C.append({'check_id':i,'status':'PASS' if o else 'FAIL','detail':str(d)})
def tsv(p):
 with p.open(encoding='utf-8') as f:return list(csv.DictReader(f,delimiter='\t'))
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
add('primary_validation',(ROOT/'independent_validation_v2/PHASE71_INDEPENDENT_VALIDATION_PASS_v2').is_file(),'9/9')
add('posthoc_validation',(ROOT/'posthoc_independent_validation_v1/POSTHOC_INDEPENDENT_VALIDATION_PASS_v1').is_file(),'6/6')
bridge=tsv(ROOT/'results_v2/phase71_primary_bridge_decision_v1.tsv');add('bridge_row',len(bridge)==1 and bridge[0]['decision']=='TRANSCRIPTIONAL_BRIDGE_NOT_FULLY_SUPPORTED',bridge)
ev=tsv(ROOT/'PHASE71_IMMUNE_EVIDENCE_REGISTER_v1.tsv');add('evidence_register',len(ev)==13 and len({r['evidence_id'] for r in ev})==13,len(ev))
add('status_addendum',(ROOT/'PROJECT_STATUS_PHASE71_ADDENDUM_v1.md').stat().st_size>0,'present')
add('manifest_addendum',len(tsv(ROOT/'RESULTS_MANIFEST_PHASE71_ADDENDUM_v1.tsv'))==4,'4 rows')
core=['results_v2/PHASE71_IMMUNE_ADVANCED_SUMMARY_v1.md','results_v2/phase71_cross_cohort_decisions_v1.tsv','posthoc_state_results_v1/PHASE71_POSTHOC_STATE_AUDIT_SUMMARY_v1.md','PHASE71_IMMUNE_EVIDENCE_REGISTER_v1.tsv','PROJECT_STATUS_PHASE71_ADDENDUM_v1.md','RESULTS_MANIFEST_PHASE71_ADDENDUM_v1.tsv']
hashes=[{'file':x,'sha256':sha(ROOT/x),'size_bytes':(ROOT/x).stat().st_size} for x in core];add('core_files',all((ROOT/x).is_file() for x in core),len(core))
out=ROOT/'delivery_validation_v1';out.mkdir(exist_ok=True)
with (out/'delivery_validation_checks_v1.tsv').open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=['check_id','status','detail'],delimiter='\t');w.writeheader();w.writerows(C)
with (out/'CORE_OUTPUT_SHA256_v1.tsv').open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=['file','sha256','size_bytes'],delimiter='\t');w.writeheader();w.writerows(hashes)
x={'checks':len(C),'passed':sum(r['status']=='PASS' for r in C),'failed':sum(r['status']=='FAIL' for r in C)};(out/'delivery_validation_summary_v1.json').write_text(json.dumps(x,indent=2)+'\n')
if x['failed']==0:(out/'PHASE71_DELIVERY_VALIDATION_PASS_v1').write_text('PASS\n')
print(json.dumps(x));raise SystemExit(0 if x['failed']==0 else 1)
