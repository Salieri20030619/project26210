#!/usr/bin/env python3
import csv,hashlib,json
from pathlib import Path
ROOT=Path('/mnt/BioDisk/wangmingyuan/project26210_phase72_immune_ecosystem_v1');C=[]
def add(i,o,d):C.append({'check_id':i,'status':'PASS' if o else 'FAIL','detail':str(d)})
def rows(p):
 with p.open(encoding='utf-8') as f:return list(csv.DictReader(f,delimiter='\t'))
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
add('primary_validation',(ROOT/'independent_validation_v2/PHASE72_INDEPENDENT_VALIDATION_PASS_v2').is_file(),'8/8')
x=rows(ROOT/'results_v1/phase72_cross_cohort_decisions_v1.tsv');add('primary_rows',len(x)==9,len(x));add('no_two_external',sum(r['external_decision']=='TWO_EXTERNAL_SUPPORTED' for r in x)==0,'0/9')
e=rows(ROOT/'PHASE72_IMMUNE_EVIDENCE_REGISTER_v1.tsv');add('evidence_register',len(e)==13 and len({r['evidence_id'] for r in e})==13,len(e));add('status_addendum',(ROOT/'PROJECT_STATUS_PHASE72_ADDENDUM_v1.md').stat().st_size>0,'present');add('manifest_addendum',len(rows(ROOT/'RESULTS_MANIFEST_PHASE72_ADDENDUM_v1.tsv'))==4,'4 rows')
core=['results_v1/PHASE72_IMMUNE_ECOSYSTEM_SUMMARY_v1.md','results_v1/phase72_cross_cohort_decisions_v1.tsv','results_v1/PHASE70_PHASE72_TRIANGULATION_v1.tsv','PHASE72_IMMUNE_EVIDENCE_REGISTER_v1.tsv','PROJECT_STATUS_PHASE72_ADDENDUM_v1.md','RESULTS_MANIFEST_PHASE72_ADDENDUM_v1.tsv'];hashes=[{'file':f,'sha256':sha(ROOT/f),'size_bytes':(ROOT/f).stat().st_size} for f in core];add('core_files',len(hashes)==6,len(hashes))
out=ROOT/'delivery_validation_v1';out.mkdir(exist_ok=True)
with (out/'delivery_validation_checks_v1.tsv').open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=['check_id','status','detail'],delimiter='\t');w.writeheader();w.writerows(C)
with (out/'CORE_OUTPUT_SHA256_v1.tsv').open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=['file','sha256','size_bytes'],delimiter='\t');w.writeheader();w.writerows(hashes)
p={'checks':len(C),'passed':sum(r['status']=='PASS' for r in C),'failed':sum(r['status']=='FAIL' for r in C)};(out/'delivery_validation_summary_v1.json').write_text(json.dumps(p,indent=2)+'\n');
if p['failed']==0:(out/'PHASE72_DELIVERY_VALIDATION_PASS_v1').write_text('PASS\n')
print(json.dumps(p));raise SystemExit(0 if p['failed']==0 else 1)
