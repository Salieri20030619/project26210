from pathlib import Path
import json,hashlib,datetime,csv,io,os,subprocess
P=Path(__file__).resolve().parents[1];R=P.parent
def sha(f):return hashlib.sha256(Path(f).read_bytes()).hexdigest()
def writejson(f,o):Path(f).write_text(json.dumps(o,ensure_ascii=False,indent=2))
q=subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True)
(P/'logs/final_check_memory.log').write_text(q.stdout+q.stderr);assert q.returncode==0
for n in ['analyze','stats_v2','validate_v4','geo_sample_metadata']:assert (P/'logs'/f'{n}.exit').read_text().strip()=='0'
assert (P/'logs/finalize.exit').read_text().strip()=='1'
assert 'keys must be str, int, float, bool or None, not tuple' in (P/'logs/finalize.log').read_text()
roots=json.loads((P/'qc/ROOT_RECORD_UPDATES.json').read_text())
assert len(roots)==4
for x in roots:assert sha(x['path'])==x['sha256']
# Finalizer completed all writes; its only failure was a tuple-keyed console summary.
(P/'qc/FINALIZE_CONSOLE_RECOVERY.txt').write_text('finalize.py completed report, manifests, backups, and all four root-record writes. Only its final json.dumps console summary failed because pandas groupby created tuple dictionary keys. No analysis/result failure. This final check verifies recorded root hashes and delivery hashes rather than rerunning or duplicating record updates. Library date prefixes in GSE204809 are metadata-derived batch-risk proxies; no claim of experimentally proven batch effect is made.\n')
# Clarify metadata inference in the report before final package hashing.
report=P/'delivery/Phase292_计划执行与首轮结果.md'
original=report.read_bytes()
(P/'qc/REPORT_BEFORE_METADATA_CLARIFICATION.md').write_bytes(original)
text=original.decode()
needle='未声明供者配对，不能直接把跨批差异解释为独立TYMP抑制效应。'
replacement='批次风险依据官方文库标题中的日期标签识别，尚未证明实际批次效应或跨批可比性。未声明供者配对，不能直接把两套文库差异解释为独立TYMP抑制效应。'
assert needle in text
report.write_text(text.replace(needle,replacement))
# Keep manifest report SHA synchronized without duplicating entries.
mf=R/'RESULTS_MANIFEST.tsv';old=mf.read_bytes()
(P/'qc/RESULTS_MANIFEST_before_report_clarification.tsv').write_bytes(old)
reader=csv.DictReader(io.StringIO(old.decode()),delimiter='\t');fields=reader.fieldnames;rows=list(reader);count=0
for row in rows:
 if row['result_id']=='phase292_report':
  before=row['notes'].split(';SHA256=')[0];row['notes']=before+';SHA256='+sha(report);count+=1
assert count==1
out=io.StringIO();w=csv.DictWriter(out,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
assert mf.read_bytes()==old
mf.write_text(out.getvalue())
delivery=P/'delivery/DELIVERY_MANIFEST.tsv'
records=[]
for f in sorted((P/'delivery').iterdir()):
 if f.is_file() and f.name!='DELIVERY_MANIFEST.tsv':records.append(dict(file=f.name,bytes=f.stat().st_size,sha256=sha(f)))
with delivery.open('w') as f:
 w=csv.DictWriter(f,fieldnames=['file','bytes','sha256'],delimiter='\t');w.writeheader();w.writerows(records)
for x in records:assert (P/'delivery'/x['file']).stat().st_size==x['bytes'] and sha(P/'delivery'/x['file'])==x['sha256']
for name in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md']:assert (R/name).read_text().startswith('## Phase292')
v=json.loads((P/'qc/INDEPENDENT_VALIDATION_v2.json').read_text());assert v['checks']==364 and v['failed']==0
stamp=datetime.datetime.now().astimezone().isoformat()
writejson(P/'qc/ROOT_RECORD_UPDATES_FINAL.json',[dict(path=str(R/n),sha256=sha(R/n)) for n in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md','RESULTS_MANIFEST.tsv']])
writejson(P/'qc/COMPLETION.json',dict(status='COMPLETE_FIRST_EXECUTION_AND_FEASIBILITY',at=stamp,checks=364,analysis_exit=0,statistics_v2_exit=0,independent_validation_exit=0,terminal_artifact_integrity='PASS',console_only_failure_retained=True,delivery_files=len(records)+1,external_biological_validation='not_yet_performed',no_local_progress_records=True))
(P/'RUN_COMPLETE').write_text(stamp+'\n')
print(json.dumps(dict(status='PASS',checks=364,delivery_files=len(records)+1,root_records=4,report=str(report)),ensure_ascii=False))
