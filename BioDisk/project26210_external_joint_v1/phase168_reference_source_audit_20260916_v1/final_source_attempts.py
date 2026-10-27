from pathlib import Path
import json,subprocess,gzip,os,signal
O=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase168_reference_source_audit_20260916_v1')
target=str(O/'fetch_sources.py');stopped=[]
for p in Path('/proc').iterdir():
 if not p.name.isdigit():continue
 try:
  if target.encode() in (p/'cmdline').read_bytes().split(b'\0'):
   os.kill(int(p.name),signal.SIGTERM);stopped.append(int(p.name))
 except (FileNotFoundError,PermissionError,ProcessLookupError):pass
(O/'SLOW_FETCH_STOPPED.json').write_text(json.dumps(dict(pids=stopped,reason='Repository API obtained required text/code files; stop long pending Zenodo archive/publisher retrieval, preserve completed sources')))
urls={
 'openalex_original.json':'https://api.openalex.org/works/doi:10.1016/j.immuni.2022.02.001',
 'original_mmc1.pdf':'https://ars.els-cdn.com/content/image/1-s2.0-S1074761322000814-mmc1.pdf',
 'publisher_supplement_retry.pdf':'https://www.science.org/doi/suppl/10.1126/sciadv.aed0733/suppl_file/sciadv.aed0733_sm.pdf'}
def fetch(item):
 name,url=item
 r=subprocess.run(['curl','-L','-sS','--connect-timeout','8','--max-time','20','-w','%{http_code}',url,'-o',str(O/name)],capture_output=True,text=True,timeout=25)
 return dict(name=name,url=url,exit=r.returncode,http=r.stdout,error=r.stderr)
import concurrent.futures
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:logs=list(ex.map(fetch,urls.items()))
(O/'FINAL_SOURCE_ATTEMPTS.json').write_text(json.dumps(logs,indent=2));print(json.dumps(logs))
try:
 j=json.loads((O/'openalex_original.json').read_text())
 print('LOCATIONS',json.dumps(j.get('locations',[])))
except Exception:pass
D=O/'author_api'
for f in D.glob('*.R'):
 if 'Figure2' in f.name or 'Figure 3' in f.name:
  print('CODE',f.name,'\n','\n'.join(f.read_text().splitlines()[:48]))
