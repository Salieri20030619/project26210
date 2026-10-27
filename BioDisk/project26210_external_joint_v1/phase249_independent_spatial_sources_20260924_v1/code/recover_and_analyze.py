from pathlib import Path
import os,json,datetime,hashlib,subprocess,requests,time,ctypes,select,shutil,traceback
P=Path(__file__).resolve().parents[1];AN=P.parent/'phase250_independent_sixpatient_mainline_20260924_v1'
def now():return datetime.datetime.now().astimezone().isoformat()
def dump(f,o):Path(f).write_text(json.dumps(o,ensure_ascii=False,indent=2))
pid=int((P/'DOWNLOADER_PID').read_text())
proc=Path('/proc')/str(pid)/'cmdline'
if proc.exists() and b'download_emt.py' in proc.read_bytes():
 fd=ctypes.CDLL(None,use_errno=True).syscall(434,pid,0)
 if fd>=0:select.select([fd],[],[]);os.close(fd)
plan=json.loads((P/'DOWNLOAD_PLAN.json').read_text());out=[];events=[]
session=requests.Session();session.trust_env=False;session.proxies={'http':'http://127.0.0.1:17899','https':'http://127.0.0.1:17899'}
for item in plan['files']:
 f=P/'downloads'/item['name'];part=f.with_suffix(f.suffix+'.part')
 if f.exists() and f.stat().st_size==item['bytes']:
  out.append(dict(**item,status='PASS',sha256=hashlib.sha256(f.read_bytes()).hexdigest(),method='prior_completed_preserved'));continue
 success=False
 for attempt in range(3):
  url=item['url'] if attempt!=1 else 'https://www.ebi.ac.uk/biostudies/files/E-MTAB-12767/'+item['name']
  offset=part.stat().st_size if part.exists() else 0
  try:
   started=time.monotonic()
   with session.get(url,headers={'Range':'bytes='+str(offset)+'-'} if offset else {},stream=True,timeout=(35,60)) as response:
    response.raise_for_status()
    if offset and response.status_code==206:
     assert response.headers.get('Content-Range','').startswith('bytes '+str(offset)+'-');target=part;mode='ab'
    elif offset:
     target=f.with_suffix(f.suffix+'.recovery_full.part');mode='wb'
    else:target=part;mode='wb'
    with open(target,mode) as handle:
     for chunk in response.iter_content(1024*1024):
      if chunk:handle.write(chunk)
      if time.monotonic()-started>180:raise TimeoutError('file transfer deadline')
    assert target.stat().st_size==item['bytes'],(target.stat().st_size,item['bytes'])
   os.rename(target,f);success=True;events.append(dict(name=item['name'],attempt=attempt,at=now(),status='PASS',url=url));break
  except Exception as exc:
   events.append(dict(name=item['name'],attempt=attempt,at=now(),status='FAILED',error=str(exc),url=url))
   session.close();session=requests.Session();session.trust_env=False;session.proxies={'http':'http://127.0.0.1:17899','https':'http://127.0.0.1:17899'}
   if attempt<2:time.sleep(3)
  finally:dump(P/'RECOVERY_EVENTS.json',events)
 out.append(dict(**item,status='PASS' if success else 'FAILED',sha256=hashlib.sha256(f.read_bytes()).hexdigest() if success else None))
 dump(P/'RECOVERY_PROGRESS.json',dict(at=now(),completed=len(out),passed=sum(x['status']=='PASS' for x in out),total=len(plan['files'])))
 print(item['name'],out[-1]['status'],flush=True)
session.close()
if (P/'DOWNLOADED_MANIFEST.json').exists():shutil.copy2(P/'DOWNLOADED_MANIFEST.json',P/'DOWNLOADED_MANIFEST.initial_attempt.json')
dump(P/'DOWNLOADED_MANIFEST.json',out)
ok=all(x['status']=='PASS' for x in out)
(P/'DOWNLOAD_EXIT_STATUS').write_text('0\n' if ok else '1\n')
dump(P/'STATUS.json',dict(status='DOWNLOAD_VERIFIED_ANALYSIS_PENDING' if ok else 'DOWNLOAD_INCOMPLETE_AFTER_RECOVERY',at=now(),passed=sum(x['status']=='PASS' for x in out)))
if ok:
 with open(AN/'logs/analysis_recovery.log','w') as log:r=subprocess.run(['/mnt/BioDisk/wangmingyuan/project26210_envs/part4_advanced_py_v1/bin/python','-u',str(AN/'code/analyze_mainline.py')],stdout=log,stderr=subprocess.STDOUT)
 (AN/'RECOVERY_RUNNER_EXIT_STATUS').write_text(str(r.returncode)+'\n')
 dump(P/'RECOVERY_FINAL.json',dict(at=now(),status='PASS' if r.returncode==0 else 'ANALYSIS_FAILED',download_files=38,analysis_exit=r.returncode))
else:dump(P/'RECOVERY_FINAL.json',dict(at=now(),status='INCOMPLETE',passed=sum(x['status']=='PASS' for x in out)))

