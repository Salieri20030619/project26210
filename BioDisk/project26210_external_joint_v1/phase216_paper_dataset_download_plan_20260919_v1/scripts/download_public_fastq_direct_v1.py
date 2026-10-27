from pathlib import Path
import csv,json,hashlib,subprocess,concurrent.futures,datetime,sys,os,traceback,shutil,time,fcntl
P=Path(__file__).resolve().parents[1]
lock=open(P/'logs/download.lock','w'); fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
def now(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def save(path,obj):
 tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(obj,indent=2,ensure_ascii=False));tmp.replace(path)
def md5(path):
 h=hashlib.md5()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''): h.update(b)
 return h.hexdigest()
rows=list(csv.DictReader((P/'metadata/PUBLIC_FASTQ_MANIFEST.tsv').open(),delimiter='\t'))
rows.sort(key=lambda x:(0 if 'candidate' in x['group'] else 1,x['run_accession']))
tasks=[]
for r in rows:
 urls=r['fastq_urls'].split(';'); sizes=r['fastq_bytes'].split(';'); checks=r['fastq_md5'].split(';')
 assert len(urls)==len(sizes)==len(checks)==2
 for u,s,h in zip(urls,sizes,checks): tasks.append(dict(run=r['run_accession'],sample_alias=r['sample_alias'],group=r['group'],url=u,bytes=int(s),md5=h))
assert len(tasks)==120 and len({t['url'] for t in tasks})==120
save(P/'STATUS.json',dict(status='DOWNLOADING',started_at=now(),pid=os.getpid(),expected_files=len(tasks),expected_bytes=sum(t['bytes'] for t in tasks),treatment_pairing_confirmed=False,biological_analysis_started=False))
(P/'logs/PID').write_text(str(os.getpid()))
def one(t):
 name=t['url'].rsplit('/',1)[1]; dest=P/'fastq'/name; partial=dest.with_suffix(dest.suffix+'.part'); evidence=P/'verification'/(name+'.json')
 result=dict(t,filename=name,started_at=now())
 try:
  if dest.exists():
   assert dest.stat().st_size==t['bytes'] and md5(dest)==t['md5'],'existing completed file invalid'
  else:
   assert shutil.disk_usage(P).free>50*10**9,'less than 50GB reserve'
   success=False
   for attempt in range(1,4):
    with (P/'logs'/(name+'.curl.log')).open('ab') as log:
     cp=subprocess.run(['curl','--fail','--location','--silent','--show-error','--connect-timeout','60','--speed-limit','1024','--speed-time','300','--retry','5','--retry-delay','15','--retry-connrefused','--continue-at','-','--output',str(partial),t['url']],stdout=log,stderr=log)
    if cp.returncode==0 and partial.exists() and partial.stat().st_size==t['bytes']:
     got=md5(partial)
     if got==t['md5']: success=True;break
     partial.rename(partial.with_name(partial.name+'.badmd5.'+str(int(time.time()))))
    if attempt<3:time.sleep(15)
   assert success,'download/size/MD5 failed after retries'
   partial.replace(dest)
  result.update(status='VERIFIED',observed_bytes=dest.stat().st_size,observed_md5=t['md5'],completed_at=now())
 except Exception as e: result.update(status='FAILED',error=repr(e),completed_at=now())
 save(evidence,result)
 print(json.dumps(result,ensure_ascii=False),flush=True)
 return result
try:
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool: results=list(pool.map(one,tasks))
 ok=[r for r in results if r['status']=='VERIFIED'];bad=[r for r in results if r['status']!='VERIFIED']
 summary=dict(status='COMPLETE_VERIFIED' if not bad else 'INCOMPLETE_DOWNLOAD',completed_at=now(),expected_files=120,verified_files=len(ok),verified_bytes=sum(r['bytes'] for r in ok),failed_files=len(bad),failures=bad,treatment_pairing_confirmed=False,biological_analysis_started=False)
 save(P/'DOWNLOAD_SUMMARY.json',summary);save(P/'STATUS.json',summary)
 fields=['run','sample_alias','group','filename','bytes','md5','observed_bytes','observed_md5','status','completed_at']
 with (P/'VERIFIED_FILES.tsv').open('w') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(results)
 code=0 if not bad and len(ok)==120 else 1
 (P/'EXIT_STATUS').write_text(str(code)+'\n')
 if code==0:(P/'COMPLETE_VERIFIED').write_text(now()+'\n')
 sys.exit(code)
except Exception:
 (P/'logs/EXCEPTION.txt').write_text(traceback.format_exc());(P/'EXIT_STATUS').write_text('2\n');raise
