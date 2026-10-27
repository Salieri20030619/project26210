import pathlib,datetime,json,hashlib,subprocess,concurrent.futures,csv,os,fcntl,sys,time,traceback,shutil
C=pathlib.Path(__file__).resolve().parent
P=C.parents[1]
ROOT=P.parent
PROXY='http://127.0.0.1:7897'
def now():return datetime.datetime.now().astimezone().isoformat()
def save(path,obj):
 t=path.with_name(path.name+'.tmp');t.write_text(json.dumps(obj,ensure_ascii=False,indent=2));t.replace(path)
def emit(obj):
 with (C/'events.jsonl').open('a') as f:f.write(json.dumps(obj,ensure_ascii=False)+'\n')
 print(json.dumps(obj,ensure_ascii=False),flush=True)
def md5(q):
 h=hashlib.md5()
 with q.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def record(text,status):
 note='\n\n## '+now()+' Phase216 proxy continuation and full integrity verification\n'+text+'\nEvidence: '+str(C)+'\n'
 for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
  q=ROOT/name;q.write_text(note+q.read_text())
 with (ROOT/'RESULTS_MANIFEST.tsv').open('a') as f:
  f.write('\n'+'\t'.join([C.name+'_'+status,'phase216_download_integrity','json',str(C/'QA.json'),now(),status,'120 FASTQ expected; size, official MD5 and gzip CRC/EOF','No biological analysis; treatment identity remains unconfirmed'])+'\n')

def verify(t):
 q=P/'fastq'/t['filename'];s=q.stat();r=dict(t,started_at=now(),observed_bytes=s.st_size)
 r['observed_md5']=md5(q)
 cp=subprocess.run(['gzip','-t','--',str(q)],capture_output=True,text=True)
 r.update(gzip_exit_code=cp.returncode,gzip_error=cp.stderr.strip(),unchanged_during_verification=(q.stat().st_size==s.st_size and q.stat().st_mtime_ns==s.st_mtime_ns),finished_at=now())
 r['status']='PASS' if r['observed_bytes']==t['bytes'] and r['observed_md5']==t['md5'] and cp.returncode==0 and r['unchanged_during_verification'] else 'FAIL'
 save(C/'files'/(t['filename']+'.json'),r)
 return r

def main():
 lock=(P/'logs/download.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 rows=list(csv.DictReader((P/'metadata/PUBLIC_FASTQ_MANIFEST.tsv').open(),delimiter='\t'))
 tasks=[]
 for row in rows:
  triples=list(zip(row['fastq_urls'].split(';'),row['fastq_bytes'].split(';'),row['fastq_md5'].split(';')))
  assert len(triples)==2
  for u,s,h in triples:tasks.append(dict(run=row['run_accession'],sample_alias=row['sample_alias'],group=row['group'],url=u,bytes=int(s),md5=h,filename=u.rsplit('/',1)[1]))
 assert len(rows)==60 and len(tasks)==120 and len({t['filename'] for t in tasks})==120
 (C/'files').mkdir()
 (C/'backups').mkdir()
 for name in ['STATUS.json','DOWNLOAD_SUMMARY.json','VERIFIED_FILES.tsv','EXIT_STATUS','DOWNLOAD_TRANSPORT_POLICY.json','DIRECT_DOWNLOAD_POLICY.json']:
  q=P/name
  if q.exists():shutil.copy2(q,C/'backups'/name)
 for name in ['PAUSED_BY_SCHEDULE','PAUSED_BY_USER']:
  q=P/name
  if q.exists():q.rename(C/'backups'/name)
 policy=dict(status='PROXY_RESUME_AND_FULL_INTEGRITY_AUTHORIZED',proxy=PROXY,proxy_only=True,direct_fallback=False,scheduled_stop_at=None,control_dir=str(C),authorized_request='Configure proxy, continue downloading, then verify all data integrity',retain_complete_and_partial_files=True)
 save(P/'DOWNLOAD_TRANSPORT_POLICY.json',policy)
 save(P/'DIRECT_DOWNLOAD_POLICY.json',dict(status='SUPERSEDED_BY_EXPLICIT_PROXY_AUTHORIZATION',active_policy=str(P/'DOWNLOAD_TRANSPORT_POLICY.json'),updated_at=now()))
 save(C/'CONFIG.json',dict(policy,verification_workers=4,validation=['official_size','fresh_official_MD5','gzip_CRC_and_EOF'],expected_files=120,expected_runs=60,manifest_sha256=hashlib.sha256((P/'metadata/PUBLIC_FASTQ_MANIFEST.tsv').read_bytes()).hexdigest()))
 state=dict(status='PROXY_DOWNLOADING',pid=os.getpid(),control_dir=str(C),started_at=now(),expected_files=120)
 save(P/'STATUS.json',state);save(C/'STATUS.json',state)
 record('User explicitly authorized proxy continuation and verification of all data. Prior one-time 07:00 pause superseded; no new deadline. Preserve completed files and resume partial data. Fresh size/MD5 and gzip integrity verification of all 120 FASTQs will follow. Started PID '+str(os.getpid())+'.','RUNNING')
 for t in tasks:
  dest=P/'fastq'/t['filename']
  if dest.exists():
   assert dest.stat().st_size==t['bytes'],'Existing file size mismatch; original preserved: '+str(dest)
   continue
  part=dest.with_name(dest.name+'.part')
  for attempt in range(1,13):
   emit(dict(event='DOWNLOAD_ATTEMPT',filename=t['filename'],attempt=attempt,partial_bytes=part.stat().st_size if part.exists() else 0,time=now()))
   with (C/'curl.log').open('ab') as log:
    cp=subprocess.run(['curl','-q','--fail','--location','--silent','--show-error','--connect-timeout','60','--speed-limit','1024','--speed-time','300','--retry','5','--retry-delay','15','--retry-connrefused','--continue-at','-','--proxy',PROXY,'--noproxy','','--output',str(part),t['url']],stdout=log,stderr=log)
   if cp.returncode==0 and part.exists() and part.stat().st_size==t['bytes']:
    got=md5(part)
    if got==t['md5']:
     part.replace(dest);emit(dict(event='DOWNLOAD_VERIFIED',filename=t['filename'],bytes=t['bytes'],time=now()));break
    part.rename(part.with_name(part.name+'.badmd5.'+str(int(time.time()))))
   if attempt==12:raise RuntimeError('Proxy download failed after retries: '+t['filename'])
   time.sleep(15)
 state.update(status='VERIFYING_ALL_FILES',updated_at=now(),verified_files=0)
 save(P/'STATUS.json',state);save(C/'STATUS.json',state)
 emit(dict(event='FULL_VERIFICATION_STARTED',files=120,bytes=sum(t['bytes'] for t in tasks),time=now()))
 results=[]
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
  fs={pool.submit(verify,t):t for t in tasks}
  for future in concurrent.futures.as_completed(fs):
   result=future.result();results.append(result)
   state.update(verified_files=len(results),failed_files=sum(r['status']!='PASS' for r in results),updated_at=now())
   save(C/'STATUS.json',state);save(P/'STATUS.json',state)
   emit(dict(event='FILE_VERIFIED',filename=result['filename'],status=result['status'],completed=len(results),total=120,time=now()))
 results.sort(key=lambda r:r['filename'])
 fields=['run','sample_alias','group','filename','bytes','observed_bytes','md5','observed_md5','gzip_exit_code','unchanged_during_verification','status','started_at','finished_at','gzip_error']
 with (C/'ALL_FILES_INTEGRITY.tsv').open('w') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(results)
 bad=[r for r in results if r['status']!='PASS']
 groups={}
 for g in sorted({r['group'] for r in results}):
  rr=[r for r in results if r['group']==g];groups[g]=dict(files=len(rr),passed=sum(r['status']=='PASS' for r in rr),runs=len({r['run'] for r in rr}))
 qa=dict(status='COMPLETE_VERIFIED' if not bad else 'INTEGRITY_FAILED',finished_at=now(),expected_files=120,verified_files=len(results)-len(bad),expected_runs=60,complete_runs=sum(all(r['status']=='PASS' for r in results if r['run']==run) for run in {r['run'] for r in results}),expected_bytes=sum(t['bytes'] for t in tasks),verified_bytes=sum(r['bytes'] for r in results if r['status']=='PASS'),failed_files=len(bad),failures=bad,groups=groups,all_files_fresh_MD5=True,all_files_gzip_test=True,control_dir=str(C),biological_analysis_started=False,treatment_pairing_confirmed=False)
 save(C/'QA.json',qa)
 if bad:raise RuntimeError('Integrity verification failed for '+str(len(bad))+' files')
 for r in results:
  save(P/'verification'/(r['filename']+'.json'),dict(r,status='VERIFIED',completed_at=r['finished_at'],full_integrity_evidence=str(C/'files'/(r['filename']+'.json'))))
 with (P/'VERIFIED_FILES.tsv').open('w') as f:
  fields2=['run','sample_alias','group','filename','bytes','md5','observed_bytes','observed_md5','status','completed_at'];w=csv.DictWriter(f,fields2,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(dict(r,status='VERIFIED',completed_at=r['finished_at']) for r in results)
 save(P/'DOWNLOAD_SUMMARY.json',qa);save(P/'STATUS.json',qa);save(C/'STATUS.json',qa)
 (P/'EXIT_STATUS').write_text('0\n');(P/'COMPLETE_VERIFIED').write_text(now()+'\n'+str(C/'QA.json')+'\n')
 (C/'COMPLETE_VERIFIED').write_text(now()+'\n')
 record('All 120/120 FASTQs (60/60 paired runs; '+str(qa['verified_bytes'])+' compressed bytes) passed fresh official size and MD5 checks and gzip CRC/EOF integrity tests. Download complete. No biological analysis performed; candidate intervention/treatment identity remains unconfirmed. Full evidence in QA.json and ALL_FILES_INTEGRITY.tsv.','PASS')
 emit(qa)

code=1
try:main();code=0
except BaseException:
 error=traceback.format_exc();(C/'ERROR.txt').write_text(error);emit(dict(event='FAILED',error=error,time=now()));save(C/'STATUS.json',dict(status='FAILED',error=error,control_dir=str(C)))
finally:(C/'EXIT_STATUS').write_text(str(code)+'\n')
sys.exit(code)
