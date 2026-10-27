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


lock=(P/'logs/full_integrity.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
(C/'files').mkdir()
rows=list(csv.DictReader((P/'metadata/PUBLIC_FASTQ_MANIFEST.tsv').open(),delimiter='\t'))
tasks=[]
for row in rows:
 for u,s,h in zip(row['fastq_urls'].split(';'),row['fastq_bytes'].split(';'),row['fastq_md5'].split(';')):tasks.append(dict(run=row['run_accession'],sample_alias=row['sample_alias'],group=row['group'],filename=u.rsplit('/',1)[1],bytes=int(s),md5=h))
assert len(tasks)==120 and len({t['filename'] for t in tasks})==120
results=[]
present=[]
for t in tasks:
 if (P/'fastq'/t['filename']).exists():present.append(t)
 else:
  r=dict(t,status='MISSING',reason='Official ENA file URL currently resolves to an empty directory; HTTPS/HTTP/FTP checked; original partial retained',finished_at=now());results.append(r);save(C/'files'/(t['filename']+'.json'),r)
state=dict(status='VERIFYING_AVAILABLE_FILES',expected_files=120,available_files=len(present),missing_files=120-len(present),checked_files=0,control_dir=str(C),started_at=now())
save(C/'STATUS.json',state);save(P/'STATUS.json',state)
emit(state)
code=1
try:
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
  futures=[pool.submit(verify,t) for t in present]
  for f in concurrent.futures.as_completed(futures):
   r=f.result();results.append(r);state.update(checked_files=sum(x['status']!='MISSING' for x in results),passed_files=sum(x['status']=='PASS' for x in results),updated_at=now());save(C/'STATUS.json',state)
   emit(dict(event='FILE_VERIFIED',filename=r['filename'],status=r['status'],checked=state['checked_files'],available=len(present),time=now()))
 results.sort(key=lambda r:r['filename'])
 fields=['run','sample_alias','group','filename','bytes','observed_bytes','md5','observed_md5','gzip_exit_code','unchanged_during_verification','status','started_at','finished_at','gzip_error','reason']
 with (C/'ALL_FILES_INTEGRITY.tsv').open('w') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(results)
 good=[r for r in results if r['status']=='PASS'];bad=[r for r in results if r['status']!='PASS']
 qa=dict(status='COMPLETE_VERIFIED' if not bad else 'INCOMPLETE_SOURCE_UNAVAILABLE' if all(r['status']=='MISSING' for r in bad) else 'INTEGRITY_FAILED',finished_at=now(),expected_files=120,verified_files=len(good),failed_files=sum(r['status']=='FAIL' for r in results),missing_files=sum(r['status']=='MISSING' for r in results),complete_runs=sum(all(r['status']=='PASS' for r in results if r['run']==run) for run in {r['run'] for r in results}),verified_bytes=sum(r['bytes'] for r in good),expected_bytes=sum(t['bytes'] for t in tasks),failures=bad,all_available_files_fresh_MD5=True,all_available_files_gzip_test=True,control_dir=str(C),biological_analysis_started=False,treatment_pairing_confirmed=False)
 save(C/'QA.json',qa);save(C/'STATUS.json',qa);save(P/'STATUS.json',qa)
 for r in good:
  save(P/'verification'/(r['filename']+'.json'),dict(r,status='VERIFIED',completed_at=r['finished_at'],full_integrity_evidence=str(C/'files'/(r['filename']+'.json'))))
 save(P/'DOWNLOAD_SUMMARY.json',qa)
 note=str(len(good))+'/120 FASTQs passed fresh official size, MD5 and gzip CRC/EOF checks; '+str(qa['complete_runs'])+'/60 paired runs complete. Missing '+str(qa['missing_files'])+', corrupt '+str(qa['failed_files'])+'. Remaining official ENA URL currently resolves to an empty directory. Partial retained; not marked complete. Proxy source recovery evidence: '+str(C.with_name(C.name.removesuffix('_verification')))+'. No biological analysis.'
 record(note,qa['status'])
 code=0 if not bad else 1
 (P/'EXIT_STATUS').write_text(str(code)+'\n')
 if code==0:(P/'COMPLETE_VERIFIED').write_text(now()+'\n')
 (C/'VERIFICATION_FINISHED').write_text(now()+'\n')
 emit(qa)
except BaseException:
 (C/'ERROR.txt').write_text(traceback.format_exc());emit(dict(event='ERROR',error=traceback.format_exc()))
finally:(C/'EXIT_STATUS').write_text(str(code)+'\n')
sys.exit(code)
