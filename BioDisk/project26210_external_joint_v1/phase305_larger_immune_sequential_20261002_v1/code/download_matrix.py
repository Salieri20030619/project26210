from pathlib import Path
import sys,json,requests,hashlib,time,datetime,traceback,zipfile,os
P=Path(__file__).resolve().parents[1];job=json.loads(Path(sys.argv[1]).read_text());name=job['name'];out=P/'downloads'/job['filename'];partial=out.with_suffix(out.suffix+'.partial');state=P/'qc'/(name+'_DOWNLOAD_STATE.json')
def status(s,**kw):state.write_text(json.dumps(dict(status=s,time=datetime.datetime.now().astimezone().isoformat(),pid=os.getpid(),**kw),indent=2))
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  while b:=f.read(8*1024*1024):h.update(b)
 return h.hexdigest()
try:
 assert (job['route']=='direct')==(job['size']>2000000000)
 for attempt in range(1,6):
  try:
   if out.exists():break
   sess=requests.Session();sess.trust_env=False
   if job['route']=='proxy':sess.proxies={'http':'http://127.0.0.1:7897','https':'http://127.0.0.1:7897'}
   have=partial.stat().st_size if partial.exists() else 0
   if have==job['size']:break
   headers={'Range':f'bytes={have}-'} if have else {}
   status('CONNECTING',attempt=attempt,bytes=have,expected=job['size'],route=job['route']);print('CONNECTING',name,attempt,have,job['route'],flush=True)
   with sess.get(job['url'],headers=headers,stream=True,timeout=(45,120)) as r:
    r.raise_for_status()
    if have and r.status_code!=206:raise RuntimeError('Resume request not honored; original partial retained')
    if r.status_code==206:assert r.headers.get('Content-Range','').startswith(f'bytes {have}-')
    start=time.monotonic();last=start
    with partial.open('ab' if have else 'wb') as f:
     for block in r.iter_content(1024*1024):
      if block:f.write(block);have+=len(block)
      if time.monotonic()-last>=60:
       status('DOWNLOADING',attempt=attempt,bytes=have,expected=job['size'],route=job['route']);print('BYTES',have,'OF',job['size'],flush=True);last=time.monotonic()
   assert have==job['size'],(have,job['size']);break
  except Exception as e:
   print('ATTEMPT_FAILED',attempt,repr(e),flush=True)
   if attempt==5:raise
   time.sleep(20*attempt)
 target=out if out.exists() else partial;status('VERIFYING',bytes=target.stat().st_size,expected=job['size'],route=job['route'])
 assert target.stat().st_size==job['size'];digest=sha(target);assert digest==job['sha256'],(digest,job['sha256'])
 with zipfile.ZipFile(target) as z:
  bad=z.testzip();assert bad is None,bad
  entries=[dict(path=f.filename,compressed=f.compress_size,uncompressed=f.file_size) for f in z.infolist()]
 (P/'metadata'/(name+'_ZIP_CONTENTS.json')).write_text(json.dumps(entries,indent=2))
 if target==partial:partial.rename(out)
 status('COMPLETE',bytes=out.stat().st_size,sha256=digest,route=job['route'],archive_entries=len(entries));(P/'logs'/(name+'_download.exit')).write_text('0');print('COMPLETE',name,flush=True)
except BaseException:
 traceback.print_exc();status('FAILED',partial_bytes=partial.stat().st_size if partial.exists() else 0,route=job['route']);(P/'logs'/(name+'_download.exit')).write_text('1');sys.exit(1)
