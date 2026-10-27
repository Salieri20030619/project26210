from pathlib import Path
import os,sys,json,time,hashlib,zipfile,requests,traceback,concurrent.futures,datetime,threading,shutil
P=Path(__file__).resolve().parents[1];j=json.loads((P/'resources/OBRADOVIC_DOWNLOAD_JOB.json').read_text());assert j['size']>2000000000
url=json.loads((P/'metadata/OBRADOVIC_HEAD_PROXY.json').read_text())['url'];out=P/'downloads'/j['filename'];parts=P/'downloads/OBRADOVIC_direct_parts';parts.mkdir(exist_ok=True)
state=P/'qc/OBRADOVIC_DOWNLOAD_STATE.json';N=8;step=(j['size']+N-1)//N;progress={};lock=threading.Lock();last=[0]
def update(status,force=False):
 with lock:
  if force or time.monotonic()-last[0]>60:
   n=sum(progress.values());state.write_text(json.dumps(dict(status=status,route='direct',pid=os.getpid(),bytes=n,expected=j['size'],connections=N,time=datetime.datetime.now().astimezone().isoformat()),indent=2));print(status,n,flush=True);last[0]=time.monotonic()
def part(i):
 start=i*step;end=min(j['size'],start+step)-1;size=end-start+1;f=parts/f'part{i:02d}'
 for attempt in range(6):
  try:
   have=f.stat().st_size if f.exists() else 0;assert have<=size;progress[i]=have
   if have==size:return f
   s=requests.Session();s.trust_env=False
   with s.get(url,headers={'Range':f'bytes={start+have}-{end}'},stream=True,timeout=(45,120)) as r:
    r.raise_for_status();assert r.status_code==206;assert r.headers.get('Content-Range','').startswith(f'bytes {start+have}-{end}/')
    with f.open('ab' if have else 'wb') as o:
     for b in r.iter_content(262144):
      if b:o.write(b);have+=len(b);progress[i]=have;update('DOWNLOADING')
   assert f.stat().st_size==size;return f
  except Exception as e:
   print('PART_RETRY',i,attempt,repr(e),flush=True)
   if attempt==5:raise
   time.sleep(20*(attempt+1))
try:
 update('CONNECTING',True)
 with concurrent.futures.ThreadPoolExecutor(N) as ex:files=list(ex.map(part,range(N)))
 update('ASSEMBLING',True);tmp=out.with_suffix('.assembling');h=hashlib.sha256()
 with tmp.open('wb') as o:
  for f in files:
   with f.open('rb') as z:
    while b:=z.read(8*1024*1024):o.write(b);h.update(b)
 assert tmp.stat().st_size==j['size'];assert h.hexdigest()==j['sha256'];update('VERIFYING_ZIP',True)
 with zipfile.ZipFile(tmp) as z:
  assert z.testzip() is None;entries=[dict(path=f.filename,compressed=f.compress_size,uncompressed=f.file_size) for f in z.infolist()]
 (P/'metadata/OBRADOVIC_ZIP_CONTENTS.json').write_text(json.dumps(entries,indent=2));assert not out.exists();tmp.rename(out)
 update('COMPLETE',True);(P/'logs/OBRADOVIC_download.exit').write_text('0')
 (P/'qc/OBRADOVIC_DOWNLOAD_VERIFIED.json').write_text(json.dumps(dict(bytes=j['size'],sha256=h.hexdigest(),route='direct',archive_entries=len(entries)),indent=2))
except BaseException:
 traceback.print_exc();update('FAILED',True);(P/'logs/OBRADOVIC_download.exit').write_text('1');sys.exit(1)
