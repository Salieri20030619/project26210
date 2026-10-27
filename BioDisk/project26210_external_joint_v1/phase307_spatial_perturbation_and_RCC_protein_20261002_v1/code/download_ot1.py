from pathlib import Path
import requests,json,hashlib,time,concurrent.futures,datetime,os
P=Path(__file__).resolve().parents[1]
def session():
 s=requests.Session();s.trust_env=False;s.proxies={'http':'http://127.0.0.1:7897','https':'http://127.0.0.1:7897'};return s
def download(f):
 j=json.loads(f.read_text());items=j['items'];assert len(items)==1;item=items[0];name=item['name'];size=int(item['size']);assert size<2_000_000_000,'Direct background route required';dest=P/'downloads'/name;tmp=dest.with_name(name+'.partial');s=session()
 if dest.exists():assert dest.stat().st_size==size
 else:
  r=s.get(item['download_url'],stream=True,timeout=(25,90));r.raise_for_status()
  with tmp.open('wb') as out:
   for block in r.iter_content(4*1024*1024):
    if block:out.write(block)
  assert tmp.stat().st_size==size;os.replace(tmp,dest)
 h=hashlib.sha256();md=hashlib.md5()
 with dest.open('rb') as inf:
  for block in iter(lambda:inf.read(8*1024*1024),b''):h.update(block);md.update(block)
 rec={'file':name,'bytes':size,'sha256':h.hexdigest(),'md5':md.hexdigest(),'official_hash':item.get('content_hash'),'official_hash_type':item.get('content_hash_name'),'metadata':str(f.relative_to(P)),'route':'proxy_each_file_below2GB','time':datetime.datetime.now().astimezone().isoformat()};(P/'qc'/('download_'+name+'.json')).write_text(json.dumps(rec,indent=2));print('VERIFIED_SIZE',name,size,flush=True);return rec
def retry_download(f):
 for attempt in range(1,6):
  try:return download(f)
  except Exception as e:
   print('RETRY',f.name,attempt,type(e).__name__,flush=True)
   if attempt==5:raise RuntimeError('Download failed: '+f.name+'; '+type(e).__name__) from None
   time.sleep(2)
files=sorted((P/'metadata').glob('share_Day*.json'),key=lambda f:int(json.loads(f.read_text())['items'][0]['size']))
try:
 with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:result=list(ex.map(retry_download,files))
 (P/'qc/OT1_DOWNLOAD_COMPLETE.json').write_text(json.dumps({'status':'COMPLETE_SIZE_SHA_RECORDED','files':len(result),'bytes':sum(r['bytes'] for r in result),'records':result},indent=2));(P/'logs/download_ot1.exit').write_text('0\n');print('OT1_DOWNLOAD_COMPLETE',len(result),flush=True)
except Exception as e:
 (P/'logs/download_ot1.exit').write_text('1\n');(P/'qc/OT1_DOWNLOAD_ERROR.json').write_text(json.dumps({'error':str(e)},indent=2));raise
