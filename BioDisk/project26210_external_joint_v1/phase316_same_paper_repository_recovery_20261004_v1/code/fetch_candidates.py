from pathlib import Path
import requests,json,subprocess,hashlib,concurrent.futures,io,zipfile
P=Path(__file__).resolve().parents[1]
api=json.loads((P/'resources/zenodo12634581.json').read_text());entries={x['key']:x for x in api['files']}
def download_example():
 e=entries['Example_Data_Human_CRC_ST_data.rds'];url=e['links']['self'];f=P/'downloads'/e['key']
 q=subprocess.run(['curl','--http1.1','-fL','--retry','1','--retry-all-errors','--connect-timeout','15','--max-time','180','--proxy','http://127.0.0.1:7897','-o',str(f),url],capture_output=True,text=True)
 (P/'logs/download_example.log').write_text(q.stderr);row={'file':f.name,'url':url,'exit':q.returncode,'expected_bytes':e['size'],'expected_md5':e['checksum']}
 if q.returncode==0:row.update(bytes=f.stat().st_size,md5=hashlib.md5(f.read_bytes()).hexdigest(),sha256=hashlib.sha256(f.read_bytes()).hexdigest());row['valid']=row['bytes']==e['size'] and 'md5:'+row['md5']==e['checksum']
 (P/'qc/EXAMPLE_DOWNLOAD.json').write_text(json.dumps(row,indent=2));print('EXAMPLE',json.dumps(row),flush=True)
def get_samples():
 url='https://www.spatialtme.yelab.site/api/samples?organ_flask=&cancer_type_flask=&accessible_ID_flask=HRA000979&sample_id_flask=&platform='
 s=requests.Session();s.trust_env=False;s.proxies={'http':'http://127.0.0.1:7897','https':'http://127.0.0.1:7897'}
 try:
  r=s.get(url,timeout=30);r.raise_for_status();(P/'resources/spatialtme_samples.json').write_text(r.text);print('SAMPLES',r.text[:6000],flush=True)
 except Exception as e:print('SAMPLES_ERROR',str(e),flush=True)
def get_zip_index():
 e=entries['Figure34.zip'];url=e['links']['self'];N=e['size'];ses=requests.Session();ses.trust_env=False
 class Remote(io.RawIOBase):
  def __init__(self):self.pos=0;self.total=0
  def seekable(self):return True
  def readable(self):return True
  def tell(self):return self.pos
  def seek(self,n,whence=0):self.pos=n if whence==0 else (self.pos+n if whence==1 else N+n);return self.pos
  def read(self,n=-1):
   n=N-self.pos if n<0 else min(n,N-self.pos)
   if n<=0:return b''
   if n>8*1024*1024:raise RuntimeError('Metadata range exceeds bound')
   start=self.pos
   with ses.get(url,headers={'Range':f'bytes={start}-{start+n-1}'},timeout=(15,30),stream=True) as r:
    if r.status_code!=206:raise RuntimeError(f'Range not supported status={r.status_code}, refusing full download')
    data=r.raw.read(n+1)
    if len(data)!=n:raise RuntimeError(f'Range length {len(data)} != {n}')
   self.pos+=n;self.total+=n;return data
 f=Remote()
 try:
  with zipfile.ZipFile(f) as z:
   rows=[{'name':x.filename,'compressed':x.compress_size,'bytes':x.file_size,'offset':x.header_offset,'CRC':x.CRC,'compress_type':x.compress_type} for x in z.infolist()]
  (P/'resources/FIGURE34_ZIP_INDEX.json').write_text(json.dumps({'url':url,'archive_bytes':N,'range_bytes_read':f.total,'members':rows},indent=2));print('ZIP_INDEX',json.dumps(rows)[:16000],flush=True)
 except Exception as x:
  (P/'qc/ZIP_RANGE_ATTEMPT.json').write_text(json.dumps({'url':url,'error':str(x),'bytes_read':f.total},indent=2));print('ZIP_RANGE_ERROR',str(x),flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
 for f in [ex.submit(download_example),ex.submit(get_samples),ex.submit(get_zip_index)]:
  try:f.result()
  except Exception as e:print('ERROR',repr(e),flush=True)
