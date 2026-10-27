from pathlib import Path
import requests,json,io,zipfile,concurrent.futures,hashlib,xml.etree.ElementTree as E
P=Path(__file__).resolve().parents[1];R=P.parent;B=next(R.glob('phase316*'))
proxy={'http':'http://127.0.0.1:44581','https':'http://127.0.0.1:44581'}
def index():
 j=json.loads((B/'resources/zenodo18482631.json').read_text());e=j['files'][0];url=e['links']['self'];N=e['size'];s=requests.Session();s.trust_env=False
 class Remote(io.RawIOBase):
  def __init__(self):self.pos=0;self.total=0
  def seekable(self):return True
  def readable(self):return True
  def tell(self):return self.pos
  def seek(self,n,whence=0):self.pos=n if whence==0 else self.pos+n if whence==1 else N+n;return self.pos
  def read(self,n=-1):
   n=N-self.pos if n<0 else min(n,N-self.pos)
   if n<=0:return b''
   assert n<8000000
   q=s.get(url,headers={'Range':f'bytes={self.pos}-{self.pos+n-1}'},timeout=(15,45),proxies=proxy);q.raise_for_status();assert q.status_code==206;data=q.content;assert len(data)==n;self.pos+=n;self.total+=n;return data
 f=Remote()
 with zipfile.ZipFile(f) as z:rows=[{'name':x.filename,'bytes':x.file_size,'compressed':x.compress_size,'offset':x.header_offset,'CRC':x.CRC,'compress_type':x.compress_type} for x in z.infolist()]
 (P/'resources/ADDITIONAL_SOURCE_ZIP_INDEX.json').write_text(json.dumps({'url':url,'archive_bytes':N,'range_bytes':f.total,'members':rows},indent=2));print('INDEX',len(rows),'bytes',f.total,flush=True)
 for x in rows:
  if any(w in x['name'].lower() for w in ['crc','qi','spp1','fap','2022','marker']):print(x['name'],x['bytes'],flush=True)
def download(name,url):
 s=requests.Session();s.trust_env=False;q=s.get(url,timeout=(15,50),proxies=proxy);q.raise_for_status();(P/'resources'/name).write_bytes(q.content);print('FETCH',name,len(q.content),hashlib.sha256(q.content).hexdigest(),flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
 jobs=[ex.submit(index),ex.submit(download,'alternative_article.html','https://www.nature.com/articles/s41588-025-02380-2'),ex.submit(download,'author_github_current.json','https://api.github.com/repos/youqiongye/CRC_scRNAseq/git/trees/main?recursive=1')]
 for f in jobs:
  try:f.result()
  except Exception as e:print('FETCH_ERROR',repr(e),flush=True)
