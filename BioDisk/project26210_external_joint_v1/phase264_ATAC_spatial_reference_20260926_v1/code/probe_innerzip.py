from pathlib import Path
import requests,tarfile,struct,json
P=Path(__file__).resolve().parents[1]
u='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207493/suppl/GSE207493_scATAC-seq.tar.gz'
with requests.get(u,stream=True,timeout=(10,40)) as r:
 r.raise_for_status()
 with tarfile.open(fileobj=r.raw,mode='r|gz') as t:
  e=next(iter(t));f=t.extractfile(e);b=f.read(8192)
  (P/'metadata/GSE207493_innerzip_prefix.bin').write_bytes(b)
  pos=0
  for _ in range(10):
   if b[pos:pos+4]!=b'PK\x03\x04':break
   h=struct.unpack('<IHHHHHIIIHH',b[pos:pos+30]);fn=b[pos+30:pos+30+h[-2]].decode(errors='replace')
   print(json.dumps(dict(name=fn,compressed=h[7],uncompressed=h[8],flags=h[2],method=h[3])),flush=True)
   pos+=30+h[-2]+h[-1]+h[7]
   if pos+30>len(b):break

