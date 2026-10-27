from pathlib import Path
import requests,hashlib,json,time,tarfile,concurrent.futures,datetime
P=Path(__file__).resolve().parents[1]
def get(acc):
 url=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/suppl/{acc}_RAW.tar'
 expected=int(next(x.split('\t')[3] for x in (P/'metadata'/f'{acc}_filelist.txt').read_text().splitlines() if x.startswith('Archive\t')))
 target=P/'downloads'/f'{acc}_RAW.tar';part=target.with_suffix('.tar.part')
 sess=requests.Session();sess.trust_env=False;sess.proxies={'http':'http://127.0.0.1:7897','https':'http://127.0.0.1:7897'}
 if not target.exists():
  with sess.get(url,timeout=(20,90),stream=True) as res:
   res.raise_for_status()
   with part.open('wb') as out:
    for chunk in res.iter_content(1024*1024):out.write(chunk)
  assert part.stat().st_size==expected,(acc,part.stat().st_size,expected)
  part.rename(target)
 assert target.stat().st_size==expected
 with tarfile.open(target) as tar:
  members=tar.getmembers();names=[x.name for x in members]
  assert all(x.isfile() and '/' not in x.name and '..' not in x.name for x in members)
  dest=P/'downloads'/acc;dest.mkdir(exist_ok=True)
  for member in members:
   f=dest/member.name
   if not f.exists():
    with tar.extractfile(member) as src,f.open('wb') as dst:
     while data:=src.read(1024*1024):dst.write(data)
   assert f.stat().st_size==member.size
 h=hashlib.sha256()
 with target.open('rb') as f:
  while chunk:=f.read(1024*1024):h.update(chunk)
 record=dict(accession=acc,url=url,bytes=expected,sha256=h.hexdigest(),members=names,complete_at=datetime.datetime.now().astimezone().isoformat())
 (P/'qc'/f'{acc}_DOWNLOAD.json').write_text(json.dumps(record,indent=2))
 print('COMPLETE',acc,expected,len(names),flush=True)
 return record
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
 result=list(ex.map(get,['GSE210038','GSE222703','GSE224630']))
(P/'qc/DOWNLOADS_COMPLETE.json').write_text(json.dumps(result,indent=2))
