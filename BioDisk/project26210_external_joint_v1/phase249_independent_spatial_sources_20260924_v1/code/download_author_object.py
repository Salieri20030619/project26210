from pathlib import Path
import json,subprocess,datetime,hashlib,os
P=Path(__file__).resolve().parents[1]
record=json.loads((P/'sources/zenodo7619249.json').read_text());v=record['files'][0]
f=P/'downloads'/v['key'];part=f.with_suffix(f.suffix+'.part')
with open(P/'logs/author_object_download.log','w') as log:
 r=subprocess.run(['curl','-q','-fL','--proxy','http://127.0.0.1:17899','--noproxy','','--connect-timeout','20','--max-time','600','--retry','1','-C','-','-o',str(part),v['links']['self']],stdout=log,stderr=subprocess.STDOUT,timeout=1250)
assert r.returncode==0 and part.stat().st_size==v['size']
h=hashlib.md5()
with open(part,'rb') as z:
 for b in iter(lambda:z.read(1024*1024),b''):h.update(b)
assert 'md5:'+h.hexdigest()==v['checksum']
os.rename(part,f)
(P/'AUTHOR_OBJECT_VALIDATION.json').write_text(json.dumps(dict(status='PASS',at=datetime.datetime.now().astimezone().isoformat(),path=str(f),bytes=v['size'],official_md5=v['checksum']),indent=2))
print('AUTHOR_OBJECT_MD5_PASS',flush=True)

