from pathlib import Path
import json,subprocess,datetime,concurrent.futures,hashlib
P=Path(__file__).resolve().parents[1]
def now():return datetime.datetime.now().astimezone().isoformat()
urls=[('zenodo16833780.json','https://zenodo.org/api/records/16833780'),('zenodo7619249.json','https://zenodo.org/api/records/7619249'),('E-MTAB-12767.json','https://www.ebi.ac.uk/biostudies/api/v1/studies/E-MTAB-12767'),('mRCC_paper.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12496100/fullTextXML')]
def fetch(x):
 name,url=x
 r=subprocess.run(['curl','-q','-fL','--proxy','http://127.0.0.1:17899','--noproxy','','--connect-timeout','10','--max-time','35',url],capture_output=True,timeout=40)
 if r.returncode:return dict(name=name,url=url,status='FAILED',error=r.stderr.decode(errors='replace')[-300:],at=now())
 (P/'sources'/name).write_bytes(r.stdout)
 return dict(name=name,url=url,status='PASS',bytes=len(r.stdout),sha256=hashlib.sha256(r.stdout).hexdigest(),at=now())
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:res=list(ex.map(fetch,urls))
(P/'SOURCE_ACCESS.json').write_text(json.dumps(res,ensure_ascii=False,indent=2))
print(json.dumps(res),flush=True)
for f in (P/'sources').glob('zenodo*.json'):
 d=json.loads(f.read_text());print(f.name,d.get('metadata',{}).get('title'))
 for x in d.get('files',[]):print(x.get('key'),x.get('size'),x.get('checksum'),x.get('links',{}).get('self'))
(P/'STATUS.json').write_text(json.dumps(dict(status='METADATA_REVIEW_PENDING',at=now()),indent=2))

