from pathlib import Path
import json,subprocess,concurrent.futures,datetime,hashlib,urllib.parse,os
import pandas as pd
P=Path(__file__).resolve().parents[1]
def now():return datetime.datetime.now().astimezone().isoformat()
def dump(f,o):Path(f).write_text(json.dumps(o,ensure_ascii=False,indent=2))
d=json.loads((P/'sources/E-MTAB-12767.json').read_text());assets={}
def walk(x):
 if isinstance(x,dict):
  if x.get('type')=='file' and 'path'in x:assets[x['path']]=x
  for v in x.values():walk(v)
 elif isinstance(x,list):
  for v in x:walk(v)
walk(d)
base=json.loads((P/'sources/ebi_info.json').read_text())['httpLink']+'/Files/'
manifest=[dict(name=k,bytes=int(v['size']),url=base+urllib.parse.quote(k)) for k,v in sorted(assets.items())]
dump(P/'DOWNLOAD_PLAN.json',dict(at=now(),source='E-MTAB-12767 public processed files',files=manifest,total_bytes=sum(x['bytes'] for x in manifest),samples='nine capture arrays representing six individuals; patient grouping from SDRF',local_tunnel_pid=3196,proxy='127.0.0.1:17899'))
s=pd.read_csv(P/'sources/sdrf_fire.txt',sep='\t')
cols=[c for c in s if c.startswith('Source Name') or c.startswith('Characteristics[') or c.startswith('Derived Array Data File')]
s[cols].drop_duplicates().to_csv(P/'tables/OFFICIAL_SAMPLE_METADATA.tsv',sep='\t',index=False)
print('PLAN',len(manifest),sum(x['bytes'] for x in manifest),flush=True)
print(s[[c for c in ['Source Name','Characteristics[individual]','Characteristics[cell type]','Characteristics[clinical information]','Characteristics[disease]'] if c in s]].drop_duplicates().to_string(index=False),flush=True)
def fetch(v):
 f=P/'downloads'/v['name'];part=f.with_suffix(f.suffix+'.part')
 if f.exists() and f.stat().st_size==v['bytes']:return dict(**v,status='PASS',sha256=hashlib.sha256(f.read_bytes()).hexdigest())
 with open(P/'logs'/(v['name']+'.curl.log'),'w') as log:
  r=subprocess.run(['curl','-q','-fL','--proxy','http://127.0.0.1:17899','--noproxy','','--connect-timeout','12','--max-time','180','--retry','1','--retry-delay','5','-C','-','-o',str(part),v['url']],stdout=log,stderr=subprocess.STDOUT,timeout=375)
 if r.returncode or not part.exists() or part.stat().st_size!=v['bytes']:return dict(**v,status='FAILED',exit=r.returncode,actual_bytes=part.stat().st_size if part.exists() else 0)
 os.rename(part,f);return dict(**v,status='PASS',sha256=hashlib.sha256(f.read_bytes()).hexdigest())
dump(P/'STATUS.json',dict(status='DOWNLOADING_EMTAB',at=now(),files=len(manifest)))
out=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
 futures=[ex.submit(fetch,v) for v in manifest]
 for fu in concurrent.futures.as_completed(futures):
  out.append(fu.result());dump(P/'DOWNLOAD_PROGRESS.json',dict(at=now(),completed=len(out),passed=sum(x['status']=='PASS' for x in out),total=len(manifest)));print(out[-1]['name'],out[-1]['status'],flush=True)
dump(P/'DOWNLOADED_MANIFEST.json',out)
ok=all(x['status']=='PASS' for x in out)
dump(P/'STATUS.json',dict(status='DOWNLOAD_VERIFIED_ANALYSIS_PENDING' if ok else 'DOWNLOAD_INCOMPLETE',at=now(),passed=sum(x['status']=='PASS' for x in out)))
(P/'DOWNLOAD_EXIT_STATUS').write_text('0\n' if ok else '1\n')
raise SystemExit(0 if ok else 1)

