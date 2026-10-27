from pathlib import Path
import requests,json,hashlib,concurrent.futures
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase308_SPP1_direct_perturbation_20261003_v1')
for d in ['code','metadata','downloads','resources','tables','qc','logs','delivery','inputs']:(P/d).mkdir(parents=True,exist_ok=True)
def session():
 s=requests.Session();s.trust_env=False;s.proxies={'http':'http://127.0.0.1:7897','https':'http://127.0.0.1:7897'};return s
urls={
'GSE274609.soft':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274609&targ=all&view=full&form=text',
'GSE316918.soft':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE316918&targ=all&view=full&form=text',
'PDAC_SPP1_article.html':'https://pmc.ncbi.nlm.nih.gov/articles/PMC12479095/',
'JITC_macSPP1_article.html':'https://pmc.ncbi.nlm.nih.gov/articles/PMC12829364/',
'LUAD_barrier_article.html':'https://link.springer.com/article/10.1186/s12967-026-08746-2',
'GSE274609_files.html':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274609/suppl/',
}
def fetch(kv):
 name,u=kv
 for attempt in range(3):
  try:
   r=session().get(u,timeout=45);f=P/'metadata'/name;f.write_bytes(r.content);record={'file':name,'url':u,'status':r.status_code,'bytes':len(r.content),'sha256':hashlib.sha256(r.content).hexdigest()};print(json.dumps(record),flush=True);return record
  except Exception as e:
   if attempt==2:return {'file':name,'url':u,'error_type':type(e).__name__}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:records=list(ex.map(fetch,urls.items()))
(P/'metadata/FETCH_MANIFEST.json').write_text(json.dumps(records,indent=2));print('METADATA_COMPLETE',flush=True)
