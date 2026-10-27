from pathlib import Path
import requests,json,datetime,re,hashlib,concurrent.futures
P=Path(__file__).resolve().parents[1]
urls={
 'geo_GSE204809':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE204809&targ=self&form=text&view=full',
 'xenium_ccrcc':'https://www.10xgenomics.com/datasets/xenium-protein-ffpe-human-renal-ccrcc',
 'zenodo_imc':'https://zenodo.org/api/records/18701672',
 'rcc_imc_readme':'https://raw.githubusercontent.com/BodenmillerGroup/RCC_publication/main/README.md'}
def get(kv):
 k,u=kv
 try:
  r=requests.get(u,timeout=25);body=r.content
  rec=dict(key=k,url=u,http=r.status_code,at=datetime.datetime.now().astimezone().isoformat(),bytes=len(body),sha256=hashlib.sha256(body).hexdigest())
  (P/'resources'/('PUBLIC_'+k+'.txt')).write_bytes(body)
  if k=='xenium_ccrcc':
   links=sorted(set(re.findall(r'https?[^\s<>"\\]+',r.text)))
   candidates=[x for x in links if any(y in x.lower() for y in ['cf.10x','gene_panel','json','csv','parquet','h5','rcc','xenium'])]
   (P/'resources/XENIUM_LINKS.json').write_text(json.dumps(candidates,indent=2))
   rec['candidate_links']=len(candidates)
  if k=='geo_GSE204809':
   rec['samples']=re.findall(r'GSM\d+',r.text)
   print(r.text[:7500],flush=True)
  return rec
 except Exception as e:return dict(key=k,url=u,error=repr(e))
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:out=list(ex.map(get,urls.items()))
(P/'resources/PUBLIC_METADATA_FETCH.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2),flush=True)
