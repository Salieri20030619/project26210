from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import requests,json,traceback,tarfile,io
P=Path(__file__).resolve().parents[1]
jobs={
'cptac_filelist':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE227nnn/GSE227898/suppl/filelist.txt',
'ccrcc_figshare':'https://api.figshare.com/v2/articles/19603750',
'ccrcc_code_tree':'https://api.github.com/repos/lessonskit/Single-cell-multi-omics-profiling-of-ccRCC/git/trees/main?recursive=1',
'long_code_tree':'https://api.github.com/repos/Dragonlongzhilin/RenalTumor/git/trees/main?recursive=1',
'heart_methods':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9364862/fullTextXML',
'ccrcc_methods':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9978887/fullTextXML'
}
def run(k,u):
 try:
  r=requests.get(u,timeout=(10,40));r.raise_for_status()
  (P/'metadata'/f'{k}.txt').write_bytes(r.content)
  return dict(key=k,url=u,status=r.status_code,bytes=len(r.content))
 except Exception as e:return dict(key=k,url=u,error=str(e))
with ThreadPoolExecutor(max_workers=4) as pool:
 for x in as_completed([pool.submit(run,k,u) for k,u in jobs.items()]):print(json.dumps(x.result()),flush=True)
u='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207493/suppl/GSE207493_scATAC-seq.tar.gz'
try:
 with requests.get(u,stream=True,timeout=(10,40)) as r:
  r.raise_for_status()
  with tarfile.open(fileobj=r.raw,mode='r|gz') as t:
   rows=[]
   for entry in t:
    rows.append(dict(name=entry.name,size=entry.size,type=str(entry.type)))
    if entry.isfile():break
  (P/'metadata/GSE207493_first_members.json').write_text(json.dumps(rows,indent=2))
  print('ATAC_FIRST_MEMBERS',json.dumps(rows),flush=True)
except Exception as e:print('ATAC_PROBE_ERROR',str(e),flush=True)

