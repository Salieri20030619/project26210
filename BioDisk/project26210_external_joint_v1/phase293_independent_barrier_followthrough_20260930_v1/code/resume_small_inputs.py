from pathlib import Path
import subprocess,json,hashlib,concurrent.futures
P=Path(__file__).resolve().parents[1]
jobs=[
('inputs/xenium_rcc_core.h5ad','https://raw.githubusercontent.com/digvijayky/nicheverse/main/examples/data/xenium_rcc_core.h5ad'),
('inputs/GSE204809_HMVECS_Control_vs_RIG-I_Counts.xlsx','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE204nnn/GSE204809/suppl/GSE204809_HMVECS_Control_vs_RIG-I_Counts.xlsx'),
('resources/GSE294109_series_only.txt','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294109&targ=self&form=text&view=full'),
('resources/S-BIAD2470.json','https://www.ebi.ac.uk/biostudies/api/v1/studies/S-BIAD2470')]
def run(j):
 rel,u=j;f=P/rel;part=Path(str(f)+'.part')
 cmd=['curl','-L','--fail','--silent','--show-error','--connect-timeout','15','--max-time','180','--max-filesize','104857600']
 if part.exists():cmd+=['-C','-']
 cmd +=[u,'-o',str(part)]
 q=subprocess.run(cmd,capture_output=True,text=True,timeout=190)
 r=dict(path=str(f),url=u,exit=q.returncode,stderr=q.stderr)
 if not q.returncode:
  part.rename(f);r.update(bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest())
  if rel.endswith('.h5ad'):
   tree=json.loads((P/'resources/nicheverse_tree.json').read_text())
   e=[e for e in tree['tree'] if e['path']=='examples/data/xenium_rcc_core.h5ad'][0]
   blob=hashlib.sha1(b'blob '+str(f.stat().st_size).encode()+b'\0'+f.read_bytes()).hexdigest()
   assert blob==e['sha'] and f.stat().st_size==e['size'];r['official_git_blob_verified']=blob
 (P/'qc'/('resume_'+f.name+'.json')).write_text(json.dumps(r,indent=2))
 print(json.dumps(r),flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:list(ex.map(run,jobs))
