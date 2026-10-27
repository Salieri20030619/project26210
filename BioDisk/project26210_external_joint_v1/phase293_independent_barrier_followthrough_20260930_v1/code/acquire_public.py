from pathlib import Path
import subprocess,json,hashlib,datetime,concurrent.futures
P=Path(__file__).resolve().parents[1]
jobs=[
('GSE264586_metadata.txt','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE264586&targ=all&form=text&view=full','resources'),
('GSE294109_metadata.txt','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294109&targ=all&form=text&view=full','resources'),
('nicheverse_tree.json','https://api.github.com/repos/digvijayky/nicheverse/git/trees/main?recursive=1','resources'),
('xenium_rcc_core.h5ad','https://raw.githubusercontent.com/digvijayky/nicheverse/main/examples/data/xenium_rcc_core.h5ad','inputs'),
('10x_old_rcc.html','https://www.10xgenomics.com/datasets/xenium-protein-ffpe-human-renal-carcinoma','resources'),
('10x_current_rcc.html','https://www.10xgenomics.com/datasets/xenium-protein-ffpe-human-renal-ccrcc','resources')]
base='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE204nnn/GSE204809/suppl/'
for fn in ['GSE204809_HMVECS_Control_vs_RIG-I_Counts.xlsx','GSE204809_HUVECS_Control_vs_RIG_Counts.xlsx','GSE204809_HUVECS_Control_vs_TPI_Counts.xlsx']:jobs.append((fn,base+fn,'inputs'))
def run(j):
 fn,u,sub=j;dst=P/sub/fn
 cmd=['curl','-L','--fail','--silent','--show-error','--connect-timeout','15','--max-time','70','--max-filesize','104857600',u,'-o',str(dst)+'.part']
 q=subprocess.run(cmd,capture_output=True,text=True,timeout=80)
 r=dict(file=fn,url=u,exit=q.returncode,stderr=q.stderr,at=datetime.datetime.now().astimezone().isoformat())
 part=Path(str(dst)+'.part')
 if q.returncode==0 and part.exists():
  part.rename(dst);r.update(bytes=dst.stat().st_size,sha256=hashlib.sha256(dst.read_bytes()).hexdigest())
 (P/'qc'/('fetch_'+fn.replace('/','_')+'.json')).write_text(json.dumps(r,indent=2))
 print(json.dumps(r),flush=True)
 return r
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
 out=list(ex.map(run,jobs))
(P/'qc/PUBLIC_FETCH_MANIFEST.json').write_text(json.dumps(out,indent=2))
