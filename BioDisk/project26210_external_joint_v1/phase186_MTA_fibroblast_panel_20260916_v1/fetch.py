from pathlib import Path
import subprocess,concurrent.futures,gzip,json
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase186_MTA_fibroblast_panel_20260916_v1')
(r/'sources').mkdir(parents=True,exist_ok=True)
jobs=[('counts.txt.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE249nnn/GSE249082/suppl/GSE249082_AC-59_RawCounts.txt.gz'),('zenodo_metadata.json','https://zenodo.org/api/records/8130216')]
def run(j):
 n,u=j
 z=subprocess.run(['curl','-fLsS','--connect-timeout','10','--max-time','60',u,'-o',str(r/'sources'/n)],capture_output=True,text=True)
 out={'file':n,'url':u,'exit':z.returncode,'error':z.stderr}
 if z.returncode==0 and n.endswith('.gz'):
  try:
   b=gzip.open(r/'sources'/n,'rt');print('HEADER',b.readline()[:2000]);print('FIRST',b.readline()[:500]);b.close()
  except Exception as e:out['parse_error']=str(e)
 return out
out=list(concurrent.futures.ThreadPoolExecutor(max_workers=2).map(run,jobs))
(r/'RETRIEVAL.json').write_text(json.dumps(out,indent=2));print(json.dumps(out))
