from pathlib import Path
import subprocess,json,concurrent.futures
P=Path(__file__).resolve().parents[1]
u={'spatialtme_shiny_app.R':'https://raw.githubusercontent.com/Yelab2020/SpatialTME/main/SpatialTME_shiny/app.R','spatialtme_shiny_UI.R':'https://raw.githubusercontent.com/Yelab2020/SpatialTME/main/SpatialTME_shiny/UI.R','uca_example.md':'https://raw.githubusercontent.com/BIGHanLab/UCASpatial/main/docs/documentation/04_UCASpatial_Example.md','uca_data.md':'https://raw.githubusercontent.com/BIGHanLab/UCASpatial/main/docs/documentation/03_data.md'}
def get(x):
 n,url=x;f=P/'resources'/n;q=subprocess.run(['curl','--http1.1','-fL','--retry','1','--retry-all-errors','--max-time','45','--proxy','http://127.0.0.1:7897','-o',str(f),url],capture_output=True,text=True);(P/'logs'/('fetch_'+n+'.log')).write_text(q.stderr);return {'name':n,'url':url,'exit':q.returncode}
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:out=list(ex.map(get,u.items()))
(P/'resources/FETCH_MANIFEST3.json').write_text(json.dumps(out,indent=2));print(json.dumps(out),flush=True)
for n in ['spatialtme_shiny_app.R','spatialtme_shiny_UI.R','uca_example.md','uca_data.md']:
 f=P/'resources'/n
 if f.exists():
  print('FILE',n)
  for i,l in enumerate(f.read_text().splitlines()):
   if any(k in l.lower() for k in ['download','readrds','read.csv','p7','p6','example_data','load(']):print(i+1,l[:300])
