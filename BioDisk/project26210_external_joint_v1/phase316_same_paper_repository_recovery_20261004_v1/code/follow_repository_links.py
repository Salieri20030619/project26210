from pathlib import Path
import subprocess,json,hashlib,concurrent.futures
P=Path(__file__).resolve().parents[1]
u={'spatialtme_app.js':'https://www.spatialtme.yelab.site/static/eRNA_QTL/js/app.js','spatialtme_tree.json':'https://api.github.com/repos/Yelab2020/SpatialTME/git/trees/main?recursive=1','uca_tutorial.html':'https://bighanlab.github.io/UCASpatial/','uca_repro_readme.md':'https://raw.githubusercontent.com/BIGHanLab/UCASpatial/main/Reproduction/Readme.md','author_issues.json':'https://api.github.com/repos/youqiongye/CRC_scRNAseq/issues?state=all&per_page=100'}
def get(x):
 n,url=x;f=P/'resources'/n
 z=subprocess.run(['curl','--http1.1','-fL','--retry','1','--retry-all-errors','--connect-timeout','12','--max-time','40','--proxy','http://127.0.0.1:7897','-o',str(f),url],capture_output=True,text=True)
 (P/'logs'/('fetch_'+n+'.log')).write_text(z.stderr)
 row={'file':n,'url':url,'exit':z.returncode,'bytes':f.stat().st_size if f.exists() else 0};print(json.dumps(row),flush=True);return row
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:rows=list(ex.map(get,u.items()))
(P/'resources/REPOSITORY_FETCH_MANIFEST2.json').write_text(json.dumps(rows,indent=2))
