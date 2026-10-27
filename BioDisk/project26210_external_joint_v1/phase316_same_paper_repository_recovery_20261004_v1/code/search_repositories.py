from pathlib import Path
import subprocess,json,hashlib,concurrent.futures
P=Path(__file__).resolve().parents[1]
urls={'zenodo12634581.json':'https://zenodo.org/api/records/12634581','zenodo18482631.json':'https://zenodo.org/api/records/18482631','zenodo14602110.json':'https://zenodo.org/api/records/14602110','spatialtme.html':'https://www.spatialtme.yelab.site/','github_author_main.json':'https://api.github.com/repos/youqiongye/CRC_scRNAseq/git/trees/main?recursive=1','github_uca.json':'https://api.github.com/repos/BIGHanLab/UCASpatial/git/trees/main?recursive=1'}
def get(item):
 name,url=item;f=P/'resources'/name
 q=subprocess.run(['curl','--http1.1','-fL','--retry','1','--retry-all-errors','--connect-timeout','12','--max-time','45','--proxy','http://127.0.0.1:7897','-o',str(f),url],capture_output=True,text=True)
 (P/'logs'/('fetch_'+name+'.log')).write_text(q.stderr)
 row={'name':name,'url':url,'exit':q.returncode,'bytes':f.stat().st_size if f.exists() else 0}
 if q.returncode==0:row['sha256']=hashlib.sha256(f.read_bytes()).hexdigest()
 print(json.dumps(row),flush=True);return row
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:rows=list(ex.map(get,urls.items()))
(P/'resources/REPOSITORY_FETCH_MANIFEST.json').write_text(json.dumps(rows,indent=2))
for f in (P/'resources').glob('zenodo*.json'):
 try:
  j=json.loads(f.read_text());print('RECORD',f.name,j.get('metadata',{}).get('title'));print('DESCRIPTION',j.get('metadata',{}).get('description','')[:2000]);print('FILES',[(a.get('key'),a.get('size'),a.get('checksum')) for a in j.get('files',[])])
 except Exception as e:print('PARSE',str(e))
