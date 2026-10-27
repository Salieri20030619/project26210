from pathlib import Path
import subprocess,json,hashlib,concurrent.futures
P=Path(__file__).resolve().parents[1]
b='https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-022-29366-6/MediaObjects/'
tasks=[('source.xlsx',b+'41467_2022_29366_MOESM9_ESM.xlsx',False),('processed_headers.txt',b+'41467_2022_29366_MOESM6_ESM.gz',True),('author_readme.md','https://raw.githubusercontent.com/youqiongye/CRC_scRNAseq/V1.0.0/README.md',False),('source_file_descriptions.pdf',b+'41467_2022_29366_MOESM3_ESM.pdf',False)]
def get(t):
 name,url,head=t;out=P/'downloads'/name
 cmd=['curl','--http1.1','--fail','--location','--retry','2','--retry-all-errors','--connect-timeout','15','--max-time','80','--proxy','http://127.0.0.1:7897','--output',str(out)]
 if head:cmd+=['--head']
 q=subprocess.run(cmd+[url],capture_output=True,text=True);(P/'logs'/('fetch2_'+name+'.log')).write_text(q.stdout+q.stderr)
 row={'file':name,'url':url,'exit':q.returncode,'bytes':out.stat().st_size if out.exists() else 0}
 if q.returncode==0:row['sha256']=hashlib.sha256(out.read_bytes()).hexdigest()
 print(json.dumps(row),flush=True);return row
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex: rows=list(ex.map(get,tasks))
(P/'resources/DOWNLOAD_MANIFEST_v2.json').write_text(json.dumps(rows,indent=2))
