from pathlib import Path
import subprocess,json,hashlib,time
P=Path(__file__).resolve().parents[1]
base='https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-022-29366-6/MediaObjects/'
files=['41467_2022_29366_MOESM9_ESM.xlsx','41467_2022_29366_MOESM7_ESM.xlsx','41467_2022_29366_MOESM2_ESM.pdf']
rows=[]
for name in files:
 url=base+name;out=P/'downloads'/name
 cmd=['curl','--fail','--location','--retry','1','--connect-timeout','20','--max-time','120','--proxy','http://127.0.0.1:7897','--output',str(out),url]
 q=subprocess.run(cmd,capture_output=True,text=True)
 (P/'logs'/('download_'+name+'.log')).write_text(q.stdout+q.stderr)
 row={'file':name,'url':url,'exit':q.returncode,'bytes':out.stat().st_size if out.exists() else 0}
 if q.returncode==0:row['sha256']=hashlib.sha256(out.read_bytes()).hexdigest()
 rows.append(row);print(json.dumps(row),flush=True)
(P/'resources/DOWNLOAD_MANIFEST_v1.json').write_text(json.dumps(rows,indent=2))
