from pathlib import Path
import subprocess,json,hashlib,concurrent.futures
P=Path(__file__).resolve().parents[1]
b='https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-022-29366-6/MediaObjects/'
u={'author_subtype_markers.xlsx':b+'41467_2022_29366_MOESM5_ESM.xlsx','author_function_score.R':'https://raw.githubusercontent.com/youqiongye/CRC_scRNAseq/V1.0.0/5.Function_score.R','author_stromal.R':'https://raw.githubusercontent.com/youqiongye/CRC_scRNAseq/V1.0.0/2.Subset_Stromal.R','author_myeloid.R':'https://raw.githubusercontent.com/youqiongye/CRC_scRNAseq/V1.0.0/2.Subset_myeloid.R'}
def get(x):
 n,url=x;f=P/'resources'/n;q=subprocess.run(['curl','--http1.1','-fL','--retry','1','--retry-all-errors','--max-time','60','--proxy','http://127.0.0.1:7897','-o',str(f),url],capture_output=True,text=True);(P/'logs'/('fetch_'+n+'.log')).write_text(q.stderr);d={'file':n,'url':url,'exit':q.returncode}
 if q.returncode==0:d.update(bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest())
 print(json.dumps(d),flush=True);return d
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:rows=list(ex.map(get,u.items()))
(P/'resources/DOWNLOAD_MANIFEST.json').write_text(json.dumps(rows,indent=2))
