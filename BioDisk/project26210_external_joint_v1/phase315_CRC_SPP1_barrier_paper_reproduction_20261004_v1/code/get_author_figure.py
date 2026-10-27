from pathlib import Path
import subprocess,json,hashlib
P=Path(__file__).resolve().parents[1]
url='https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41467-022-29366-6/MediaObjects/41467_2022_29366_Fig6_HTML.png'
f=P/'resources/author_Fig6.png'
x=subprocess.run(['curl','--http1.1','-fL','--retry','2','--retry-all-errors','--max-time','60','--proxy','http://127.0.0.1:7897','-o',str(f),url],capture_output=True,text=True)
(P/'logs/author_figure_download.log').write_text(x.stderr)
print(x.returncode,f.stat().st_size if f.exists() else 0)
if x.returncode==0:(P/'resources/AUTHOR_FIGURE_SOURCE.json').write_text(json.dumps({'url':url,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()},indent=2))
