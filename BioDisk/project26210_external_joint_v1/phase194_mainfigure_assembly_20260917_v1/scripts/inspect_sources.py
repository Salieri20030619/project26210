from pathlib import Path
import subprocess,json
from PIL import Image,ImageDraw
from pypdf import PdfReader
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1')
d=r/'phase194_mainfigure_assembly_20260917_v1'
for s in ['qc','scripts','figures','tables']:(d/s).mkdir(parents=True,exist_ok=True)
for f in sorted((r/'phase191_discovery_order_report_20260917_v1/figures').glob('*.pdf')):
 rr=PdfReader(f);print(f.name,[(float(p.mediabox.width),float(p.mediabox.height)) for p in rr.pages])
 subprocess.run(['pdftoppm','-scale-to','1600','-singlefile','-jpeg',str(f),str(d/'qc'/f.stem)],check=True)
blocks=json.loads((r/'phase191_discovery_order_report_20260917_v1/report/report_materialized.json').read_text())
for b in blocks:
 if b['type'] in ['h2','figure','caption']: print(json.dumps(b,ensure_ascii=False)[:1800])
