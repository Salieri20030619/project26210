from pathlib import Path
import json,subprocess,xml.etree.ElementTree as E
from PIL import Image,ImageDraw
from pypdf import PdfReader
d=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase190_mechanism_focused_report_20260917_v1')
p=d/'report/26210_TYMP_MTA_SPP1_macrophage_integrated_report_v2.pdf'
pages=PdfReader(p).pages;txt=[x.extract_text() or '' for x in pages]
(d/'qc/PAGE_TEXT.json').write_text(json.dumps(txt,ensure_ascii=False,indent=2))
assert len(pages)>=30 and all(t.strip() for t in txt)
for s in ['2.8','2.9','2.10','2.11','2.12','3.6','0.09755','1.44865']:
 assert any(s in t for t in txt),s
print('NEW_SECTION_PAGES',[i+1 for i,t in enumerate(txt) if any(k in t for k in ['2.8 多','2.9 公','2.10 无','2.11 C4','2.12 SPP1','3.2 MTA','3.6 TYMP'])],flush=True)
png=d/'qc/pages_final';png.mkdir(exist_ok=True)
subprocess.run(['pdftoppm','-scale-to','1200','-png',str(p),str(png/'page')],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
files=sorted(png.glob('page-*.png'));assert len(files)==len(pages)
contacts=d/'qc/contact_final';contacts.mkdir(exist_ok=True)
for start in range(0,len(files),4):
 canvas=Image.new('RGB',(1740,2460),'#cccccc')
 for k,fp in enumerate(files[start:start+4]):
  im=Image.open(fp).convert('RGB');im.thumbnail((850,1200))
  x=(k%2)*870+10;y=(k//2)*1230+25
  canvas.paste(im,(x,y));ImageDraw.Draw(canvas).text((x,y-20),fp.stem,fill='black')
 canvas.save(contacts/f'contact_{start//4+1:02}.jpg',quality=90)
print('RENDERED',len(files),flush=True)
