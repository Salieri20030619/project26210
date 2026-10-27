from pathlib import Path
import json,hashlib,subprocess,xml.etree.ElementTree as ET
from pypdf import PdfReader
from PIL import Image,ImageOps,ImageDraw
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=r/'phase201_three_candidates_C4_main_report_20260918_v4';old=r/'phase191_discovery_order_report_20260917_v1'
p=d/'report/26210_TYMP_MTA_SPP1_macrophage_C4_primary_report_v4.pdf'
reader=PdfReader(p);texts=[z.extract_text() for z in reader.pages]
assert len(reader.pages)==34
a=json.loads((old/'report/report_materialized.json').read_text());b=json.loads((d/'report/report_materialized.json').read_text())
assert [v for v in a if v['type']=='table']==[v for v in b if v['type']=='table']
assert sum(v['type']=='figure' for v in b)==8
checks={'pages':34,'original_tables_identical':31,'new_figures':1}
target=[]
for i,t in enumerate(texts,1):
 if any(s in t for s in ['2.9 三个','图8 三个','2.10 C4','3.3 SPP1','1.10 髓系']):target.append(i)
 for term in ['4,502','9.76%','59.63%','54.59%','97.64%']:
  checks['present_'+term]=any(term in tx for tx in texts)
assert all(v for k,v in checks.items() if k.startswith('present_'))
for f in (old/'figures').glob('*.pdf'):
 assert hashlib.sha256(f.read_bytes()).hexdigest()==hashlib.sha256((d/'figures'/f.name).read_bytes()).hexdigest()
checks['original_figure_hashes_equal']=7
pl=json.loads((d/'qc/PLACEMENTS.json').read_text())
for pos in pl:
 page=reader.pages[pos['page']-1];w=float(page.mediabox.width);h=float(page.mediabox.height)
 assert pos['x']>=0 and pos['y']>=30 and pos['x']+pos['width']<=w and pos['y']+pos['height']<=h
checks['vector_figures_inside_page']=8
f8=next(z for z in pl if z['figure']==8)
assert '图8 三个' in texts[f8['page']-1]
checks['figure8_caption_same_page']=f8['page']
# Full-page visual contact sheets and detailed pages around new material.
pages=sorted((d/'qc/pages').glob('page-*.jpg'))
assert len(pages)==34
for start in range(0,len(pages),6):
 canvas=Image.new('RGB',(1040,2265),'#DADADA');draw=ImageDraw.Draw(canvas)
 for j,f in enumerate(pages[start:start+6]):
  im=Image.open(f).convert('RGB');im.thumbnail((500,708))
  x=10+(j%2)*520;y=30+(j//2)*755
  canvas.paste(im,(x,y));draw.text((x,y-20),f'Page {start+j+1}',fill='black')
 canvas.save(d/'qc'/f'contact_{start//6+1:02d}.jpg',quality=90)
for n in sorted(set(target+[f8['page'],f8['page']-1,f8['page']+1,1,34])):
 if not 1<=n<=34:continue
 subprocess.run(['pdftoppm','-f',str(n),'-l',str(n),'-r','135','-singlefile','-jpeg',str(p),str(d/'qc'/f'detail_{n:02d}')],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
checks['detail_pages']=sorted(set(target+[f8['page'],f8['page']-1,f8['page']+1,1,34]))
(d/'qc/TEXT_AND_LAYOUT_CHECKS.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2))
(d/'qc/PAGE_TEXT.json').write_text(json.dumps(texts,ensure_ascii=False,indent=2))
print('QA_READY',json.dumps(checks,ensure_ascii=False),flush=True)
