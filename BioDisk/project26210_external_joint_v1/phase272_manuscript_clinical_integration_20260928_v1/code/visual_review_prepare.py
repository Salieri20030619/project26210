from pathlib import Path
import pymupdf as fitz,json,re
from PIL import Image,ImageDraw
p=Path('phase272_manuscript_clinical_integration_20260928_v1');d=fitz.open(p/'report/26210_manuscript_with_clinical_section_v1.pdf')
v=p/'qc/visual';v.mkdir(exist_ok=True)
for start in range(0,len(d),9):
 canvas=Image.new('RGB',(1320,1935),'#cccccc');dr=ImageDraw.Draw(canvas)
 for k in range(start,min(start+9,len(d))):
  im=Image.open(p/f'qc/render/page-{k+1:02}.png').convert('RGB');im.thumbnail((420,594));x=10+440*((k-start)%3);y=25+645*((k-start)//3)
  canvas.paste(im,(x,y));dr.text((x,y-18),f'Page {k+1}',fill='black')
 canvas.save(v/f'overview_{start+1:02}.png')
for pg in [1,5,6,19,24,25,26,27,39,40,41]:
 d[pg-1].get_pixmap(matrix=fitz.Matrix(125/72,125/72)).save(v/f'detail_{pg:02}.png')
b=json.loads((p/'report/report_materialized.json').read_text())
norm=lambda s:re.sub(r'\s+','',s).replace('−','-').replace('–','-').replace('—','-')
txt=norm(''.join(page.get_text() for page in d))
checks=[]
for i,x in enumerate(b):
 if x['type']=='table' and x.get('_new_clinical'):
  for row in [x['headers']]+x['rows']:
   for cell in row:assert norm(cell) in txt,(i,cell)
for x in b:
 if x['type']=='p':
  ms=re.findall(r'(?:第|见)\s*[12]\.\d+\s*节',x['text'])
  if ms:checks.append(ms)
(p/'qc/NEW_TABLE_TEXT_CHECK.json').write_text(json.dumps({'new_table_cells_present':True,'explicit_section_references':checks,'pages':len(d)},ensure_ascii=False,indent=2))
print('VISUAL_READY',len(d),len(list(v.iterdir())),checks)
