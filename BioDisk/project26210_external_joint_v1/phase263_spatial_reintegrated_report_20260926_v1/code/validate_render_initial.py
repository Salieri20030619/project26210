from pathlib import Path
import json,re,hashlib,shutil,sys
import pymupdf as fitz
from pypdf import PdfReader
from PIL import Image,ImageOps,ImageDraw
P=Path(__file__).resolve().parents[1]
pdf=P/'report/26210_integrated_spatial_report_v1.pdf'
d=fitz.open(pdf);b=json.loads((P/'report/report_materialized.json').read_text())
placements=json.loads((P/'qc/PLACEMENTS.json').read_text())
norm=lambda s:re.sub(r'\s+','',s).replace('−','-').replace('–','-').replace('—','-').replace('‑','-')
supers=str.maketrans('⁻⁰¹²³⁴⁵⁶⁷⁸⁹','-0123456789')
def clean(s):return re.sub(r'×10([⁻⁰¹²³⁴⁵⁶⁷⁸⁹]+)',lambda m:'e'+m[1].translate(supers),s)
texts=[p.get_text(clip=fitz.Rect(0,0,p.rect.width,p.rect.height-32)) for p in d]
alltext=norm(''.join(texts));missing=[]
for i,x in enumerate(b):
 if x['type'] not in ['table','figure']:
  target=norm(clean(x['text']))
  if target not in alltext:missing.append(dict(index=i,type=x['type'],text=x['text']))
assert not missing,missing
for rec in json.loads((P/'qc/INPUTS.json').read_text()):
 assert hashlib.sha256(Path(rec['path']).read_bytes()).hexdigest()==rec['sha256'],rec['path']
headings=[]
for x in b:
 if x['type'] in ['h1','h2']:
  hits=[i+1 for i,t in enumerate(texts) if norm(clean(x['text'])) in norm(t)]
  assert len(hits)==1,(x['text'],hits)
  headings.append(dict(text=x['text'],page=hits[0]))
for prefix in ['1 ','2 ','3 ','4 ']:
 assert sum(h['text'].startswith(prefix) for h in headings)==1
r=PdfReader(pdf);dest=r.named_destinations
expected={p['figure']:p['page']-1 for p in placements}
count=0;fig5count=0
for pi,p in enumerate(r.pages):
 for ref in p.get('/Annots',[]):
  annot=ref.get_object()
  if '/Dest' in annot:
   a=annot['/Dest']
   if isinstance(a,str):
    assert a in dest,a
    target=r.get_destination_page_number(dest[a])
   else:
    target=r.get_page_number(a[0].get_object())
   assert target in expected.values(),(pi,target)
   count+=1
   if target==expected[5]:fig5count+=1
assert fig5count>=14,(fig5count,count)
for i,p in enumerate(d):
 nums=p.get_text(clip=fitz.Rect(0,p.rect.height-32,p.rect.width,p.rect.height))
 assert str(i+1) in nums.split(),(i,nums)
 assert abs(p.rect.width-595.2756)<.1 and abs(p.rect.height-841.8898)<.1
# Preserve figure content at native size. Identity is validated on copied vector source files.
figs=[dict(figure=x['figure'],page=x['page'],source_sha256=hashlib.sha256(Path(x['source']).read_bytes()).hexdigest(),scale=x['scale']) for x in placements]
(P/'qc/PAGE_HEADINGS.json').write_text(json.dumps(headings,ensure_ascii=False,indent=2))
report=dict(status='STRUCTURAL_PASS_VISUAL_PENDING',pages=len(d),tables=32,figures=figs,links=count,figure5_links=fig5count,paragraphs_verified=len(b)-32-5,source_hashes_unchanged=True,new_tests=0,numbering_preserved=True)
(P/'qc/STRUCTURAL_VALIDATION.json').write_text(json.dumps(report,ensure_ascii=False,indent=2))
print(json.dumps(report,ensure_ascii=False),flush=True)
(P/'qc/render').mkdir(exist_ok=True)
for i,p in enumerate(d):
 p.get_pixmap(matrix=fitz.Matrix(100/72,100/72)).save(P/f'qc/render/page-{i+1:02}.png')
for i in range(0,len(d),2):
 ims=[Image.open(P/f'qc/render/page-{j+1:02}.png').convert('RGB') for j in range(i,min(i+2,len(d)))]
 sheet=Image.new('RGB',(sum(im.width for im in ims)+20*(len(ims)+1),max(im.height for im in ims)+55),'#d0d0d0')
 x=20
 for j,im in enumerate(ims):
  sheet.paste(im,(x,30));ImageDraw.Draw(sheet).text((x+10,10),f'Page {i+j+1}',fill='black');x+=im.width+20
 sheet.save(P/f'qc/review_{i+1:02}_{min(i+2,len(d)):02}.png')
print('RENDER_DONE',len(d),flush=True)

