from pathlib import Path
import json,re,hashlib
import pymupdf as fitz
from pypdf import PdfReader
from PIL import Image,ImageDraw
P=Path(__file__).resolve().parents[1]
pdf=P/'report/26210_manuscript_immune_integrated_v2.pdf'
d=fitz.open(pdf);b=json.loads((P/'report/report_materialized.json').read_text())
def clean(s):
 s=str(s).replace('−','-').replace('–','-').replace('—','-').replace('‑','-')
 supers=str.maketrans('⁻⁰¹²³⁴⁵⁶⁷⁸⁹','-0123456789')
 return re.sub(r'×10([⁻⁰¹²³⁴⁵⁶⁷⁸⁹]+)',lambda m:'e'+m[1].translate(supers),s)
def norm(s):return re.sub(r'\s+','',clean(s))
texts=[p.get_text(clip=fitz.Rect(0,0,p.rect.width,p.rect.height-32)) for p in d]
alltext=norm(''.join(texts));missing=[]
for i,x in enumerate(b):
 if x['type'] not in ['table','figure','supplementary_figure'] and x.get('text')!='补充图':
  if norm(x['text']) not in alltext:missing.append({'i':i,'text':x['text']})
 if x['type']=='table':
  for row in [x['headers']]+x['rows']:
   for v in row:
    if norm(v) and norm(v) not in alltext:missing.append({'table':x['caption'],'value':v})
(P/'qc/MISSING_TEXT.json').write_text(json.dumps(missing,ensure_ascii=False,indent=2))
assert not missing,missing[:5]
for s in json.loads((P/'inputs/SOURCE_MANIFEST.json').read_text()):
 assert hashlib.sha256(Path(s['path']).read_bytes()).hexdigest()==s['sha256']
headings=[];toc=d.get_toc()
for x in b:
 if x['type'] in ['h1','h2'] and x.get('text')!='补充图':
  hits=[z[2] for z in toc if norm(z[1])==norm(x['text'])]
  assert len(hits)==1,(x['text'],hits)
  headings.append({'text':x['text'],'page':hits[0]})
for c in ['1','2','3','4']:
 nums=[int(re.match(r'^\d+\.(\d+)',x['text'])[1]) for x in headings if re.match('^'+c+r'\.\d+',x['text'])]
 assert nums==list(range(1,len(nums)+1)),(c,nums)
captions=[x['caption'] for x in b if x['type']=='table' and x['caption'].startswith('表')]
assert [int(re.match(r'^表(\d+)',x)[1]) for x in captions]==list(range(1,14)),captions
refids=[int(re.match(r'^\[(\d+)\]',x.get('text',''))[1]) for x in b if re.match(r'^\[\d+\]',x.get('text',''))]
assert refids==list(range(1,36)),refids
citation_ids=set()
for x in b:
 if x.get('text','').startswith('['):continue
 for match in re.finditer(r'\[(\d+(?:[-–]\d+)?(?:,\d+(?:[-–]\d+)?)*)\]',x.get('text','')):
  for z in match[1].replace('–','-').split(','):
   citation_ids.update(range(int(z.split('-')[0]),int(z.split('-')[-1])+1))
assert citation_ids.issubset(set(refids)),citation_ids-set(refids)
assert not any('表免疫' in x.get('text','') or '[[ITABLE_' in x.get('caption','') for x in b)
placements=json.loads((P/'qc/PLACEMENTS.json').read_text());r=PdfReader(pdf);dest=r.named_destinations
expected={x['figure']:x['page']-1 for x in placements};links=0
for pi,pg in enumerate(r.pages):
 for rr in pg.get('/Annots',[]):
  a=rr.get_object()
  if '/Dest' in a:
   q=a['/Dest'];target=r.get_destination_page_number(dest[q]) if isinstance(q,str) else r.get_page_number(q[0].get_object())
   assert target in expected.values(),(pi,target)
   links+=1
for i,pg in enumerate(d):
 footer=pg.get_text(clip=fitz.Rect(0,pg.rect.height-32,pg.rect.width,pg.rect.height))
 assert str(i+1) in footer.split(),(i,footer)
 assert abs(pg.rect.width-595.2756)<.1 and abs(pg.rect.height-841.8898)<.1
# Check figures preserve every byte of original source and their native physical size.
base=P.parent/'phase273_clinical_suppfig_wording_20260928_v1'
for x in placements:
 f=Path(x['source']);assert hashlib.sha256(f.read_bytes()).digest()==hashlib.sha256((base/'figures'/f.name).read_bytes()).digest()
 assert x['scale']==1.
(P/'qc/PAGE_HEADINGS.json').write_text(json.dumps(headings,ensure_ascii=False,indent=2))
(P/'qc/FULL_TEXT.txt').write_text('\n\n'.join(f'PAGE {i+1}\n{x}' for i,x in enumerate(texts)))
(P/'qc/RENDER_VALIDATION.json').write_text(json.dumps({'status':'STRUCTURAL_PASS_VISUAL_PENDING','pages':len(d),'tables':41,'main_tables':13,'references':35,'native_figures':7,'links':links,'all_paragraphs_and_table_values_present':True,'continuous_sections':True,'sources_unchanged':True,'page_numbers_correct':True},indent=2))
out=P/'qc/render';out.mkdir(exist_ok=True)
for i,pg in enumerate(d):pg.get_pixmap(matrix=fitz.Matrix(110/72,110/72)).save(out/f'page-{i+1:02}.png')
# Two-page review sheets retain full page pixels; no downscaling.
for i in range(0,len(d),2):
 ims=[Image.open(out/f'page-{j+1:02}.png').convert('RGB') for j in range(i,min(i+2,len(d)))]
 sheet=Image.new('RGB',(sum(im.width for im in ims)+30,max(im.height for im in ims)+35),'#dddddd');xx=10
 for j,im in enumerate(ims):
  sheet.paste(im,(xx,25));ImageDraw.Draw(sheet).text((xx,5),f'Page {i+j+1}',fill='black');xx+=im.width+10
 sheet.save(P/f'qc/review_{i+1:02}_{min(i+2,len(d)):02}.png')
print('STRUCTURAL_PASS',len(d),'pages',links,'figure links; full preview generated')
print('MAIN_SECTION_PAGES',[(x['text'],x['page']) for x in headings if x['text'].startswith(('1.11','1.20','2.13','2.14','2.17','2.22','2.23','3.5','3.10','3.12','4.16'))])
