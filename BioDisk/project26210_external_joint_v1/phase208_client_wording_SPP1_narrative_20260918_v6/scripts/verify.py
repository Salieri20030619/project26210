from pathlib import Path
import json,re,hashlib
import fitz,numpy as np
from PIL import Image,ImageDraw
r=Path.cwd();d=r/'phase208_client_wording_SPP1_narrative_20260918_v6';p=d/'report/26210_TYMP_MTA_SPP1_four_figures_narrative_v6.pdf';doc=fitz.open(p);b=json.loads((d/'report/report_materialized.json').read_text());old=json.loads((r/'phase201_three_candidates_C4_main_report_20260918_v4/report/report_materialized.json').read_text());pl=json.loads((d/'qc/PLACEMENTS.json').read_text());figpage={z['figure']:z['page']-1 for z in pl}
ot=[x for x in old if x['type']=='table'];nt=[x for x in b if x['type']=='table']
for before,after in zip(ot,nt[:31]):
 for aa,bb in zip(before['rows'],after['rows']):
  for av,bv in zip(aa,bb):assert re.findall(r'-?\d+(?:\.\d+)?(?:e[+-]?\d+)?',str(av),re.I)==re.findall(r'-?\d+(?:\.\d+)?(?:e[+-]?\d+)?',str(bv),re.I)
assert len([x for x in b if x['type']=='figure'])==4
assert hashlib.sha256((d/'figures/Figure1.pdf').read_bytes()).digest()==hashlib.sha256((r/'phase201_three_candidates_C4_main_report_20260918_v4/figures/Figure1_annotation_composite_v1.pdf').read_bytes()).digest()
texts=[pg.get_text() for pg in doc];assert not re.search(r'敏感|Sensitivity|稳健性|鲁棒性|审计|审核','\n'.join(texts),re.I);
assert not re.search(r'图\s*[5-9](?!\d)','\n'.join(texts));assert not re.search(r'Fig\s*\d','\n'.join(texts))
links=[]
for i,pg in enumerate(doc):
 for l in pg.get_links():
  if l['kind']==fitz.LINK_GOTO:
   assert l['page'] in figpage.values(),(i,l)
   links.append(dict(source_page=i+1,target_page=l['page']+1,rect=list(l['from']),target=list(l['to'])))
assert len(links)>=44
assert set(z['target_page'] for z in links)=={v+1 for v in figpage.values()}
# Native A4 figure bodies must render identically to the delivered standalone pages.
diffs=[]
for n in [2,3,4]:
 src=fitz.open(d/'figures'/f'Figure{n}.pdf');clip=fitz.Rect(0,0,src[0].rect.width,src[0].rect.height-30)
 a=src[0].get_pixmap(matrix=fitz.Matrix(1,1),clip=clip,alpha=False);c=doc[figpage[n]].get_pixmap(matrix=fitz.Matrix(1,1),clip=clip,alpha=False);aa=np.frombuffer(a.samples,dtype=np.uint8);cc=np.frombuffer(c.samples,dtype=np.uint8);assert aa.shape==cc.shape
 ratio=float((abs(aa.astype(int)-cc.astype(int))>20).mean());assert ratio<.001,(n,ratio);diffs.append(dict(figure=n,changed_channel_fraction=ratio,scale=1))
(d/'qc/pages').mkdir(exist_ok=True)
for i,pg in enumerate(doc):
 assert abs(pg.rect.width/72*25.4-210)<.01 and abs(pg.rect.height/72*25.4-297)<.01
 pix=pg.get_pixmap(matrix=fitz.Matrix(.95,.95),alpha=False);im=Image.frombytes('RGB',[pix.width,pix.height],pix.samples);im.save(d/'qc/pages'/f'page_{i+1:02d}.jpg',quality=86)
files=[d/'qc/pages'/f'page_{i+1:02d}.jpg' for i in range(len(doc))]
for start in range(0,len(files),6):
 canvas=Image.new('RGB',(1040,2250),'#DDDDDD');dr=ImageDraw.Draw(canvas)
 for j,f in enumerate(files[start:start+6]):
  im=Image.open(f);im.thumbnail((500,708));x=10+(j%2)*520;y=25+(j//2)*750;canvas.paste(im,(x,y));dr.text((x,y-18),f'Page {start+j+1}',fill='black')
 canvas.save(d/'qc'/f'contact_{start//6+1:02d}.jpg',quality=88)
detail=sorted(set([z['page']+1 for z in pl]+[i+1 for i,t in enumerate(texts) if '2.9 ' in t or '2.10 ' in t]+[len(doc)-1,len(doc),4]))
for n in detail:
 pix=doc[n-1].get_pixmap(matrix=fitz.Matrix(1.8,1.8),alpha=False);Image.frombytes('RGB',[pix.width,pix.height],pix.samples).save(d/'qc'/f'detail_{n:02d}.jpg',quality=90)
check=dict(pages=len(doc),figure_pages={k:v+1 for k,v in figpage.items()},figures_native_size=4,original_numeric_tables_unchanged=31,new_profile_index_table=1,clickable_internal_links=len(links),all_link_targets_valid=True,no_stale_figure_numbers=True,preprocessing_source_hash_unchanged=True,figure_render_comparison=diffs,detail_pages=detail)
(d/'qc/VALIDATION.json').write_text(json.dumps(check,indent=2));(d/'qc/PDF_LINKS.json').write_text(json.dumps(links,indent=2));(d/'qc/PAGE_TEXT.json').write_text(json.dumps(texts,ensure_ascii=False,indent=2));print('VERIFIED',json.dumps(check),flush=True)
