from pathlib import Path
import json,re,hashlib,collections
import pymupdf as fitz
import numpy as np
from PIL import Image,ImageDraw
P=Path(__file__).resolve().parents[1];O=P.parent/'phase327_abstract_clarity_revision_20261006_v1'
n=fitz.open(P/'report/26210_manuscript_immune_figure_integrated_v4.pdf');o=fitz.open(O/'report/26210_manuscript_immune_integrated_v3.pdf')
def content(pg):return ''.join(w[4] for w in pg.get_text('words') if w[3]<pg.rect.height-30)
lookup=collections.defaultdict(list)
for i,pg in enumerate(o):lookup[content(pg)].append(i)
changed=[];same=[];links=[]
for i,pg in enumerate(n):
 candidates=lookup[content(pg)];match=None
 for j in candidates:
  a=pg.get_pixmap(matrix=fitz.Matrix(1,1),clip=fitz.Rect(0,0,pg.rect.width,pg.rect.height-30),alpha=False)
  b=o[j].get_pixmap(matrix=fitz.Matrix(1,1),clip=fitz.Rect(0,0,o[j].rect.width,o[j].rect.height-30),alpha=False)
  if a.samples==b.samples:match=j;break
 if match is not None:same.append({'new':i+1,'old':match+1})
 else:changed.append(i+1)
 for l in pg.get_links():
  if l['kind']==fitz.LINK_GOTO:
   assert 0<=l['page']<len(n),(i,l);links.append({'source':i+1,'target':l['page']+1})
headings=[]
for i,pg in enumerate(n):
 for line in pg.get_text().splitlines():
  if re.match(r'^[1-4]\.\d+\s',line):headings.append({'page':i+1,'heading':line})
(P/'qc/PAGE_REVIEW_MAP.json').write_text(json.dumps({'pages':len(n),'changed_pages':changed,'same_content_and_pixels_excluding_footer':same,'valid_internal_links':len(links),'links':links,'headings':headings},ensure_ascii=False,indent=2))
for i in changed:
 pg=n[i-1];pg.get_pixmap(matrix=fitz.Matrix(1.4,1.4),alpha=False).save(P/f'qc/page_{i:02d}.png')
for k in range(0,len(changed),2):
 chunk=changed[k:k+2];ims=[Image.open(P/f'qc/page_{i:02d}.png').convert('RGB') for i in chunk]
 out=Image.new('RGB',(sum(im.width for im in ims),max(im.height for im in ims)+25),'#d5dbe0');x=0;dr=ImageDraw.Draw(out)
 for i,im in zip(chunk,ims):out.paste(im,(x,25));dr.text((x+12,6),f'PAGE {i}',fill='black');x+=im.width
 out.save(P/f'qc/review_{k//2+1:02d}.png')
print('REVIEW',json.dumps({'pages':len(n),'changed':changed,'pixel_identical_count':len(same),'links':len(links)}))
