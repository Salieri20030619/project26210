from pathlib import Path
import pymupdf as f,json,hashlib
from PIL import Image,ImageOps,ImageDraw
w=Path(__file__).resolve().parent;q=w/'qc';q.mkdir(exist_ok=True)
old=f.open(r'C:\Users\31937\Downloads\26210\Phase332_摘要屏障主线微调_v5\26210_manuscript_immune_figure_integrated_v5.pdf')
new=f.open(r'C:\Users\31937\Downloads\26210\Phase335_完整稿Fig2更新_v6\26210_manuscript_Fig2_updated_v6.pdf')
assert len(old)==len(new)==54
changed=[];same=[];invalid=[];links=[]
for n in range(len(new)):
 a=old[n].get_pixmap(matrix=f.Matrix(1,1));b=new[n].get_pixmap(matrix=f.Matrix(1,1))
 if a.samples!=b.samples:
  changed.append(n+1);new[n].get_pixmap(matrix=f.Matrix(1.5,1.5)).save(q/f'page_{n+1:02}.png')
 else:same.append(n+1)
 for l in new[n].get_links():
  if l['kind']==f.LINK_GOTO:
   links.append({'source':n+1,'target':l['page']+1,'to':list(l['to'])})
   if not 0<=l['page']<len(new):invalid.append(l)
assert not invalid
assert '图2G' not in ''.join(p.get_text() for p in new)
for st in range(0,len(changed),4):
 imgs=[]
 for n in changed[st:st+4]:
  im=Image.open(q/f'page_{n:02}.png').convert('RGB');im.thumbnail((595,842));tile=Image.new('RGB',(615,880),'#e8e8e8');tile.paste(im,((615-im.width)//2,30));ImageDraw.Draw(tile).text((10,8),f'Page {n}',fill='black');imgs.append(tile)
 sheet=Image.new('RGB',(1230,1760),'white')
 for k,im in enumerate(imgs):sheet.paste(im,((k%2)*615,(k//2)*880))
 sheet.save(q/f'contact_{st//4+1}.png')
(q/'PAGE_AND_LINK_CHECKS.json').write_text(json.dumps({'changed_pages':changed,'identical_pages':same,'links':links,'invalid_links':invalid,'pages':54},indent=2),encoding='utf-8')
print('Changed',changed,'Identical',len(same),'Links',len(links))
