from pathlib import Path
import json,subprocess,xml.etree.ElementTree as ET
from PIL import Image,ImageDraw
w=Path(__file__).resolve().parents[1];q=w/'qc';rows=[]
for section in ['2_6','2_7']:
 thumbs=[]
 for p in sorted((w/('figures_'+section)).glob('*.pdf')):
  info=subprocess.run(['pdfinfo',str(p)],capture_output=True,text=True,check=True).stdout
  assert any(l.startswith('Pages:') and l.split(':')[1].strip()=='1' for l in info.splitlines())
  out=q/(section+'_'+p.stem+'_pdf')
  subprocess.run(['pdftoppm','-r','110','-singlefile','-png',str(p),str(out)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
  xml=subprocess.run(['pdftotext','-bbox',str(p),'-'],capture_output=True,text=True,check=True).stdout
  tree=ET.fromstring(xml);bad=[];n=0
  for page in tree.iter('{http://www.w3.org/1999/xhtml}page'):
   width=float(page.attrib['width']);height=float(page.attrib['height'])
   for word in page.iter('{http://www.w3.org/1999/xhtml}word'):
    a=word.attrib;n+=1
    if float(a['xMin'])<-.5 or float(a['yMin'])<-.5 or float(a['xMax'])>width+.5 or float(a['yMax'])>height+.5:bad.append(word.text)
  rows.append({'section':section,'name':p.name,'pages':1,'outside_page_text':bad,'word_count':n})
  im=Image.open(str(out)+'.png').convert('RGB');im.thumbnail((680,510));tile=Image.new('RGB',(720,550),'white');tile.paste(im,((720-im.width)//2,10));ImageDraw.Draw(tile).text((15,530),p.stem,fill='black');thumbs.append(tile)
 sheet=Image.new('RGB',(1440,550*((len(thumbs)+1)//2)),'#DFE5EB')
 for i,im in enumerate(thumbs):sheet.paste(im,((i%2)*720,(i//2)*550))
 sheet.save(q/('CONTACT_'+section+'.jpg'),quality=92)
(q/'PDF_VALIDATION.json').write_text(json.dumps(rows,indent=2));print(json.dumps(rows),flush=True)
