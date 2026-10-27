from pathlib import Path
import io,json,hashlib
import pdfplumber
from pypdf import PdfReader,PdfWriter,Transformation
from reportlab.pdfgen import canvas
base=Path(r'C:\Users\31937\Downloads\26210')
out=Path('C:/codex_test/output/26210/Figure2_CNV');out.mkdir(exist_ok=True,parents=True)
mm=72/25.4
sources=[base/'02_2.1恶性上皮亚克隆_v1/01_CNV_chromosome_heatmap_v2.pdf',base/'02_2.1恶性上皮亚克隆_v1/10_CNV_five_class_density_reference_H_horizontal_v3.pdf',base/'02_2.1_G_I参考图重绘_v3/G_CNV_reference_complete_v3.pdf',base/'02_2.1_I_K补图_v1/K2_malignant_vs_nonmalignant_radar_v2.pdf']
writer=PdfWriter();page=writer.add_blank_page(190*mm,275*mm)
boxes=[(15,171,160,99),(3,88,184,77),(3,5,108,77),(115,5,72,77)]
manifest=[]
for i,(path,box) in enumerate(zip(sources,boxes)):
 p=PdfReader(path).pages[0];sw,sh=float(p.mediabox.width),float(p.mediabox.height)
 # Hide only obsolete panel letters; retain all data, axes and legends.
 masks=[]
 with pdfplumber.open(path) as d:
  for q in d.pages[0].extract_words():
   if (i==1 and q['text'] in list('ABCDE')) or (i==2 and q['text']=='G'):
    masks.append((q['x0']-1,sh-q['bottom']-1,q['x1']-q['x0']+2,q['bottom']-q['top']+2))
 if masks:
  buf=io.BytesIO();c=canvas.Canvas(buf,pagesize=(sw,sh));c.setFillColorRGB(1,1,1)
  for rect in masks:c.rect(*rect,stroke=0,fill=1)
  c.save();p.merge_page(PdfReader(buf).pages[0])
 crop=[0,0,sw,sh]
 if i==1:crop=[0,15,sw,sh-25]
 if i==3:crop=[55,sh-428,sw,sh-45]
 p.cropbox.lower_left=crop[:2];p.cropbox.upper_right=crop[2:]
 x,y,w,h=[v*mm for v in box];cw,ch=crop[2]-crop[0],crop[3]-crop[1];scale=min(w/cw,h/ch)
 tx=x+(w-cw*scale)/2-crop[0]*scale;ty=y+(h-ch*scale)/2-crop[1]*scale
 page.merge_transformed_page(p,Transformation().scale(scale).translate(tx,ty))
 manifest.append(dict(panel=chr(65+i),source=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),scale=scale,crop=crop))
buf=io.BytesIO();c=canvas.Canvas(buf,pagesize=(190*mm,275*mm));c.setFont('Helvetica-Bold',8)
for x,y,t in [(3,272,'A  Chromosome-level CNV signal'),(3,168,'B  Five CNV classes in the epithelial UMAP'),(3,85,'C  Patient-specific CNV subclones'),(115,85,'D  Epithelial expression profiles')]:c.drawString(x*mm,y*mm,t)
c.setFont('Helvetica',6);c.drawString(119*mm,81*mm,'3 paired patients; exploratory')
c.save();page.merge_page(PdfReader(buf).pages[0])
writer.add_metadata({'/Title':'Figure 2: epithelial CNV classification and subclones'})
with (out/'Figure2_CNV_composite_v1.pdf').open('wb') as f:writer.write(f)
(out/'composition_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
print('COMPOSED')
