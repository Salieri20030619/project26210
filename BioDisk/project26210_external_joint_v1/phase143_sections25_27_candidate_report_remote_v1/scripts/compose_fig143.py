from pathlib import Path
import io,json,hashlib
from pypdf import PdfReader,PdfWriter,Transformation
from reportlab.pdfgen import canvas
b=Path(r'C:\Users\31937\Downloads\26210\06_2.5_density_reference_v2');out=Path('C:/codex_test/output/26210/Figure5_protein');out.mkdir(parents=True,exist_ok=True)
sources=[b/'01_MR_effect_activity_density_contours_v2.pdf',b/'02_F1_protein_activity_density_contours_v2.pdf'];mm=72/25.4
writer=PdfWriter();page=writer.add_blank_page(190*mm,255*mm)
def put(path,box,top=None,bottom=None):
 p=PdfReader(path).pages[0];sw,sh=float(p.mediabox.width),float(p.mediabox.height)
 lo=0 if bottom is None else sh-bottom;hi=sh if top is None else sh-top
 p.cropbox.lower_left=(0,lo);p.cropbox.upper_right=(sw,hi)
 x,y,w,h=[v*mm for v in box];scale=min(w/sw,h/(hi-lo))
 page.merge_transformed_page(p,Transformation().scale(scale).translate(x+(w-sw*scale)/2,y+(h-(hi-lo)*scale)/2-lo*scale))
put(sources[0],(20,112,150,138))
for j,(top,bottom) in enumerate([(68,281),(289,502),(510,723),(731,946)]):
 r,c=divmod(j,2);put(sources[1],(2+c*94,57-r*52,92,49),top,bottom)
buf=io.BytesIO();c=canvas.Canvas(buf,pagesize=(190*mm,255*mm));c.setFont('Helvetica-Bold',8)
c.drawString(3*mm,252*mm,'A  MR effects and inferred protein activity')
c.drawString(3*mm,109*mm,'B  F1 usage and inferred activity across 20 patients')
c.setFont('Helvetica',6);c.drawCentredString(95*mm,2*mm,'All four panels: x = F1 normalized usage; y = metaVIPER NES')
c.save();page.merge_page(PdfReader(buf).pages[0])
with (out/'Figure5_protein_density_composite_v1.pdf').open('wb') as f:writer.write(f)
(out/'composition_manifest.json').write_text(json.dumps([dict(source=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sources],ensure_ascii=False,indent=2),encoding='utf-8')
