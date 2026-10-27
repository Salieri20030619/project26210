from pathlib import Path
import io,json,hashlib
import pdfplumber
from pypdf import PdfReader,PdfWriter,Transformation
from reportlab.pdfgen import canvas
base=Path(r'C:\Users\31937\Downloads\26210');out=Path('C:/codex_test/output/26210/Figure4_pySCENIC');out.mkdir(parents=True,exist_ok=True)
sources=[base/'05_2.4_pySCENIC_v1/01_A_RSS_GENIE3_F1_F7_v1.pdf',base/'05_2.4_pySCENIC_v1/02_A_RSS_GRNBoost2_F1_F7_v1.pdf',base/'05_2.4_B_E_correction_v2/03_B_cNMF_patient_balanced_regulon_heatmap_annotations_v2.pdf',base/'05_2.4_B_E_correction_v2/06_E_gene_to_gene_regulon_network_v2.pdf']
mm=72/25.4;writer=PdfWriter();page=writer.add_blank_page(190*mm,275*mm)
boxes=[(3,215,90,52),(97,215,90,52),(3,98,184,109),(3,2,184,88)];manifest=[]
for i,(path,box) in enumerate(zip(sources,boxes)):
 p=PdfReader(path).pages[0];sw,sh=float(p.mediabox.width),float(p.mediabox.height)
 buf=io.BytesIO();c=canvas.Canvas(buf,pagesize=(sw,sh));c.setFillColorRGB(1,1,1)
 with pdfplumber.open(path) as d:
  for q in d.pages[0].extract_words():
   if q['text'] in ['B','E'] and q['x0']<sw*.08 and q['top']<sh*.1:c.rect(q['x0']-1,sh-q['bottom']-1,q['x1']-q['x0']+2,q['bottom']-q['top']+2,fill=1,stroke=0)
 c.save();p.merge_page(PdfReader(buf).pages[0])
 x,y,w,h=[v*mm for v in box];scale=min(w/sw,h/sh)
 page.merge_transformed_page(p,Transformation().scale(scale).translate(x+(w-sw*scale)/2,y+(h-sh*scale)/2))
 manifest.append(dict(panel=chr(65+i),source=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),scale=scale))
buf=io.BytesIO();c=canvas.Canvas(buf,pagesize=(190*mm,275*mm));c.setFont('Helvetica-Bold',8)
for x,y,t in [(3,271,'A  GENIE3 regulon specificity'),(97,271,'B  GRNBoost2 regulon specificity'),(3,211,'C  Patient-balanced regulon activity'),(3,94,'D  Candidate TF-to-target gene network')]:c.drawString(x*mm,y*mm,t)
c.save();page.merge_page(PdfReader(buf).pages[0])
with (out/'Figure4_pySCENIC_candidate_composite_v1.pdf').open('wb') as f:writer.write(f)
(out/'composition_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
