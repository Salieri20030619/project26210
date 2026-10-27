from pathlib import Path
import io,json,hashlib
import pdfplumber
from pypdf import PdfReader,PdfWriter,Transformation
from reportlab.pdfgen import canvas
base=Path(r'C:\Users\31937\Downloads\26210');out=Path('C:/codex_test/output/26210/Figure3_cNMF');out.mkdir(parents=True,exist_ok=True)
sources=[base/'03_2.2_cNMF模块与TYMP_v1/02_cNMF_F1_TYMP_UMAP_v1.pdf',base/'03_2.2_cNMF模块与TYMP_v1/05_I_gene_contribution_reference_blue_v2.pdf',base/'03_2.2_EH重绘_v2/E_F1_ordered_smooth_gradient_v2.pdf',base/'03_2.2_EH重绘_v2/H_patient_program_expression_blocks_v2.pdf',base/'03_2.2_cNMF模块与TYMP_v1/06_J_program_study_composition_v1.pdf']
mm=72/25.4;writer=PdfWriter();page=writer.add_blank_page(190*mm,275*mm)
boxes=[(3,207,184,62),(2,79,50,120),(54,79,56,120),(112,79,76,120),(3,3,184,68)];manifest=[]
for i,(path,box) in enumerate(zip(sources,boxes)):
 p=PdfReader(path).pages[0];sw,sh=float(p.mediabox.width),float(p.mediabox.height)
 buf=io.BytesIO();c=canvas.Canvas(buf,pagesize=(sw,sh));c.setFillColorRGB(1,1,1)
 with pdfplumber.open(path) as d:
  for q in d.pages[0].extract_words():
   if q['text'] in ['I','E','H','J'] and q['x0']<sw*.08 and q['top']<sh*.1:
    c.rect(q['x0']-1,sh-q['bottom']-1,q['x1']-q['x0']+2,q['bottom']-q['top']+2,fill=1,stroke=0)
 c.save();p.merge_page(PdfReader(buf).pages[0])
 x,y,w,h=[v*mm for v in box];scale=min(w/sw,h/sh)
 page.merge_transformed_page(p,Transformation().scale(scale).translate(x+(w-sw*scale)/2,y+(h-sh*scale)/2))
 manifest.append(dict(panel=chr(65+i),source=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),scale=scale))
buf=io.BytesIO();c=canvas.Canvas(buf,pagesize=(190*mm,275*mm));c.setFont('Helvetica-Bold',7)
for x,y,t in [(3,272,'A  cNMF programs, F1 usage and TYMP'),(3,202,'B  Gene contributions'),(55,202,'C  F1-ordered RNA trends'),(113,202,'D  Patient-program profiles'),(3,74,'E  Program counts and study composition')]:c.drawString(x*mm,y*mm,t)
c.save();page.merge_page(PdfReader(buf).pages[0])
with (out/'Figure3_cNMF_candidate_composite_v1.pdf').open('wb') as f:writer.write(f)
(out/'composition_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
