from pathlib import Path
import copy,io,json,subprocess,hashlib
from pypdf import PdfReader,PdfWriter,Transformation
from pypdf.generic import RectangleObject
from reportlab.pdfgen import canvas
from reportlab.lib.units import mm
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=r/'phase194_mainfigure_assembly_20260917_v1'
o=r/'phase191_discovery_order_report_20260917_v1/figures';n=r/'phase192_sections28_212_figures_20260917_v1/figures'
src={str(i):next(o.glob(f'Figure{i}_*.pdf')) for i in range(1,8)}
for i in range(1,9):src['n'+str(i)]=next(n.glob(f'{i:02d}_*.pdf'))
src['n6']=r/'phase193_section211_contours_20260917_v2/figures/06_C4_expression_contours_v2.pdf'
# Each crop uses normalized source coordinates: left, top, right, bottom.
plan=[
('Main_Figure_1_TYMP_programs','TYMP-associated malignant epithelial programs','Atlas and malignancy  /  Transcriptional programs  /  Regulatory and protein activity',[
('A','Cell atlas and major compartments','1',(.03,.028,.98,.24),(10,32,196,80),'2.1 / atlas'),
('B','Chromosome-level CNV evidence','2',(.075,.020,.920,.373),(216,32,194,124),'2.1'),
('C','cNMF programs, F1 usage and TYMP','3',(.027,.019,.97,.244),(10,138,196,75),'2.2'),
('D','Program genes and F1-ordered expression','3',(.043,.323,.574,.674),(10,240,155,202),'2.2 / 2.3'),
('E','Patient-balanced regulon activity','4',(.095,.25,.89,.631),(216,192,194,132),'2.4'),
('F','F1-associated inferred protein activity','5',(.063,.602,.985,.997),(178,350,232,109),'2.5')]),
('Main_Figure_2_MTA_SPP1_discovery','Metabolic context and discovery of an SPP1-enriched state','Metabolic programs  /  MTA biochemical context  /  Macrophage discovery  /  Matched tissue evidence',[
('A','Metabolic pathway activity in malignant cells','6',(.105,.025,.915,.187),(10,32,246,85),'2.6'),
('B','Biochemical route to MTA','n1',(.018,.10,.985,.95),(267,32,143,88),'2.8'),
('C','Malignant-cell dynamics: Monocle3','7',(.205,.334,.755,.647),(10,154,119,116),'2.7'),
('D','Malignant-cell dynamics: scTour','7',(.205,.666,.755,.985),(140,154,119,116),'2.7'),
('E','Discovery and localization of C4','n2',(.023,.12,.993,.955),(269,154,141,117),'2.9'),
('F','C4 multi-gene identity','n3',(.070,.105,.936,.941),(10,294,145,163),'2.9'),
('G','Measured MTA and tissue expression','n4',(.03,.10,.97,.997),(172,306,238,150),'2.10')]),
('Main_Figure_3_C4_function_dynamics','Functional identity and state ordering of C4 macrophages','Functional programs  /  Equal-expression contours  /  Original trajectories  /  Patient-balanced dynamics',[
('A','Macrophage functional programs','n5',(.015,.09,.977,.987),(10,32,187,138),'2.11'),
('B','Functional expression contours','n6',(.018,.065,.974,.988),(213,32,197,175),'2.11; revised contour v2'),
('C','Original partitioned myeloid trajectories','n7',(.025,.09,.978,.997),(10,243,187,151),'2.12'),
('D','Patient-balanced expression dynamics','n8',(.055,.076,.995,.951),(213,243,197,210),'2.12')])]
W,H=420*mm,470*mm
rows=[]
for name,title,subtitle,panels in plan:
 writer=PdfWriter();page=writer.add_blank_page(W,H)
 for letter,label,key,crop,rect,section in panels:
  f=src[key];sp=copy.deepcopy(PdfReader(f).pages[0])
  sw,sh=float(sp.mediabox.width),float(sp.mediabox.height)
  l,t,rr,b=crop;x0,y0,x1,y1=l*sw,(1-b)*sh,rr*sw,(1-t)*sh
  sp.cropbox=RectangleObject((x0,y0,x1,y1))
  xx,yt,ww,hh=[v*mm for v in rect];scale=min(ww/(x1-x0),hh/(y1-y0))
  pw,ph=(x1-x0)*scale,(y1-y0)*scale
  # Center horizontally, align panel body to top.
  tx=xx+(ww-pw)/2;ty=H-yt-ph
  page.merge_transformed_page(sp,Transformation().translate(-x0,-y0).scale(scale).translate(tx,ty))
  rows.append({'main_figure':name,'panel':letter,'label':label,'section':section,'source':str(f),'crop_normalized':crop,'slot_mm':rect,'effective_source_scale':scale})
 overlay=io.BytesIO();c=canvas.Canvas(overlay,pagesize=(W,H))
 c.setFont('Helvetica-Bold',19);c.drawString(10*mm,H-12*mm,title)
 c.setFillColorRGB(.38,.41,.43);c.setFont('Helvetica',10);c.drawString(10*mm,H-19*mm,subtitle)
 for letter,label,key,crop,rect,section in panels:
  xx,yt,ww,hh=rect
  c.setFillColorRGB(.14,.16,.18);c.setFont('Helvetica-Bold',15);c.drawString(xx*mm,H-(yt-4)*mm,letter)
  c.setFont('Helvetica-Bold',10.5);c.drawString((xx+7)*mm,H-(yt-4)*mm,label)
 c.setFont('Helvetica',7.5);c.setFillColorRGB(.45,.45,.45)
 c.drawString(10*mm,5*mm,'Standalone assembly | Original data and panel content retained | See panel provenance and figure legends')
 c.save();overlay.seek(0);page.merge_page(PdfReader(overlay).pages[0])
 f=d/'figures'/f'{name}.pdf'
 with f.open('wb') as out:writer.write(out)
 subprocess.run(['pdftoppm','-scale-to','2200','-singlefile','-jpeg',str(f),str(d/'qc'/name)],check=True)
 print('ASSEMBLED',name,flush=True)
(d/'tables/PANEL_PROVENANCE.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2))
