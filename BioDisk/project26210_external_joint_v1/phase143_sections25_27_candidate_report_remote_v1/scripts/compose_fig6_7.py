from pathlib import Path
import io,json,hashlib
from pypdf import PdfReader,PdfWriter,Transformation
from reportlab.pdfgen import canvas
b=Path(r'C:\Users\31937\Downloads\26210');out=Path('C:/codex_test/output/26210');mm=72/25.4
def compose(folder,name,sources,boxes,labels):
 d=out/folder;d.mkdir(parents=True,exist_ok=True);w=PdfWriter();p=w.add_blank_page(190*mm,275*mm);manifest=[]
 for path,box in zip(sources,boxes):
  q=PdfReader(path).pages[0];sw,sh=float(q.mediabox.width),float(q.mediabox.height);x,y,ww,hh=[v*mm for v in box];s=min(ww/sw,hh/sh)
  p.merge_transformed_page(q,Transformation().scale(s).translate(x+(ww-sw*s)/2,y+(hh-sh*s)/2));manifest.append(dict(source=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
 buf=io.BytesIO();c=canvas.Canvas(buf,pagesize=(190*mm,275*mm));c.setFont('Helvetica-Bold',8)
 for x,y,t in labels:c.drawString(x*mm,y*mm,t)
 c.save();p.merge_page(PdfReader(buf).pages[0])
 with (d/(name+'.pdf')).open('wb') as f:w.write(f)
 (d/'composition_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
compose('Figure6_metabolism','Figure6_metabolism_candidate_composite_v1',[b/'07_2.6_reference_layout_v2/01_four_pathway_UMAP_reference_v2.pdf',b/'07_2.6_reference_layout_v2/03_all85_pathways_TYMP_dotplot_reference_v2.pdf',b/'07_2.6_reference_layout_v2/06_scFEA_proxy_UMAP_reference_v2.pdf'],[(3,223,184,45),(3,53,184,161),(3,2,184,44)],[(3,272,'A  Four AUCell pathways'),(3,218,'B  All 85 pathways by TYMP status and cNMF program'),(3,48,'C  scFEA metabolic proxies')])
compose('Figure7_trajectory','Figure7_trajectory_candidate_composite_v1',[out/'Figure7_trajectory/14_Monocle3_original_UMAP_trajectory_v1.pdf',b/'08_2.7_gradient_heatmaps_v1/11_Monocle3_gradient_heatmap_v1.pdf',b/'08_2.7_gradient_heatmaps_v1/13_scTour_gradient_heatmap_v1.pdf'],[(3,189,184,79),(3,96,184,87),(3,2,184,88)],[(3,272,'A  Original Monocle3 trajectory'),(3,185,'B  Monocle3 expression dynamics'),(3,92,'C  scTour expression dynamics')])
