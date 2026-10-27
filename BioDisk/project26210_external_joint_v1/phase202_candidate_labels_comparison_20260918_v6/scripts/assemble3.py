from pathlib import Path
import io,copy,json,subprocess
import pandas as pd
from pypdf import PdfReader,PdfWriter,Transformation
from pypdf.generic import RectangleObject
from reportlab.pdfgen import canvas
from PIL import Image,ImageChops
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1')
d=r/'phase202_candidate_labels_comparison_20260918_v6'
prev=r/'phase197_top20_network_dynamics_20260917_v4'
old=json.loads((r/'phase194_mainfigure_assembly_20260917_v1/tables/PANEL_PROVENANCE.json').read_text())
p0={str(z['main_figure'].split('_')[2])+z['panel']:(Path(z['source']),z['crop_normalized']) for z in old}
pan={
 'marker':(r/'phase192_sections28_212_figures_20260917_v1/figures/03_macrophage_marker_heatmap.pdf',[.025,.09,.965,.99]),
 'contours':(d/'panels/functional_candidate_locations.pdf',[0,.065,1,1]),
 'gsea':(r/'phase199_fig3_marker_GSEA_20260918_v5/figures/Standalone_C4_GSEA_four_programs_v5.pdf',[0,0,1,1]),
 'trajectory':(p0['3C'][0],[.025,.065,.978,.997]),
 'curves':(prev/'panels/expression_dynamics_curves.pdf',[0,0,1,1])
}
tight={}
for key,(path,box) in pan.items():
 p=copy.deepcopy(PdfReader(path).pages[0]);sw,sh=float(p.mediabox.width),float(p.mediabox.height)
 l,t,rr,b=box;x0,y0,x1,y1=l*sw,(1-b)*sh,rr*sw,(1-t)*sh
 p.cropbox=RectangleObject((x0,y0,x1,y1))
 w=PdfWriter();page=w.add_blank_page(x1-x0,y1-y0)
 page.merge_transformed_page(p,Transformation().translate(-x0,-y0))
 tmp=d/'panels'/f'raw_{key}.pdf'
 with tmp.open('wb') as f:w.write(f)
 subprocess.run(['pdftoppm','-r','100','-singlefile','-png',str(tmp),str(d/'qc'/f'bound_{key}')],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
 im=Image.open(d/'qc'/f'bound_{key}.png').convert('RGB')
 mask=ImageChops.difference(im,Image.new('RGB',im.size,'white')).convert('L').point(lambda x:255 if x>12 else 0)
 bb=mask.getbbox();assert bb
 pad=4;bb=(max(0,bb[0]-pad),max(0,bb[1]-pad),min(im.width,bb[2]+pad),min(im.height,bb[3]+pad))
 w,h=x1-x0,y1-y0
 bx0,by0,bx1,by1=bb[0]/im.width*w,(1-bb[3]/im.height)*h,bb[2]/im.width*w,(1-bb[1]/im.height)*h
 page.cropbox=RectangleObject((bx0,by0,bx1,by1));tight[key]=(page,bx0,by0,bx1-bx0,by1-by0)
def height(k,w):return w*tight[k][4]/tight[k][3]
yc=24+max(height('marker',150),height('contours',165))+12
yd=yc+height('gsea',324)+12
panels=[
 ('A','Marker programs across macrophage states','marker',8,24,150),
 ('B','Functional expression contours','contours',167,24,165),
 ('C','Patient-ranked enrichment of C4 functional programs','gsea',8,yc,324),
 ('D','Original partitioned myeloid trajectories','trajectory',8,yd,155),
 ('E','Expression dynamics across myeloid states','curves',167,yd,165)]
mm=72/25.4;W=340*mm;H=(max(y+height(k,w) for _,_,k,x,y,w in panels)+9)*mm
wr=PdfWriter();page=wr.add_blank_page(W,H);records=[]
for letter,label,key,x,y,w in panels:
 p,bx,by,pw,ph=tight[key];scale=w*mm/pw
 page.merge_transformed_page(copy.deepcopy(p),Transformation().translate(-bx,-by).scale(scale).translate(x*mm,H-y*mm-ph*scale))
 records.append(dict(panel=letter,label=label,source=str(pan[key][0]),source_crop=pan[key][1],x_mm=x,y_mm=y,width_mm=w,height_mm=height(key,w)))
buf=io.BytesIO();ca=canvas.Canvas(buf,pagesize=(W,H))
ca.setFont('Helvetica-Bold',17);ca.drawString(8*mm,H-11*mm,'Macrophage identity, functional enrichment and state ordering')
for letter,label,key,x,y,w in panels:
 ca.setFont('Helvetica-Bold',12);ca.drawString(x*mm,H-(y-3)*mm,letter)
 ca.setFont('Helvetica-Bold',8.8);ca.drawString((x+6)*mm,H-(y-3)*mm,label)
ca.save();buf.seek(0);page.merge_page(PdfReader(buf).pages[0])
f=d/'figures/Main_Figure_3_candidate_labels_v6.pdf'
with f.open('wb') as out:wr.write(out)
(d/'tables/FIG3_PANEL_PROVENANCE.json').write_text(json.dumps(records,indent=2))
for path in (d/'figures').glob('*.pdf'):
 subprocess.run(['pdftoppm','-scale-to','2400','-singlefile','-jpeg',str(path),str(d/'qc'/path.stem)],check=True)
print('ASSEMBLY_READY',round(W/mm,1),round(H/mm,1))
