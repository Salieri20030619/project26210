from pathlib import Path
import io,copy,json,subprocess
from pypdf import PdfReader,PdfWriter,Transformation
from pypdf.generic import RectangleObject
from reportlab.pdfgen import canvas
from PIL import Image,ImageChops
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=r/'phase198_fig2_contrast_corr_revision_20260918_v5';old195=r/'phase195_mainfigures_compact_20260917_v2';prev=r/'phase194_mainfigure_assembly_20260917_v1'
old=json.loads((prev/'tables/PANEL_PROVENANCE.json').read_text())
pan={str(z['main_figure'].split('_')[2])+z['panel']:(Path(z['source']),z['crop_normalized']) for z in old}
pan['3C']=(pan['3C'][0],[.025,.065,.978,.997])
# Use exact user-specified network.
pan.update({'network':(old195/'panels/network_compact.pdf',[0,0,1,1])})
for name in ['protein_compact','program_compact','functional_gene_pairs','spp1_contribution','dynamic_compact']:
 pan[name]=(old195/'panels'/f'{name}.pdf',[0,0,1,1])
for key in ['top20_pathways','gene_function_network','C4_source_share','expression_dynamics_curves','tissue_compact']:
 pan[key]=(d/'panels'/f'{key}.pdf',[0,0,1,1])
pan['epi']=(next((r/'phase191_discovery_order_report_20260917_v1/figures').glob('Figure7_*.pdf')),[.27,.027,.68,.303])
pan['2E']=(d/'panels/C4_no_nested_letters.pdf',[.023,.12,.993,.955])
pan['corr']=(d/'panels/correlation_consistency.pdf',[0,0,1,1])
pan['tissue_compact']=(r/'phase197_top20_network_dynamics_20260917_v4/panels/tissue_compact.pdf',[0,0,1,1])
pan={k:pan[k] for k in ['top20_pathways','epi','2E','corr','tissue_compact']}
# Reconstruct tightly clipped vector panels. Rasterization is used only to detect blank margins.
tight={}
for key,(f,box) in pan.items():
 p=copy.deepcopy(PdfReader(f).pages[0]);sw,sh=float(p.mediabox.width),float(p.mediabox.height)
 l,t,rr,b=box;x0,y0,x1,y1=l*sw,(1-b)*sh,rr*sw,(1-t)*sh
 p.cropbox=RectangleObject((x0,y0,x1,y1))
 wr=PdfWriter();page=wr.add_blank_page(x1-x0,y1-y0);page.merge_transformed_page(p,Transformation().translate(-x0,-y0))
 tmp=d/'panels'/f'raw_{key}.pdf'
 with tmp.open('wb') as out:wr.write(out)
 subprocess.run(['pdftoppm','-r','100','-singlefile','-png',str(tmp),str(d/'qc'/f'bound_{key}')],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
 im=Image.open(d/'qc'/f'bound_{key}.png').convert('RGB')
 mask=ImageChops.difference(im,Image.new('RGB',im.size,'white')).convert('L').point(lambda x:255 if x>12 else 0)
 bb=mask.getbbox();assert bb
 pad=4;bb=(max(0,bb[0]-pad),max(0,bb[1]-pad),min(im.width,bb[2]+pad),min(im.height,bb[3]+pad))
 w,h=x1-x0,y1-y0
 bx0,by0,bx1,by1=bb[0]/im.width*w,(1-bb[3]/im.height)*h,bb[2]/im.width*w,(1-bb[1]/im.height)*h
 page.cropbox=RectangleObject((bx0,by0,bx1,by1))
 tight[key]=(page,bx0,by0,bx1-bx0,by1-by0)
def height(k,w):return w*tight[k][4]/tight[k][3]
figs=[]
yd=24+max(height('epi',59),height('corr',104))+12
ye=max(24+height('top20_pathways',150),yd+height('2E',165))+12
pp=[('A','Top 20 TYMP-associated metabolic pathways','top20_pathways',8,24,150),
('B','Epithelial trajectory','epi',167,24,59),
('C','F1-pathway associations','corr',230,24,104),
('D','Discovery and localization of C4','2E',167,yd,165),
('E','Measured MTA and tissue expression','tissue_compact',8,ye,324)]
figs.append(('Main_Figure_2_contrast_consistency_v5','Metabolic context and discovery of an SPP1-enriched state',pp))
records=[]
for name,title,panels in figs:
 mm=72/25.4;W=340*mm;H=(max(y+height(k,w) for _,_,k,x,y,w in panels)+9)*mm
 writer=PdfWriter();page=writer.add_blank_page(W,H)
 for letter,label,key,x,y,w in panels:
  p,bx,by,pw,ph=tight[key];scale=w*mm/pw
  page.merge_transformed_page(copy.deepcopy(p),Transformation().translate(-bx,-by).scale(scale).translate(x*mm,H-y*mm-ph*scale))
  records.append({'main':name,'panel':letter,'title':label,'source':str(pan[key][0]),'source_crop':pan[key][1],'width_mm':w,'height_mm':height(key,w),'x_mm':x,'y_mm':y})
 buf=io.BytesIO();ca=canvas.Canvas(buf,pagesize=(W,H))
 ca.setFont('Helvetica-Bold',17);ca.drawString(8*mm,H-11*mm,title)
 for letter,label,key,x,y,w in panels:
  ca.setFont('Helvetica-Bold',12);ca.drawString(x*mm,H-(y-3)*mm,letter)
  ca.setFont('Helvetica-Bold',8.8);ca.drawString((x+6)*mm,H-(y-3)*mm,label)
 ca.save();buf.seek(0);page.merge_page(PdfReader(buf).pages[0])
 f=d/'figures'/f'{name}.pdf'
 with f.open('wb') as out:writer.write(out)
 subprocess.run(['pdftoppm','-scale-to','2300','-singlefile','-jpeg',str(f),str(d/'qc'/name)],check=True)
 print(name,'340 x',round(H/mm,1),'mm',flush=True)
(d/'tables/PANEL_PROVENANCE.json').write_text(json.dumps(records,indent=2))
