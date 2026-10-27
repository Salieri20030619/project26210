from pathlib import Path
import io,copy,json,subprocess
from pypdf import PdfReader,PdfWriter,Transformation
from pypdf.generic import RectangleObject
from reportlab.pdfgen import canvas
from PIL import Image,ImageChops
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=r/'phase196_selected_panels_20260917_v3';old195=r/'phase195_mainfigures_compact_20260917_v2';prev=r/'phase194_mainfigure_assembly_20260917_v1'
old=json.loads((prev/'tables/PANEL_PROVENANCE.json').read_text())
pan={str(z['main_figure'].split('_')[2])+z['panel']:(Path(z['source']),z['crop_normalized']) for z in old}
pan['3C']=(pan['3C'][0],[.025,.065,.978,.997])
# Use exact user-specified network.
pan.update({'network':(old195/'panels/network_compact.pdf',[0,0,1,1])})
for name in ['protein_compact','program_compact','functional_gene_pairs','spp1_contribution','dynamic_compact']:
 pan[name]=(old195/'panels'/f'{name}.pdf',[0,0,1,1])
pan['1D']=(d/'inputs/H_patient_program_expression_blocks_v2.pdf',[.05,0,1,1])
pan['2A']=(d/'panels/85_pathways_column.pdf',[0,0,1,1])
pan['epi']=(next((r/'phase191_discovery_order_report_20260917_v1/figures').glob('Figure7_*.pdf')),[.27,.027,.68,.303])
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
# All positions are top-aligned and derived from actual panel aspect ratios.
a=24;c=a+height('1A',157)+11
row2=max(c+height('1C',157),a+height('1B',154))+13
row3=row2+max(height('1D',142),height('1E',176))+13
pp=[('A','Cell atlas','1A',8,a,157),('B','Chromosome-level CNV','1B',178,a,154),('C','cNMF programs, F1 and TYMP','1C',8,c,157),('D','Patient-program expression structure','1D',8,row2,142),('E','Patient-balanced regulon activity','1E',156,row2,176),('F','Candidate TF-to-target network','network',8,row3,142),('G','F1-associated inferred protein activity','protein_compact',158,row3,174)]
figs.append(('Main_Figure_1_selected_v3','TYMP-associated malignant epithelial programs',pp))
# The 85-pathway panel occupies one full-height left column.
yr=24+max(height('epi',90),height('2B',95))+12
ye=yr+height('2E',191)+12
pp=[('A','All 85 pathways by TYMP status','2A',8,24,205),
('B','Malignant epithelial trajectory','epi',221,24,90),
('C','Biochemical route to MTA','2B',317,24,95),
('D','Discovery and localization of C4','2E',221,yr,191),
('E','Measured MTA and tissue expression','2G',221,ye,191)]
figs.append(('Main_Figure_2_column_v3','Metabolic context and discovery of an SPP1-enriched state',pp))
records=[]
for name,title,panels in figs:
 mm=72/25.4;W=(420 if name.startswith('Main_Figure_2') else 340)*mm;H=(max(y+height(k,w) for _,_,k,x,y,w in panels)+9)*mm
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
 print(name,round(W/mm,1),'x',round(H/mm,1),'mm',flush=True)
(d/'tables/PANEL_PROVENANCE.json').write_text(json.dumps(records,indent=2))
