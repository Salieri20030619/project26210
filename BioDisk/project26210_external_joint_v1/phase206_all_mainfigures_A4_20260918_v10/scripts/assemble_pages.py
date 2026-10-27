from pathlib import Path
import json,re,hashlib,math
import fitz,numpy as np
from PIL import Image,ImageChops
r=Path.cwd();d=r/'phase206_all_mainfigures_A4_20260918_v10';MM=72/25.4
old1=json.loads((r/'phase196_selected_panels_20260917_v3/tables/PANEL_PROVENANCE.json').read_text())[:7]
old2=json.loads((r/'phase205_fig2_concentric_20260918_v9/tables/FIG2_PANEL_PROVENANCE.json').read_text())
old3=json.loads((r/'phase202_candidate_labels_comparison_20260918_v6/tables/FIG3_PANEL_PROVENANCE.json').read_text())
# Each rectangle is an allocated panel box in final-page mm, not a scaled oversized figure.
layouts={1:[('A','Cell atlas',10,19,94,33),('B','Chromosome CNV',110,19,90,59),('C','cNMF programs',10,63,94,33),('D','Patient-program expression',10,107,94,90),('E','Regulon activity',110,90,90,101),('F','TF-target network',10,207,94,78),('G','Inferred protein activity',110,203,90,82)],2:[('A','Metabolic pathways',10,19,96,141),('C','F1-pathway associations',112,19,88,70),('B','Epithelial trajectory',112,104,39,49),('D','Macrophage candidates',155,104,45,49),('E','SPP1 expression and contribution',10,174,190,45),('F','Tissue MTA and expression',10,234,190,51)],3:[('A','Macrophage markers',10,19,90,115),('B','Expression contours',106,19,94,115),('C','C4 functional enrichment',10,146,190,72),('D','Myeloid trajectories',10,231,91,54),('E','Expression dynamics',107,231,93,54)]}
records=[];removed=[];keys=[]
remove_patterns=['258,120 cells','75 libraries','Other: erythrocyte','Patient-balanced regulon activity across','Malignant epithelial cells |','Within-partition normalized','Within-partition pseudotime','Original graph edges','Gray: no finite','Star: monocyte','Purple contour:','Frozen macrophage candidate','Mean AUCell per','Median AUCell per','C1-C4: display','Bold gene labels','Historical candidates','Within-cohort row','Colors compare','Original AUC values','cNMF program |','Within-pathway z-score','20 patients |','GPL: glycerophospholipid','F1 - Pyrimidine:','F1 - Cys','F1 - Arg','F1 - GPL:','26 paired patients','Pre-ranked GSEA','SPP1 expression and source contribution','Center: expression','Percentages below:','Colored outlines','2.11  |','2.11 |','2.9  |','2.9 |','Patient medians','Patients ordered by','31 profiles from','Display clusters of','Chromosomes','UMAP1','UMAP2','Patients','Reference','Tumour epithelial','200 cells','701 cells','Original Monocle3','C4 versus other mature macrophages','Cell type','Malignant epithelial','Regulon cluster']
def source_panel(z,fn,letter):
 src=Path(z['source']);crop=z.get('source_crop',[0,0,1,1]);doc=fitz.open(src);p=doc[0];rect=p.rect;l,t,rr,b=crop;clip=fitz.Rect(l*rect.width,t*rect.height,rr*rect.width,b*rect.height)
 # Crop only selected panel, retaining vector plots and raster heatmaps.
 c=fitz.open();pg=c.new_page(width=clip.width,height=clip.height);pg.show_pdf_page(pg.rect,doc,0,clip=clip)
 bad=[]
 for block in pg.get_text('dict')['blocks']:
  for line in block.get('lines',[]):
   txt=''.join(s['text'] for s in line['spans']).strip()
   if any(s in txt for s in remove_patterns) or len(txt)>78 or (fn==1 and letter=='D' and ('Original cNMF program' in txt or (txt.startswith('GSE') and '/' not in txt))) or (fn==1 and letter=='B' and txt.isdigit() and 1<=int(txt)<=22 and int(txt) not in [1,5,10,15,20]):
    bad.append(line['bbox']);removed.append(dict(figure=fn,panel=letter,text=txt))
 for box in bad:pg.add_redact_annot(fitz.Rect(box),fill=False)
 if bad:pg.apply_redactions(images=0,graphics=0)
 # Tight geometry excludes newly removed explanatory text.
 pix=pg.get_pixmap(matrix=fitz.Matrix(1,1),alpha=False);im=Image.frombytes('RGB',[pix.width,pix.height],pix.samples);bb=ImageChops.difference(im,Image.new('RGB',im.size,'white')).convert('L').point(lambda q:255 if q>15 else 0).getbbox()
 if bb:clip2=fitz.Rect(max(0,bb[0]-2),max(0,bb[1]-2),min(pg.rect.width,bb[2]+2),min(pg.rect.height,bb[3]+2))
 else:clip2=pg.rect
 return c,clip2

def retype(page,box,fn,letter):
 lines=[]
 for block in page.get_text('dict',clip=box)['blocks']:
  for line in block.get('lines',[]):
   for sp in line['spans']:
    if not sp['text'].strip():continue
    lines.append((sp,line['dir']))
 for sp,di in lines:page.add_redact_annot(fitz.Rect(sp['bbox']),fill=False)
 page.apply_redactions(images=0,graphics=0)
 for sp,di in lines:
  txt=sp['text'].replace('−','-');br=fitz.Rect(sp['bbox']);ox,oy=sp['origin'];fs=max(6.,min(8.,sp['size']))
  if fn==1 and letter=='D' and ('/' in txt and ('GSE' in txt or 'RCC' in txt or 'SS_' in txt)):
   label=f'P{len([x for x in keys if x["figure"]==1])+1:02d}';keys.append(dict(figure=1,panel='D',display=label,original=txt));txt=label
  txt=txt.replace('Normalized usage','F1 usage').replace('Within-method row z-score','Row z').replace('Relative density','Density').replace('Antigen processing and presentation','Antigen presentation')
  if len(txt)>40 and fn==2 and letter=='A':fs=6
  font='hebo' if ('Bold' in sp['font'] or 'bold' in sp['font']) else ('heit' if 'Oblique' in sp['font'] or 'Italic' in sp['font'] else 'helv')
  rgb=tuple(((sp['color']>>sh)&255)/255 for sh in [16,8,0]);tw=fitz.get_text_length(txt,fontname=font,fontsize=fs)
  angle=math.degrees(math.atan2(-di[1],di[0]));anchor=fitz.Point(ox,oy)
  if abs(angle)<1:
   x=(br.x0+br.x1-tw)/2;y=(br.y0+br.y1)/2+fs*.32
   x=max(box.x0+.5,min(x,box.x1-tw-.5))
   if fn==1 and letter=='D' and re.fullmatch(r'P[0-9]+',txt):x=box.x0+box.width*.70
   anchor=fitz.Point(x,y)
   page.insert_text(anchor,txt,fontsize=fs,fontname=font,color=rgb)
  else:
   # Preserve the original baseline direction of rotated gene and axis labels.
   page.insert_text(anchor,txt,fontsize=fs,fontname=font,color=rgb,morph=(anchor,fitz.Matrix(angle)))
for fn in [1,2,3]:
 doc=fitz.open();page=doc.new_page(width=210*MM,height=297*MM)
 arr={z['panel']:z for z in {1:old1,2:old2,3:old3}[fn]}
 for letter,title,x,y,w,h in layouts[fn]:
  z=arr[letter].copy();custom=d/'panels'/f'{fn}{letter}.pdf'
  if custom.exists():z.update(source=str(custom),source_crop=[0,0,1,1])
  elif fn==2 and letter=='B':z.update(source_crop=[.27,.048,.68,.282])
  elif fn==3 and letter=='C':z.update(source=str(r/'phase199_fig3_marker_GSEA_20260918_v5/panels/GSEA_curves_raw.pdf'),source_crop=[0,0,1,1])
  c,cl=source_panel(z,fn,letter);box=fitz.Rect(x*MM,y*MM,(x+w)*MM,(y+h)*MM);page.show_pdf_page(box,c,0,clip=cl)
  if not custom.exists():retype(page,box,fn,letter)
  records.append(dict(figure=fn,panel=letter,source=z['source'],crop=z.get('source_crop'),x=x,y=y,width=w,height=h,custom=custom.exists()))
 # Draw panel headings in an isolated graphics state after all panels.
 headings=fitz.open();hp=headings.new_page(width=210*MM,height=297*MM)
 from matplotlib import font_manager
 hf=font_manager.findfont(font_manager.FontProperties(family='DejaVu Sans',weight='bold'));hp.insert_font(fontname='HeaderFont',fontfile=hf)
 for letter,title,x,y,w,h in layouts[fn]:
  hp.insert_text((x*MM,(y-2)*MM),letter,fontsize=10,fontname='HeaderFont',color=(0,0,0));hp.insert_text(((x+5)*MM,(y-2)*MM),title,fontsize=7.5,fontname='HeaderFont',color=(0,0,0))
 if fn==1:
  hp.draw_rect(fitz.Rect(10*MM,195*MM,60*MM,202*MM),color=None,fill=(1,1,1))
  hp.draw_rect(fitz.Rect(183.6*MM,54*MM,192*MM,59*MM),color=None,fill=(1,1,1))
 page.show_pdf_page(page.rect,headings,0)
 # Explicit repeated letter uses a separate font resource to prevent missing first-glyph rendering.
 page.insert_font(fontname='PanelLetter',fontfile=hf)
 page.insert_text((10*MM,17*MM),'A',fontsize=10,fontname='PanelLetter',color=(0,0,0))
 out=d/'figures'/f'Main_Figure_{fn}_A4_v10.pdf';doc.save(out,garbage=4,deflate=True)
 page.get_pixmap(matrix=fitz.Matrix(1.8,1.8)).save(str(d/'qc'/f'Fig{fn}_preview.png'))
 print('PAGE_READY',fn,flush=True)
(d/'tables/PANEL_PROVENANCE.json').write_text(json.dumps(records,indent=2));(d/'tables/REMOVED_TEXT_FOR_CAPTIONS.json').write_text(json.dumps(removed,ensure_ascii=False,indent=2));(d/'tables/PATIENT_LABEL_KEY.json').write_text(json.dumps(keys,indent=2))
