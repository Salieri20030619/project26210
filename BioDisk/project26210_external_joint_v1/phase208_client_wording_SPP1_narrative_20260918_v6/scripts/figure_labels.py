from pathlib import Path
import json,io
import fitz,numpy as np
from pypdf import PdfReader,PdfWriter
from reportlab.pdfgen import canvas
r=Path.cwd();d=r/'phase208_client_wording_SPP1_narrative_20260918_v6';p=d/'report/report_materialized.json';b=json.loads(p.read_text())
for x in b:
 if 'text' in x:x['text']=x['text'].replace('因此，SPP1由分群后的标志鉴定自然进入后续研究，它对应的是分布于多个巨噬状态的表达信号。','这一共同表达特征使SPP1成为进一步解析巨噬细胞状态的重点分子，其信号来源同时涉及上述三个群。')
p.write_text(json.dumps(b,ensure_ascii=False,indent=2))
f=d/'figures/Figure4.pdf';doc=fitz.open(f);page=doc[0];recs=[]
for old,new in [('Primary','Model 1'),('Sensitivity','Model 2')]:
 boxes=page.search_for(old);assert len(boxes)==1;rect=boxes[0];recs.append((rect,new));page.add_redact_annot(rect+(-.1,-.1,.1,.1),fill=False)
page.apply_redactions(images=0,graphics=0);clean=d/'qc/Figure4_text_removed.pdf';doc.save(clean,garbage=4,deflate=True)
buf=io.BytesIO();c=canvas.Canvas(buf,pagesize=(page.rect.width,page.rect.height));c.setFont('Helvetica',6)
for rect,new in recs:c.drawString(rect.x0,page.rect.height-rect.y1+1.5,new)
c.save();buf.seek(0);wr=PdfWriter();pg=wr.add_page(PdfReader(clean).pages[0]);pg.merge_page(PdfReader(buf).pages[0]);tmp=io.BytesIO();wr.write(tmp);doc.close();f.write_bytes(tmp.getvalue())
old=fitz.open(r/'phase207_manuscript_figures_crossrefs_20260918_v5/figures/Figure4.pdf');new=fitz.open(f);a=old[0].get_pixmap(alpha=False);z=new[0].get_pixmap(alpha=False);aa=np.frombuffer(a.samples,np.uint8).reshape(a.height,a.width,3);zz=np.frombuffer(z.samples,np.uint8).reshape(z.height,z.width,3);delta=np.any(aa!=zz,axis=2);mask=np.zeros(delta.shape,bool);mask[753:775,518:558]=True;assert not np.any(delta&~mask),'Unexpected plot change outside legend'
assert not new[0].search_for('Sensitivity');new[0].get_pixmap(matrix=fitz.Matrix(2,2),clip=fitz.Rect(480,735,565,795)).save(str(d/'qc/legend.png'));print('LEGEND_READY; all pixels outside legend unchanged')
