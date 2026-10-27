from pathlib import Path
import json,html,re,io,hashlib
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.pagesizes import A4
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas
from pypdf import PdfReader,PdfWriter
P=Path(__file__).resolve().parents[1];R=P.parent
source=R/'phase263_spatial_reintegrated_report_20260926_v1/delivery/26210_integrated_spatial_report_v1.pdf'
font=R/'phase143_sections25_27_candidate_report_remote_v1/fonts'
pdfmetrics.registerFont(TTFont('Song',str(font/'simsun.ttc'),subfontIndex=0));pdfmetrics.registerFont(TTFont('Hei',str(font/'simhei.ttf')))
styles={'title':ParagraphStyle('title',fontName='Hei',fontSize=13,leading=20,spaceAfter=15),'p':ParagraphStyle('p',fontName='Song',fontSize=10.2,leading=16,spaceAfter=10,wordWrap='CJK'),'note':ParagraphStyle('note',fontName='Song',fontSize=9.2,leading=14,spaceAfter=8,wordWrap='CJK')}
caption=json.loads((P/'report/Supplementary_Figure_S1_caption.json').read_text())
def rich(text,typ):
    latin='Times-Bold' if typ=='title' else 'Times-Roman'
    return ''.join('<font name="'+latin+'">'+html.escape(x)+'</font>' if x and all(32<=ord(c)<=126 for c in x) else html.escape(x) for x in re.split(r'([\x20-\x7e]+)',text))
story=[Paragraph(rich(t,typ),styles[typ]) for typ,t in caption]
cap=P/'qc/caption.pdf'
SimpleDocTemplate(str(cap),pagesize=A4,leftMargin=45,rightMargin=45,topMargin=44,bottomMargin=42).build(story)
cp=PdfReader(cap);assert len(cp.pages)==1, len(cp.pages)
old=PdfReader(source);writer=PdfWriter();writer.clone_document_from_reader(old)
figure=PdfReader(P/'figures/Supplementary_Figure_S1.pdf');assert len(figure.pages)==1
writer.add_page(figure.pages[0]);writer.add_page(cp.pages[0])
for page_number in [40,41]:
    buf=io.BytesIO();c=canvas.Canvas(buf,pagesize=A4);c.setFont('Song',9);c.drawCentredString(A4[0]/2,23,str(page_number));c.save();buf.seek(0)
    writer.pages[page_number-1].merge_page(PdfReader(buf).pages[0])
outline=writer.add_outline_item('补充图',39)
writer.add_outline_item('补充图 S1：细胞组成、代谢与 ATAC 足迹活性空间投射',39,parent=outline)
writer.add_outline_item('补充图 S1 图注',40,parent=outline)
writer.add_named_destination('suppfigS1',39);writer.add_named_destination('suppfigS1_caption',40)
writer.add_metadata({'/Title':'26210 整合报告（含补充图 S1）','/Author':'项目26210','/Subject':'原完整报告及补充图S1和图注'})
out=P/'delivery/26210_integrated_report_with_SuppFigS1_v1.pdf'
with out.open('wb') as f:writer.write(f)
new=PdfReader(out)
assert len(new.pages)==41
for i in range(39):
    assert old.pages[i].get_contents().get_data()==new.pages[i].get_contents().get_data(),i
    assert old.pages[i].extract_text()==new.pages[i].extract_text(),i
oldlinks=sum(len(pg.get('/Annots',[])) for pg in old.pages);newlinks=sum(len(pg.get('/Annots',[])) for pg in new.pages[:39]);assert oldlinks==newlinks
assert 'Supplementary Fig. S1' in new.pages[39].extract_text()
captext=new.pages[40].extract_text()
for token in ['HNF4A','CEBPD','JUNB','STAT3','80%','362/3199','1701/3880','500']:assert token in captext,token
oldb=json.loads((R/'phase263_spatial_reintegrated_report_20260926_v1/report/report_materialized.json').read_text())
newb=oldb+[dict(type='supplementary_figure',figure='S1',source=str(P/'figures/Supplementary_Figure_S1.pdf'))]+[dict(type='caption' if typ=='title' else typ,text=t) for typ,t in caption]
(P/'report/report_materialized.json').write_text(json.dumps(newb,ensure_ascii=False,indent=2))
oldmd=(R/'phase263_spatial_reintegrated_report_20260926_v1/report/26210_integrated_spatial_report_v1.md').read_text()
(P/'report/26210_integrated_report_with_SuppFigS1_v1.md').write_text(oldmd+'\n\n## 补充图\n\n![补充图 S1](../figures/Supplementary_Figure_S1.png)\n\n'+'\n\n'.join(t for _,t in caption))
qc=dict(status='PASS_VISUAL_PENDING',pages=41,original_pages_unchanged=39,main_figures=5,supplementary_figures=1,supplementary_figure_page=40,caption_page=41,original_annotations_retained=newlinks,caption_tokens_verified=True,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),output_sha256=hashlib.sha256(out.read_bytes()).hexdigest())
(P/'qc/BUILD_VALIDATION.json').write_text(json.dumps(qc,indent=2))
print(json.dumps(qc),flush=True)
