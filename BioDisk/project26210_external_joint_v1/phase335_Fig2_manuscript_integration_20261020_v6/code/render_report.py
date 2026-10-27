from pathlib import Path
import json,re,html,hashlib,shutil
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,Flowable,CondPageBreak,PageBreak
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.enums import TA_CENTER,TA_JUSTIFY
from reportlab.lib.pagesizes import A4
from reportlab.lib import colors
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from pypdf import PdfReader,PdfWriter,Transformation
import logging
logging.getLogger('pypdf').setLevel(logging.ERROR)
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=r/'phase335_Fig2_manuscript_integration_20261020_v6';base=r/'phase143_sections25_27_candidate_report_remote_v1'
b=json.loads((d/'report/report_materialized.json').read_text())
pdfmetrics.registerFont(TTFont('Song',str(base/'fonts/simsun.ttc'),subfontIndex=0));pdfmetrics.registerFont(TTFont('Hei',str(base/'fonts/simhei.ttf')))
S={}
S['title']=ParagraphStyle('title',fontName='Hei',fontSize=16,leading=23,alignment=TA_CENTER,spaceAfter=9,keepWithNext=True)
S['date']=ParagraphStyle('date',fontName='Song',fontSize=10,leading=15,alignment=TA_CENTER,spaceAfter=15,keepWithNext=True)
S['h1']=ParagraphStyle('h1',fontName='Hei',fontSize=13,leading=19,spaceBefore=12,spaceAfter=6,keepWithNext=True)
S['h2']=ParagraphStyle('h2',fontName='Hei',fontSize=11.2,leading=16.5,spaceBefore=9,spaceAfter=4,keepWithNext=True)
S['p']=ParagraphStyle('p',fontName='Song',fontSize=10.4,leading=15.8,alignment=TA_JUSTIFY,firstLineIndent=20.8,spaceAfter=5,wordWrap='CJK',allowWidows=0,allowOrphans=0)
S['caption']=ParagraphStyle('caption',fontName='Hei',fontSize=9.7,leading=14,spaceBefore=7,spaceAfter=5)
S['tablecap']=ParagraphStyle('tablecap',parent=S['caption'],keepWithNext=False)
S['note']=ParagraphStyle('note',fontName='Song',fontSize=8.6,leading=12.3,spaceAfter=7,wordWrap='CJK')
S['ref']=ParagraphStyle('ref',parent=S['p'],fontSize=9.5,leading=14,alignment=0,firstLineIndent=0,spaceAfter=6)
S['cell']=ParagraphStyle('cell',fontName='Song',fontSize=8.3,leading=11.8,wordWrap='CJK')
S['cellc']=ParagraphStyle('cellc',parent=S['cell'],alignment=TA_CENTER)
S['cellh']=ParagraphStyle('cellh',parent=S['cell'],fontName='Hei')
S['cellhc']=ParagraphStyle('cellhc',parent=S['cellh'],alignment=TA_CENTER)
def clean(t):
 t=str(t).replace('−','-').replace('–','-').replace('—','-').replace('‑','-')
 supers=str.maketrans('⁻⁰¹²³⁴⁵⁶⁷⁸⁹','-0123456789')
 return re.sub(r'×10([⁻⁰¹²³⁴⁵⁶⁷⁸⁹]+)',lambda m:'e'+m.group(1).translate(supers),t)

panel_defs=json.loads((r/'phase206_all_mainfigures_A4_20260918_v10/tables/PANEL_PROVENANCE.json').read_text())
figpaths={n:d/'figures'/f'Figure{n}.pdf' for n in [1,2,3,4,5,6]};figpaths['S1']=d/'figures/Supplementary_Figure_S1.pdf';figpaths['S2']=d/'figures/Supplementary_Figure_S2.pdf'
placements=[]
valid={1:set('ABC'),2:set('ABCDEF'),3:set('ABCDEF'),4:set('ABCDE'),5:set('ABCDEFGH'),6:set('ABCDEFGHI'),'S1':set('ABC'),'S2':set('ABCD')}
def rich_plain(t,latin):
 return ''.join('<font name="'+latin+'">'+html.escape(x)+'</font>' if x and all(32<=ord(c)<=126 for c in x) else html.escape(x).replace('\n','<br/>') for x in re.split(r'([\x20-\x7e]+)',t))
def para(t,sty='p'):
 if t=='公开组织bulk PRJNA1223824':t='公开组织bulk\nPRJNA1223824'
 if re.match(r'^\[\d+\]',t):
  sty='ref';t=re.sub(r'\s+(https?://)',r'\n\1',t)
 t=clean(t);latin='Times-Bold' if sty in ['title','h1','h2','caption','tablecap','cellh','cellhc'] else 'Times-Roman'
 parts=re.split(r'((?:补充图\s*S[12]|图\s*[1-6])[A-I]?)',t);rich=[]
 for part in parts:
  m=re.fullmatch(r'(?:补充图\s*(S[12])|图\s*([1-6]))([A-I]?)',part)
  ss=rich_plain(part,latin)
  if m and sty not in ['caption','tablecap','cellh','cellhc']:
   n=m[1] if m[1] else int(m[2]);p=m[3];assert not p or p in valid[n]
   ss='<link href="#fig'+str(n)+p+'" color="#245C7D">'+ss+'</link>'
  rich.append(ss)
 p=Paragraph(''.join(rich),S[sty])
 if sty in ['h1','h2']:p.outline_title=t;p.outline_level=0 if sty=='h1' else 1
 return p
class FullFigure(Flowable):
 def __init__(self,n):super().__init__();self.n=n;self.width=1;self.height=700
 def draw(self):
  c=self.canv;n=self.n;c.bookmarkPage(f'fig{n}');c.addOutlineEntry(f'补充图 {n}' if isinstance(n,str) else f'图{n}',f'fig{n}',level=0 if isinstance(n,str) else 1,closed=False)
  if n=='S2':
   anchors=[('A',10,29),('B',110,29),('C',10,145),('D',110,145)]
  elif n=='S1':
   anchors=[('A',3,22),('B',3,111),('C',3,200)]
  elif n==5:
   anchors=[('A',9,17),('B',9,73),('C',9,145),('D',116,145),('E',9,197),('F',116,197),('G',9,244),('H',116,244)]
  elif n==6:
   anchors=[('A',3,3),('B',75,3),('C',150,3),('D',3,83),('E',105,83),('F',3,130),('G',3,178),('H',3,220),('I',140,220)]
  elif n==2:
   anchors=[('A',10,19),('B',10,106),('C',76,103),('D',142,103),('E',10,202),('F',120,202)]
  elif n==1:
   anchors=[('A',10,24),('B',10,85),('C',10,170)]
  else:anchors=[(z['panel'],z['x'],z['y']) for z in panel_defs if z['figure']==n-1]
  for letter,x,y in anchors:c.bookmarkPage(f'fig{n}{letter}') if n==6 else c.bookmarkHorizontalAbsolute(f'fig{n}{letter}',A4[1]-(y-5)*72/25.4,left=x*72/25.4)
  placements.append(dict(figure=n,page=c.getPageNumber()))
story=[];md=[];tables=0
for i,x in enumerate(b):
 typ=x['type']
 if typ=='h1' and x['text']=='补充图':
  md+=['','## 补充图'];continue
 if typ in ['figure','supplementary_figure']:
  n=x['figure'];story +=[PageBreak(),FullFigure(n),PageBreak()];md +=['',f'<a id="fig{n}"></a>',f'![图{n}](figures/Figure{n}.png)'];continue
 if typ=='caption' and i+1<len(b) and b[i+1]['type']=='table' and x['text']==b[i+1]['caption']:continue
 if typ=='table':
  tables+=1;headers=x['headers'];rows=x['rows'];n=len(headers)
  weights=[1]*n
  if n==2:weights=[1.25,4.75]
  elif n==3:weights=[2.0,3.3,1.8]
  elif n==4:weights=[2.0,1.7,2.5,2.7]
  elif n==5:weights=[2.6,0.6,1.0,1.0,1.2]
  elif n==6:weights=[1.5,.5,.9,1.05,.95,.95]
  elif n==7:weights=[1.4,1.0,.5,.9,.9,1,1]
  elif n==8:weights=[1.25,1.3,1.1,.5,.95,.8,1,1]
  if x.get('_col_weights'):weights=x['_col_weights']
  if x['caption'].startswith('附表 S1B'):weights=[3.8,.9,1.2,1.2]
  if x['caption'].startswith('附表S18 '):weights=[2.0,.45,.65,1.2,1.6]
  if x['caption'].startswith(('附表S19 ','附表S20 ')):weights=[2.7,.8,2.3,1.0,1.0]
  left={j for j,h in enumerate(headers) if h in ['基因','关系与模型','关系','模型','证据来源','研究对象','观察结果','证据层次','端点','有效范围','结果','统计或解释','通路','比较','类型','预测量','通路/轴','目标','蛋白','端点/方法','分析','关键参数和统计单位','TYMP敲低后的变化','公共多组学对应','证据定位','原校正口径','算法','调控子','队列','判定','方法']}
  data=[[para(v,'cellh' if j in left else 'cellhc') for j,v in enumerate(headers)]]
  data +=[[para(v,'cell' if j in left else 'cellc')for j,v in enumerate(row)] for row in rows]
  assert len(weights)==n and all(v>0 for v in weights), (x['caption'],n,weights)
  tt=Table(data,colWidths=[(A4[0]-90)*v/sum(weights) for v in weights],repeatRows=1,hAlign='LEFT')
  style=[('VALIGN',(0,0),(-1,-1),'MIDDLE'),('BACKGROUND',(0,0),(-1,0),colors.HexColor('#dfeaf0')),('GRID',(0,0),(-1,-1),.25,colors.HexColor('#D9D9D9')),('LEFTPADDING',(0,0),(-1,-1),5),('RIGHTPADDING',(0,0),(-1,-1),5),('TOPPADDING',(0,0),(-1,-1),4),('BOTTOMPADDING',(0,0),(-1,-1),4)]
  for j in range(1,len(data)):
   if j%2==0:style.append(('BACKGROUND',(0,j),(-1,j),colors.HexColor('#f2f5f6')))
  if x['caption'].startswith('附表 S12'):style += [('TOPPADDING',(0,0),(-1,-1),3),('BOTTOMPADDING',(0,0),(-1,-1),3)]
  tt.setStyle(TableStyle(style));cap=para(x['caption'],'tablecap');cap.keepWithNext=len(rows)<=10
  story += ([] if len(rows)<=10 else [CondPageBreak(80)])+[cap,tt,Spacer(1,6)]
  if x.get('_immune_table_id') and x.get('source'):story.append(para(x['source'],'note'))
  md+=['',x['caption'],'| '+' | '.join(headers)+' |','| '+' | '.join(['---']*n)+' |']+['| '+' | '.join(map(str,row))+' |'for row in rows]
 else:
  if typ=='h1' and x['text']=='1 研究设计与分析方法':story.append(PageBreak())
  if typ in ['h1','h2']:story.append(CondPageBreak(270 if x['text'].startswith('4.15 ') else 80 if x['text'].startswith('2.26 ') else 125))
  pobj=para(x['text'],typ)
  if x.get('keep_with_next'):pobj.keepWithNext=True
  story.append(pobj);md+=['',{'title':'# ','h1':'## ','h2':'### '}.get(typ,'')+clean(x['text'])]

title=b[0]['text'];stem='26210_manuscript_Fig2_updated_v6'
(d/'report'/f'{stem}.md').write_text('\n'.join(md),encoding='utf-8')
def footer(c,doc):
 c.saveState();c.setFont('Song',9);c.drawCentredString(A4[0]/2,23,str(doc.page));c.restoreState()
tmp=d/'qc/layout_with_links.pdf'
class ReportDoc(SimpleDocTemplate):
 def afterFlowable(self,flow):
  if hasattr(flow,'outline_title'):
   self._outline_serial=getattr(self,'_outline_serial',0)+1
   key='section_'+str(self._outline_serial)
   self.canv.bookmarkPage(key)
   self.canv.addOutlineEntry(flow.outline_title,key,level=flow.outline_level,closed=False)
doc=ReportDoc(str(tmp),pagesize=A4,leftMargin=45,rightMargin=45,topMargin=40,bottomMargin=39,title=title,author='项目26210')
doc.build(story,onFirstPage=footer,onLaterPages=footer)
# Clone the full PDF to retain the internal links, destinations and figure bookmarks.
reader=PdfReader(tmp);writer=PdfWriter();writer.clone_document_from_reader(reader)
import io
from reportlab.pdfgen import canvas
for pl in placements:
 pg=writer.pages[pl['page']-1];fp=PdfReader(figpaths[pl['figure']]).pages[0];w=float(fp.mediabox.width);h=float(fp.mediabox.height)
 assert w<=A4[0]+.1 and h<=A4[1]+.1
 dx=(A4[0]-w)/2;dy=(A4[1]-h)/2;pg.merge_transformed_page(fp,Transformation().translate(dx,dy))
 buf=io.BytesIO();ca=canvas.Canvas(buf,pagesize=A4);ca.setFont('Song',9);ca.drawCentredString(A4[0]/2,23,str(pl['page']));ca.save();buf.seek(0);pg.merge_page(PdfReader(buf).pages[0])
 pl.update(source=str(figpaths[pl['figure']]),scale=1.,x=dx,y=dy)
writer.add_metadata({'/Title':title,'/Author':'项目26210','/Subject':'TYMP相关代谢 SPP1巨噬免疫与空间联系 含补充图S1-S2'})
out=d/'report'/f'{stem}.pdf'
with out.open('wb') as f:writer.write(f)
(d/'qc/PLACEMENTS.json').write_text(json.dumps(placements,indent=2));(d/'qc/BUILD.json').write_text(json.dumps(dict(pages=len(writer.pages),tables=tables,figures=8,bytes=out.stat().st_size,sha256=hashlib.sha256(out.read_bytes()).hexdigest()),indent=2));print('BUILD',len(writer.pages),'pages',tables,'tables',placements,flush=True)
