from pathlib import Path
import json,re,html,hashlib,shutil
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,Flowable,CondPageBreak
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.enums import TA_CENTER,TA_JUSTIFY
from reportlab.lib.pagesizes import A4
from reportlab.lib import colors
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from pypdf import PdfReader,PdfWriter,Transformation
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=r/'phase190_mechanism_focused_report_20260917_v1';base=r/'phase143_sections25_27_candidate_report_remote_v1'
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
S['cell']=ParagraphStyle('cell',fontName='Song',fontSize=8.3,leading=11.8,wordWrap='CJK')
S['cellc']=ParagraphStyle('cellc',parent=S['cell'],alignment=TA_CENTER)
S['cellh']=ParagraphStyle('cellh',parent=S['cell'],fontName='Hei')
S['cellhc']=ParagraphStyle('cellhc',parent=S['cellh'],alignment=TA_CENTER)
def clean(t):
 t=str(t).replace('−','-').replace('–','-').replace('—','-').replace('‑','-')
 supers=str.maketrans('⁻⁰¹²³⁴⁵⁶⁷⁸⁹','-0123456789')
 return re.sub(r'×10([⁻⁰¹²³⁴⁵⁶⁷⁸⁹]+)',lambda m:'e'+m.group(1).translate(supers),t)
def para(t,sty='p'):
 t=clean(t);latin='Times-Bold' if sty in ['title','h1','h2','caption','tablecap','cellh','cellhc'] else 'Times-Roman'
 rich=''.join('<font name="'+latin+'">'+html.escape(x)+'</font>' if x and all(32<=ord(c)<=126 for c in x) else html.escape(x).replace('\n','<br/>') for x in re.split(r'([\x20-\x7e]+)',t))
 return Paragraph(rich,S[sty])
fignames=['Figure1_annotation_composite_v1.pdf','Figure2_CNV_composite_v1.pdf','Figure3_cNMF_candidate_composite_v1.pdf','Figure4_pySCENIC_candidate_composite_v1.pdf','Figure5_protein_density_composite_v1.pdf','Figure6_metabolism_candidate_composite_v1.pdf','Figure7_trajectory_candidate_composite_v1.pdf']
figpaths=[]
for fn in fignames:
 pp=next(r.glob('phase*/figures/'+fn));dest=d/'figures'/fn
 if not dest.exists():shutil.copy2(pp,dest)
 figpaths.append(dest)
placements=[]
class Fig(Flowable):
 def __init__(self,n):
  super().__init__();self.n=n;pp=PdfReader(figpaths[n-1]).pages[0];self.width=493.27559055118115;self.height=float(pp.mediabox.height)*self.width/float(pp.mediabox.width)
 def draw(self):pass
 def drawOn(self,c,x,y,_sW=0):
  placements.append({'figure':self.n,'page':c.getPageNumber(),'x':x,'y':y,'width':self.width,'height':self.height});super().drawOn(c,x,y,_sW)
story=[];md=[];figi=0;tables=0
for i,x in enumerate(b):
 typ=x['type']
 if typ=='figure':
  figi+=1;story +=[Spacer(1,8),Fig(figi),Spacer(1,5)]
  md+=['',f'图{figi} [既有图件](../figures/{fignames[figi-1]})'];continue
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
  if x['caption'].startswith('附表 S1B'):weights=[3.8,.9,1.2,1.2]
  if x['caption'].startswith('表 11 '):weights=[2.0,.45,.65,1.2,1.6]
  left={j for j,h in enumerate(headers) if h in ['基因','关系与模型','关系','模型','证据来源','研究对象','观察结果','本报告用途','端点','有效范围','结果','统计或解释','通路','比较','类型','预测量','通路/轴','目标','蛋白','端点/方法','分析','关键参数和统计单位','客户报告的敲低变化','生信对应','证据定位','原校正口径','算法','调控子','队列','判定','方法']}
  data=[[para(v,'cellh' if j in left else 'cellhc') for j,v in enumerate(headers)]]
  data +=[[para(v,'cell' if j in left else 'cellc')for j,v in enumerate(row)] for row in rows]
  tt=Table(data,colWidths=[(A4[0]-90)*v/sum(weights) for v in weights],repeatRows=1,hAlign='LEFT')
  style=[('VALIGN',(0,0),(-1,-1),'MIDDLE'),('BACKGROUND',(0,0),(-1,0),colors.HexColor('#dfeaf0')),('GRID',(0,0),(-1,-1),.25,colors.HexColor('#D9D9D9')),('LEFTPADDING',(0,0),(-1,-1),5),('RIGHTPADDING',(0,0),(-1,-1),5),('TOPPADDING',(0,0),(-1,-1),4),('BOTTOMPADDING',(0,0),(-1,-1),4)]
  for j in range(1,len(data)):
   if j%2==0:style.append(('BACKGROUND',(0,j),(-1,j),colors.HexColor('#f2f5f6')))
  tt.setStyle(TableStyle(style));cap=para(x['caption'],'tablecap');cap.keepWithNext=len(rows)<=10
  story += ([] if len(rows)<=10 else [CondPageBreak(80)])+[cap,tt,Spacer(1,6)]
  md+=['',x['caption'],'| '+' | '.join(headers)+' |','| '+' | '.join(['---']*n)+' |']+['| '+' | '.join(map(str,row))+' |'for row in rows]
 else:
  if typ in ['h1','h2']:story.append(CondPageBreak(125))
  story.append(para(x['text'],typ));md+=['',{'title':'# ','h1':'## ','h2':'### '}.get(typ,'')+clean(x['text'])]
title=b[0]['text'];stem='26210_TYMP_MTA_SPP1_macrophage_integrated_report_v2'
(d/'report'/f'{stem}.md').write_text('\n'.join(md),encoding='utf-8')
def footer(c,doc):
 c.saveState();c.setFont('Song',9);c.drawCentredString(A4[0]/2,23,str(doc.page));c.restoreState()
tmp=d/'qc/layout_before_vectors.pdf'
doc=SimpleDocTemplate(str(tmp),pagesize=A4,leftMargin=45,rightMargin=45,topMargin=40,bottomMargin=39,title=title,author='项目26210')
doc.build(story,onFirstPage=footer,onLaterPages=footer)
writer=PdfWriter()
for i,page in enumerate(PdfReader(tmp).pages):
 for pl in placements:
  if pl['page']==i+1:
   fp=PdfReader(figpaths[pl['figure']-1]).pages[0];s=pl['width']/float(fp.mediabox.width);page.merge_transformed_page(fp,Transformation().scale(s).translate(pl['x'],pl['y']))
 writer.add_page(page)
writer.add_metadata({'/Title':title,'/Author':'项目26210','/Subject':'TYMP相关代谢 MTA候选与SPP1富集巨噬状态'})
out=d/'report'/f'{stem}.pdf'
with out.open('wb')as f:writer.write(f)
(d/'qc/PLACEMENTS.json').write_text(json.dumps(placements,indent=2))
(d/'qc/BUILD.json').write_text(json.dumps({'pages':len(writer.pages),'tables':tables,'figures':len(placements),'sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'bytes':out.stat().st_size},indent=2))
print('BUILD',len(writer.pages),'pages',tables,'tables',len(placements),'figures')
