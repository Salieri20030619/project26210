from pathlib import Path
import json,re,shutil,sys,subprocess,hashlib,html,io,ast
import pymupdf as fitz
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,PageBreak,CondPageBreak,Flowable,KeepTogether
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from PIL import Image,ImageDraw
w=Path(__file__).resolve().parent;inp=w/'inputs';q=w/'qc';q.mkdir(exist_ok=True)
delivery=w/'delivery';delivery.mkdir(exist_ok=True)
mainname='26210_Main_Manuscript_v7.pdf';supname='26210_Supplement_v3.pdf'
old=json.loads(Path(r'C:\Users\31937\tmp\26210_manuscript_v6\report_v6.json').read_text(encoding='utf-8'))
b=json.loads(json.dumps(old));word_changes=[]
replacement={
'相互条件调整与程序基因敏感性分析':'相互条件调整、共享基因删除与逐基因删除',
'程序敏感性分析分别移除':'程序基因分析分别移除',
'描述性敏感性分析':'描述性逐基因及联合删基因分析',
'进行敏感性分析':'进行模型比较',
'用于敏感性分析':'用于去共享基因与附加协变量比较',
'敏感性宽度':'补充比较宽度',
'试验敏感性':'试验分层与CM025子集比较',
'糖酵解删基因敏感性':'糖酵解逐基因及联合删基因结果',
'主克隆划分对距离阈值的敏感性':'不同距离阈值下的主克隆划分变化',
'临床交互与敏感性结果':'临床交互、试验分层及CM025子集结果',
'生存关联及敏感性':'生存关联、试验分层与患者子集分析',
'及完整敏感性保留于配套clinical_statistics文件夹':'及试验分层、患者子集和逐患者遗漏结果保留于配套统计数据包E13',
'逐基因敏感性':'逐基因删除结果',
'同一队列的内部敏感性分析':'同一队列的试验分层与CM025子集分析',
'分割敏感性':'不同分割方法下的结果',
}
def rewrite(x,path=''):
 if isinstance(x,str):
  before=x
  for a,c in replacement.items():x=x.replace(a,c)
  if x!=before:word_changes.append({'path':path,'before':before,'after':x})
  return x
 if isinstance(x,list):return [rewrite(y,path+f'/{j}') for j,y in enumerate(x)]
 if isinstance(x,dict):return {k:rewrite(v,path+'/'+k) if k in ['text','caption','rows','headers'] else v for k,v in x.items()}
 return x
b=rewrite(b)
assert '敏感性' not in ''.join(str(x.get('text',x.get('rows',''))) for x in b)
start=next(j for j,x in enumerate(b) if x.get('text')=='4 统计附表与可复现信息')
refs=next(j for j,x in enumerate(b) if x.get('type')=='h1' and x.get('text')=='参考文献')
supstart=next(j for j,x in enumerate(b) if x.get('type')=='h1' and x.get('text')=='补充图')
appendix=b[start:refs];suppfig=b[supstart+1:];main=b[:start]+[{'type':'p','text':'完整统计附表S1–S20、补充图S1–S11及扩展分析收录于配套补充材料；正文中的附表与补充图引用均对应独立补充册。'}]+b[refs:supstart]
assert len([x for x in main if x['type']=='table'])==13
assert len([x for x in appendix if x['type']=='table'])==28
for a,c in zip([x for x in old if x['type']=='table'],[x for x in main if x['type']=='table']+[x for x in appendix if x['type']=='table']):
 # Numeric values and all statistics unchanged; only one descriptive source-cell wording is updated.
 assert re.findall(r'[-+]?\d*\.?\d+(?:e[-+]?\d+)?',str(a['rows']))==re.findall(r'[-+]?\d*\.?\d+(?:e[-+]?\d+)?',str(c['rows']))
mainroot=w/'main'
for s in ['report','figures','qc','code']:(mainroot/s).mkdir(parents=True,exist_ok=True)
for f in inp.glob('Figure*.pdf'):shutil.copy2(f,mainroot/'figures'/f.name)
(mainroot/'report/report_materialized.json').write_text(json.dumps(main,ensure_ascii=False,indent=2),encoding='utf-8')
(q/'MOVED_APPENDIX_BLOCKS.json').write_text(json.dumps(appendix+suppfig,ensure_ascii=False,indent=2),encoding='utf-8')
(q/'WORDING_CHANGES.json').write_text(json.dumps(word_changes,ensure_ascii=False,indent=2),encoding='utf-8')
s=(inp/'render_original.py').read_text(encoding='utf-8')
s=re.sub(r"r=Path\('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1'\);d=r/'[^']+';base=r/'[^']+'",f"r=Path({w.as_posix()!r});d=r/'main';base=r/'inputs'",s)
s=s.replace("str(base/'fonts/simsun.ttc')","r'C:\\Windows\\Fonts\\simsun.ttc'").replace("str(base/'fonts/simhei.ttf')","r'C:\\Windows\\Fonts\\simhei.ttf'")
s=s.replace("r/'phase206_all_mainfigures_A4_20260918_v10/tables/PANEL_PROVENANCE.json'","r/'inputs/PANEL_PROVENANCE.json'")
s=s.replace("if m and sty not in", "if m and not m[1] and sty not in")
s=s.replace('26210_manuscript_Fig2_updated_v6','26210_Main_Manuscript_v7').replace('figures=8','figures=len(placements)').replace('含补充图S1-S2','补充材料另册')
class TextEncoding(ast.NodeTransformer):
 def visit_Call(self,node):
  self.generic_visit(node)
  if isinstance(node.func,ast.Attribute) and node.func.attr in ['read_text','write_text'] and not any(k.arg=='encoding' for k in node.keywords):node.keywords.append(ast.keyword(arg='encoding',value=ast.Constant('utf-8')))
  return node
s=ast.unparse(ast.fix_missing_locations(TextEncoding().visit(ast.parse(s))))
(mainroot/'code/render.py').write_text(s,encoding='utf-8')
run=subprocess.run([sys.executable,str(mainroot/'code/render.py')],capture_output=True,text=True,encoding='utf-8',creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
print(run.stdout,run.stderr,flush=True);run.check_returncode()

# Supplement: original statistics + S1/S2 + expanded results + S3-S11.
pdfmetrics.registerFont(TTFont('Song',r'C:\Windows\Fonts\simsun.ttc',subfontIndex=0));pdfmetrics.registerFont(TTFont('Hei',r'C:\Windows\Fonts\simhei.ttf'))
styles={k:ParagraphStyle(k,fontName='Hei' if k in ['title','h1','h2','caption','cellh'] else 'Song',fontSize=fs,leading=ld,spaceAfter=sp,wordWrap='CJK',keepWithNext=k in ['h1','h2']) for k,fs,ld,sp in [('title',17,25,12),('h1',13,19,8),('h2',11,16,6),('body',10.2,15.5,7),('p',10.2,15.5,7),('note',8.6,12.3,6),('caption',9.7,14,5),('cell',8.3,11.8,0),('cellh',8.3,11.8,0)]}
def P(t,k='body'):
 t=str(t);assert '敏感性' not in t and '结果解释与适用范围' not in t
 pp=Paragraph(html.escape(t).replace('\n','<br/>'),styles[k]);return pp
toc=[];targets={};placements=[]
def heading(t,level=1,key=None):
 p=P(t,'h1' if level==1 else 'h2');p.tag=(t,level,key);return p
class SupDoc(SimpleDocTemplate):
 def afterFlowable(self,f):
  if hasattr(f,'tag'):
   title,level,key=f.tag;toc.append([level,title,self.page])
   if key:targets.setdefault(key,self.page-1)
class FullFigure(Flowable):
 def __init__(self,n):super().__init__();self.n=n;self.width=1;self.height=700
 def draw(self):
  placements.append({'figure':self.n,'page':self.canv.getPageNumber()-1});targets['fig'+self.n]=self.canv.getPageNumber()-1;toc.append([1,'补充图'+self.n,self.canv.getPageNumber()])
story=[P('TYMP相关代谢、SPP1巨噬状态及基质免疫联系','title'),P('补充材料 v3｜与正文v7配套','h1'),P('本册统一收录原正文的统计附表S1–S20、方法与数据说明、补充图S1/S2，以及扩展分析和补充图S3–S11。原图表编号保留，表中样本量、方向、效应量、P值和FDR不变。'),P('内容目录','h1'),P('第一部分：统计附表与方法、数据说明（原第4部分）。\n第二部分：补充图S1、S2及图注。\n第三部分：27个扩展分析模块，保留分析对象、方法、主要结果和数据来源。\n第四部分：补充图S3–S11及图注。'),PageBreak(),heading('第一部分  统计附表与方法、数据说明',1)]
def add_table(x):
 headers=x['headers'];rows=x['rows'];n=len(headers)
 weights={2:[1.25,4.75],3:[2,3.3,1.8],4:[2,1.7,2.5,2.7],5:[2.6,.6,1,1,1.2],6:[1.5,.5,.9,1.05,.95,.95],7:[1.4,1,.5,.9,.9,1,1],8:[1.25,1.3,1.1,.5,.95,.8,1,1]}.get(n,[1]*n)
 if x.get('_col_weights'):weights=x['_col_weights']
 if x['caption'].startswith('附表 S1B'):weights=[3.8,.9,1.2,1.2]
 if x['caption'].startswith('附表S18 '):weights=[2,.45,.65,1.2,1.6]
 if x['caption'].startswith(('附表S19 ','附表S20 ')):weights=[2.7,.8,2.3,1,1]
 data=[[P(t,'cellh') for t in headers]]+[[P(t,'cell') for t in row] for row in rows]
 table=Table(data,colWidths=[505*v/sum(weights) for v in weights],repeatRows=1,hAlign='LEFT')
 table.setStyle(TableStyle([('VALIGN',(0,0),(-1,-1),'TOP'),('BACKGROUND',(0,0),(-1,0),colors.HexColor('#DFEAF0')),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#F2F5F6')]),('GRID',(0,0),(-1,-1),.25,colors.HexColor('#D9D9D9')),('LEFTPADDING',(0,0),(-1,-1),5),('RIGHTPADDING',(0,0),(-1,-1),5),('TOPPADDING',(0,0),(-1,-1),4),('BOTTOMPADDING',(0,0),(-1,-1),4)]))
 cap=P(x['caption'],'caption');m=re.search(r'附表\s*(S\d+[AB]?)',x['caption']);cap.tag=(x['caption'],2,'table'+m[1] if m else None);cap.keepWithNext=True
 story.extend([cap,table,Spacer(1,7)])
for ai,x in enumerate(appendix[1:],start=1):
 if x['type']=='caption' and ai+1<len(appendix) and appendix[ai+1]['type']=='table' and x['text']==appendix[ai+1]['caption']:continue
 if x['type']=='table':add_table(x)
 elif x['type']=='h2':story.extend([CondPageBreak(100),heading(re.sub(r'^4\.(\d+)\s*',r'A\1  ',x['text']),2)])
 else:
  t=x.get('text','').replace('配套免疫统计与来源资料','配套统计数据包的E20–E25及BASE目录')
  story.append(P(t,'note' if x['type']=='note' else 'body'))
# S1/S2 begin on their graphic pages; avoid a title-only divider page.
for x in suppfig:
 if x['type']=='supplementary_figure':story +=[PageBreak(),FullFigure(x['figure']),PageBreak()]
 elif x['type']=='caption':story.append(heading(x['text'],2))
 else:story.append(P(x.get('text',''),'note' if x['type']=='note' else 'body'))
mods=json.loads(Path(r'C:\Users\31937\tmp\26210_supplement_full\qc\FINAL_MODULES.json').read_text(encoding='utf-8'))
story +=[PageBreak(),heading('第三部分  扩展分析',1)]
for m in mods:
 begin_len=len(story)
 story.extend([CondPageBreak(120),heading(m['id']+'  '+m['title'],2,m['id']),P('分析对象与方法','h2'),P(m['method']),P('主要结果','h2'),P(m['result'])])
 rows=[[P(t,'cellh') for t in ['比较/读出','实际单位或范围','原结果与统计']]]+[[P(t,'cell') for t in row] for row in m['rows']]
 tb=Table(rows,colWidths=[122,132,251],repeatRows=1);tb.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#DFEAF0')),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#F2F5F6')]),('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),5),('RIGHTPADDING',(0,0),(-1,-1),5),('TOPPADDING',(0,0),(-1,-1),4),('BOTTOMPADDING',(0,0),(-1,-1),4)]))
 story.extend([tb,Spacer(1,7),P('数据来源：Phase '+', '.join(map(str,m['phases']))+'；完整统计见配套统计数据包'+m['id']+'目录。','note')])
 story[begin_len:]=[KeepTogether(story[begin_len:])]
story +=[PageBreak(),heading('迁移附表中的文献编号',1),P('迁移内容保留原文献编号，完整来源如下；该列表与正文相同，避免独立阅读时缺少来源。','note')]
for x in b[refs+1:supstart]:story.append(P(x['text'],'note'))
def footer(c,doc):
 c.setFont('Song',8);c.setFillColor(colors.HexColor('#637582'));c.drawString(45,A4[1]-24,'26210｜补充材料 v3');c.drawCentredString(A4[0]/2,22,'SM-'+str(doc.page))
layout=w/'supplement_layout.pdf';SupDoc(str(layout),pagesize=A4,leftMargin=45,rightMargin=45,topMargin=43,bottomMargin=39).build(story,onFirstPage=footer,onLaterPages=footer)
sup=fitz.open(layout)
for pl in placements:
 source=fitz.open(inp/f"Supplementary_Figure_{pl['figure']}.pdf");sup[pl['page']].show_pdf_page(sup[pl['page']].rect,source,0)
for k in list(targets):
 if k.startswith('table'):
  baseid=re.sub('[AB]$','',k);targets.setdefault(baseid,targets[k])
# Unchanged graphic pages S3-S11, after deleting the 27 explanatory subsections from the prose.
prev=fitz.open(r'C:\Users\31937\Downloads\26210\Phase338_扩充补充材料_v2\26210_Expanded_Supplement_v2.pdf')
for idx in range(17,26):
 sup.insert_pdf(prev,from_page=idx,to_page=idx);n='S'+str(idx-14);targets['fig'+n]=len(sup)-1;toc.append([1,'补充图'+n,len(sup)])
for k,p in enumerate(sup):
 p.draw_rect(fitz.Rect(220,p.rect.height-37,375,p.rect.height-10),fill=(1,1,1),color=None);p.insert_text((278,p.rect.height-20),f'SM-{k+1}',fontsize=8)
sup.set_toc(toc)
mainpdf=fitz.open(mainroot/'report'/mainname);mainpl=json.loads((mainroot/'qc/PLACEMENTS.json').read_text());mainfig={str(x['figure']):x['page']-1 for x in mainpl}
links=[]
def add_refs(doc,kind):
 for pn,page in enumerate(doc):
  text=page.get_text();seen=set()
  patterns=[(r'附表\s*(S\d+[AB]?)','table'),(r'补充图\s*(S\d+)([A-I]?)','fig')]
  if kind=='supp':patterns.append((r'(?<!补充)图\s*([1-6])([A-I]?)','mainfig'))
  for pat,typ in patterns:
   for mt in re.finditer(pat,text):
    token=mt.group(0);target=mainfig.get(mt[1]) if typ=='mainfig' else targets.get(typ+mt[1])
    if target is None:continue
    for rect in page.search_for(token):
     key=(tuple(rect),typ,target)
     if key in seen:continue
     seen.add(key)
     if typ=='mainfig':ln={'kind':fitz.LINK_GOTOR,'from':rect,'file':mainname,'page':target,'to':fitz.Point(0,0)}
     elif kind=='main':ln={'kind':fitz.LINK_GOTOR,'from':rect,'file':supname,'page':target,'to':fitz.Point(0,0)}
     else:ln={'kind':fitz.LINK_GOTO,'from':rect,'page':target,'to':fitz.Point(0,0)}
     # Replace any obsolete link in this span, retaining main-figure links elsewhere.
     for prior in page.get_links():
      if prior['from'].intersects(rect):page.delete_link(prior)
     page.insert_link(ln);links.append({'document':kind,'source_page':pn+1,'text':token,'target_page':target+1,'target_file':ln.get('file','self')})
add_refs(mainpdf,'main');add_refs(sup,'supp')
mainpdf.save(delivery/mainname,garbage=4,deflate=True);sup.save(delivery/supname,garbage=4,deflate=True)
alltext='\n'.join(p.get_text() for p in sup);assert '结果解释与适用范围' not in alltext and '敏感性' not in alltext
assert '4 统计附表与可复现信息' not in '\n'.join(p.get_text() for p in mainpdf)
(q/'LINK_MAP.json').write_text(json.dumps(links,ensure_ascii=False,indent=2),encoding='utf-8')
(q/'TARGET_PAGES.json').write_text(json.dumps(targets,ensure_ascii=False,indent=2),encoding='utf-8')
(q/'CONTENT_VALIDATION.json').write_text(json.dumps({'main_pages':len(mainpdf),'supplement_pages':len(sup),'main_tables':13,'moved_tables':28,'original_total_tables':41,'main_figures':6,'supplement_figures':11,'removed_explanation_sections':27,'statistical_values_unchanged':True,'cross_document_links':sum(x['target_file']!='self' for x in links),'all_refs':len(links),'prohibited_terms':0},indent=2),encoding='utf-8')
for label,doc in [('main',mainpdf),('supp',sup)]:
 (q/label).mkdir(exist_ok=True)
 for k,p in enumerate(doc):p.get_pixmap(matrix=fitz.Matrix(1,1)).save(q/label/f'page_{k+1:02}.png')
 for st in range(0,len(doc),6):
  sheet=Image.new('RGB',(1230,1740),'#eeeeee')
  for j,k in enumerate(range(st,min(st+6,len(doc)))):
   im=Image.open(q/label/f'page_{k+1:02}.png').convert('RGB');im.thumbnail((400,815));x=(j%3)*410;y=(j//3)*870;sheet.paste(im,(x,y+30));ImageDraw.Draw(sheet).text((x+8,y+8),f'{label} page {k+1}',fill='black')
  sheet.save(q/label/f'contact_{st//6+1}.png')
print('COMPLETE',len(mainpdf),len(sup),'links',len(links),flush=True)
