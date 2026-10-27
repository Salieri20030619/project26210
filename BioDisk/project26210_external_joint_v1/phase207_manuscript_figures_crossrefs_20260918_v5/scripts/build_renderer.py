from pathlib import Path
r=Path.cwd();d=r/'phase207_manuscript_figures_crossrefs_20260918_v5';s=(r/'phase201_three_candidates_C4_main_report_20260918_v4/scripts/render_report.py').read_text();prefix=s[:s.index('def para(')].replace("d=r/'phase201_three_candidates_C4_main_report_20260918_v4'","d=r/'phase207_manuscript_figures_crossrefs_20260918_v5'").replace('Flowable,CondPageBreak','Flowable,CondPageBreak,PageBreak')
mid=r'''
panel_defs=json.loads((r/'phase206_all_mainfigures_A4_20260918_v10/tables/PANEL_PROVENANCE.json').read_text())
figpaths=[d/'figures'/f'Figure{n}.pdf' for n in [1,2,3,4]]
placements=[]
valid={1:set('ABC'),2:set('ABCDEFG'),3:set('ABCDEF'),4:set('ABCDE')}
def rich_plain(t,latin):
 return ''.join('<font name="'+latin+'">'+html.escape(x)+'</font>' if x and all(32<=ord(c)<=126 for c in x) else html.escape(x).replace('\n','<br/>') for x in re.split(r'([\x20-\x7e]+)',t))
def para(t,sty='p'):
 t=clean(t);latin='Times-Bold' if sty in ['title','h1','h2','caption','tablecap','cellh','cellhc'] else 'Times-Roman'
 parts=re.split(r'(图\s*[1-4][A-G]?)',t);rich=[]
 for part in parts:
  m=re.fullmatch(r'图\s*([1-4])([A-G]?)',part)
  ss=rich_plain(part,latin)
  if m and sty not in ['caption','tablecap','cellh','cellhc']:
   n=int(m[1]);p=m[2];assert not p or p in valid[n]
   ss='<link href="#fig'+str(n)+p+'" color="#245C7D">'+ss+'</link>'
  rich.append(ss)
 return Paragraph(''.join(rich),S[sty])
class FullFigure(Flowable):
 def __init__(self,n):super().__init__();self.n=n;self.width=1;self.height=700
 def draw(self):
  c=self.canv;n=self.n;c.bookmarkPage(f'fig{n}');c.addOutlineEntry(f'图{n}',f'fig{n}',level=0,closed=False)
  if n==1:
   anchors=[('A',10,24),('B',10,85),('C',10,170)]
  else:anchors=[(z['panel'],z['x'],z['y']) for z in panel_defs if z['figure']==n-1]
  for letter,x,y in anchors:c.bookmarkHorizontalAbsolute(f'fig{n}{letter}',A4[1]-(y-5)*72/25.4,left=x*72/25.4)
  placements.append(dict(figure=n,page=c.getPageNumber()))
story=[];md=[];tables=0
for i,x in enumerate(b):
 typ=x['type']
 if typ=='figure':
  n=x['figure'];story +=[PageBreak(),FullFigure(n),PageBreak()];md +=['',f'<a id="fig{n}"></a>',f'![图{n}](figures/Figure{n}.png)'];continue
'''
loop=s[s.index(" if typ=='caption' and"):s.index("title=b[0]['text'];stem=")]
tail=r'''
title=b[0]['text'];stem='26210_TYMP_MTA_SPP1_four_figures_crossrefs_v5'
(d/'report'/f'{stem}.md').write_text('\n'.join(md),encoding='utf-8')
def footer(c,doc):
 c.saveState();c.setFont('Song',9);c.drawCentredString(A4[0]/2,23,str(doc.page));c.restoreState()
tmp=d/'qc/layout_with_links.pdf'
doc=SimpleDocTemplate(str(tmp),pagesize=A4,leftMargin=45,rightMargin=45,topMargin=40,bottomMargin=39,title=title,author='项目26210')
doc.build(story,onFirstPage=footer,onLaterPages=footer)
# Clone the full PDF to retain the internal links, destinations and figure bookmarks.
reader=PdfReader(tmp);writer=PdfWriter();writer.clone_document_from_reader(reader)
import io
from reportlab.pdfgen import canvas
for pl in placements:
 pg=writer.pages[pl['page']-1];fp=PdfReader(figpaths[pl['figure']-1]).pages[0];w=float(fp.mediabox.width);h=float(fp.mediabox.height)
 assert w<=A4[0]+.1 and h<=A4[1]+.1
 dx=(A4[0]-w)/2;dy=(A4[1]-h)/2;pg.merge_transformed_page(fp,Transformation().translate(dx,dy))
 buf=io.BytesIO();ca=canvas.Canvas(buf,pagesize=A4);ca.setFont('Song',9);ca.drawCentredString(A4[0]/2,23,str(pl['page']));ca.save();buf.seek(0);pg.merge_page(PdfReader(buf).pages[0])
 pl.update(source=str(figpaths[pl['figure']-1]),scale=1.,x=dx,y=dy)
writer.add_metadata({'/Title':title,'/Author':'项目26210','/Subject':'四图整合与可点击面板交叉引用'})
out=d/'report'/f'{stem}.pdf'
with out.open('wb') as f:writer.write(f)
(d/'qc/PLACEMENTS.json').write_text(json.dumps(placements,indent=2));(d/'qc/BUILD.json').write_text(json.dumps(dict(pages=len(writer.pages),tables=tables,figures=4,bytes=out.stat().st_size,sha256=hashlib.sha256(out.read_bytes()).hexdigest()),indent=2));print('BUILD',len(writer.pages),'pages',tables,'tables',placements,flush=True)
'''
(d/'scripts/render_report.py').write_text(prefix+mid+loop+tail)
