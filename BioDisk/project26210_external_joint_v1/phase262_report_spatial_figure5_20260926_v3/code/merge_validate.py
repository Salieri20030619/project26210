
from pathlib import Path
import json,hashlib,re,zipfile
import pymupdf as fitz
from lxml import etree
P=Path(__file__).resolve().parents[1];R=P.parent
base=R/'phase208_client_wording_SPP1_narrative_20260918_v6/report/26210_TYMP_MTA_SPP1_four_figures_narrative_v6.pdf'
figure=R/'phase261_spatial_ternary_integrated_20260926_v5/delivery/26210_Spatial_A4_GSEA_Ternary_v5.pdf'
appendix=next((P/'qc/render').glob('*.pdf'))
b=fitz.open(base);a=fitz.open(appendix);f=fitz.open(figure)
assert len(b)==31 and len(a)==9 and len(f)==1
for i,page in enumerate(a):
 assert page.get_text().strip().splitlines()[-1]==str(32+i),(i,page.get_text()[-30:])
out=fitz.open();out.insert_pdf(b);out.insert_pdf(a)
# Replace rasterized Word figure with the unchanged native A4 vector source.
ip=36;pg=out[ip]
for xref in pg.get_contents():out.update_stream(xref,b'')
pg.show_pdf_page(pg.rect,f,0)
pg.insert_text((pg.rect.width/2-15,10),'Figure 5',fontname='helv',fontsize=8)
pg.insert_text((pg.rect.width/2-4,pg.rect.height-4),'37',fontname='tiro',fontsize=8)
linkcount=0
for page in out:
 for link in page.get_links():
  if link.get('kind')==fitz.LINK_GOTO and link.get('page')==37:
   link['page']=36;link['to']=fitz.Point(0,0);page.update_link(link);linkcount+=1
toc=b.get_toc()
toc += [[1,'空间研究与图5',32],[2,'空间研究方法',32],[2,'空间结果',34],[2,'TYMP与SPP1相关功能特征',36],[2,'图5 空间综合图',37],[2,'图5 图注',38],[2,'空间讨论',39],[2,'新增参考文献',40]]
out.set_toc(toc)
out.set_metadata(dict(title='26210 TYMP MTA SPP1 完整报告 含空间图5',author='',subject='原报告与空间续稿 图5及正文交叉引用'))
target=P/'delivery/26210_完整报告_含空间图5与交叉引用_v3.pdf'
out.save(target,garbage=4,deflate=True);out.close()
d=fitz.open(target)
checks=[]
for i in range(31):
 assert b[i].get_text()==d[i].get_text(),i
 assert b[i].get_pixmap(matrix=fitz.Matrix(1,1)).samples==d[i].get_pixmap(matrix=fitz.Matrix(1,1)).samples,i
checks.append('Original 31 pages retain identical text and rendered pixels')
for i in range(31,40):
 if i!=36:assert a[i-31].get_text()==d[i].get_text(),i
checks.append('All updated appendix text and page numbers retained')
assert linkcount>=14
for i in range(31,40):
 for lk in d[i].get_links():
  if lk['kind']==fitz.LINK_GOTO:assert lk['page']==36
checks.append(f'{linkcount} PDF link rectangles navigate to the native Figure5 page')
docx=next((P/'delivery').glob('*.docx'))
with zipfile.ZipFile(docx) as z:
 root=etree.fromstring(z.read('word/document.xml'))
 ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
 fields=root.xpath('//w:fldSimple/@w:instr',namespaces=ns)
 refs=[x for x in fields if 'REF fig5_spatial' in x]
 assert len(refs)==14,len(refs)
 assert len(root.xpath('//w:bookmarkStart[@w:name="fig5_spatial"]',namespaces=ns))==1
 assert any('SEQ Figure' in x and '\\r 5' in x for x in fields)
checks.append('DOCX Figure5 SEQ bookmark and 14 REF fields verified')
srcmd=next((R/'phase253_spatial_continue_numbering_20260925_v2/text').glob('*.md')).read_text()
newmd=(P/'text/26210_spatial_figure5_body_v3.md').read_text()
for para in srcmd.split('\n\n'):
 if para.strip():assert para.strip() in newmd,para[:100]
checks.append('Every prior spatial paragraph and heading retained verbatim; new text added without altering previous statistics')
build=json.loads((P/'qc/BUILD.json').read_text())
for record in build['inputs']:assert hashlib.sha256(Path(record['path']).read_bytes()).hexdigest()==record['sha256']
checks.append('All original input hashes unchanged')
for i in [36,37]:
 d[i].get_pixmap(matrix=fitz.Matrix(150/72,150/72)).save(P/f'qc/full_report_page_{i+1}.png')
(P/'qc/STRUCTURAL_VALIDATION.json').write_text(json.dumps(dict(status='PASS',pages=len(d),figure5_page=37,caption_page=38,body_ref_fields=14,pdf_ref_rectangles=linkcount,checks=checks),ensure_ascii=False,indent=2))
print(json.dumps(dict(pages=len(d),figure5_page=37,caption_page=38,links=linkcount,checks=checks),ensure_ascii=False))
