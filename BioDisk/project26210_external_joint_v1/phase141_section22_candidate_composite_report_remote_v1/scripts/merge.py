from pathlib import Path
import json,hashlib
from pypdf import PdfReader,PdfWriter,Transformation
w=Path('C:/codex_test/tmp/fig141');r=PdfReader(w/'report_before_vector.pdf')
placements=json.loads((w/'placement.json').read_text())['placements']
b=Path('C:/codex_test/output/26210')
figures={1:b/'Figure1_annotation/Figure1_annotation_composite_v1.pdf',2:b/'Figure2_CNV/Figure2_CNV_composite_v1.pdf',3:b/'Figure3_cNMF/Figure3_cNMF_candidate_composite_v1.pdf'}
writer=PdfWriter()
for i,p in enumerate(r.pages):
 for q in placements:
  if q['page']==i+1:
   fp=PdfReader(figures[q['figure_id']]).pages[0];scale=q['width']/float(fp.mediabox.width)
   p.merge_transformed_page(fp,Transformation().scale(scale).translate(q['x'],q['y']))
 writer.add_page(p)
out=b/'26210_TYMP_report_with_Fig1_Fig2_Fig3_draft_v4.pdf'
with out.open('wb') as f:writer.write(f)
check=PdfReader(out)
assert len(check.pages)==22
assert '2.2 ' in check.pages[6].extract_text() and '2.4 ' not in check.pages[6].extract_text()
assert 'cNMF programs, F1 usage and TYMP' in check.pages[7].extract_text()
assert '2.4 ' in ''.join(p.extract_text() for p in check.pages[8:10])
print(json.dumps(dict(pages=len(check.pages),placements=placements,sha256=hashlib.sha256(out.read_bytes()).hexdigest())))
