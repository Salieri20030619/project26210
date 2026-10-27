from pathlib import Path
import json,hashlib
from pypdf import PdfReader,PdfWriter,Transformation
w=Path('C:/codex_test/tmp/fig140')
r=PdfReader(w/'report_before_vector.pdf')
placements=json.loads((w/'placement.json').read_text())['placements']
figures={1:Path('C:/codex_test/output/26210/Figure1_annotation/Figure1_annotation_composite_v1.pdf'),2:Path('C:/codex_test/output/26210/Figure2_CNV/Figure2_CNV_composite_v1.pdf')}
writer=PdfWriter()
for i,p in enumerate(r.pages):
 for q in placements:
  if q['page']==i+1:
   fp=PdfReader(figures[q['figure_id']]).pages[0]
   scale=q['width']/float(fp.mediabox.width)
   assert q['y']>30 and q['y']+q['height']<float(p.mediabox.height)
   p.merge_transformed_page(fp,Transformation().scale(scale).translate(q['x'],q['y']))
 writer.add_page(p)
out=Path('C:/codex_test/output/26210/26210_TYMP_report_with_Fig1_Fig2_v3.pdf')
with out.open('wb') as f:writer.write(f)
check=PdfReader(out)
assert len(check.pages)==20
assert '2.1 ' in check.pages[4].extract_text()
assert 'Chromosome-level CNV signal' in check.pages[5].extract_text()
assert '2.2 ' not in check.pages[4].extract_text()
assert '2.2 ' in ''.join(p.extract_text() for p in check.pages[6:8])
print(json.dumps(dict(pages=len(check.pages),placements=placements,sha256=hashlib.sha256(out.read_bytes()).hexdigest())))
