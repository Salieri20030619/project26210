from pathlib import Path
import json,hashlib
from pypdf import PdfReader,PdfWriter,Transformation
w=Path('C:/codex_test/tmp/figfinal');b=Path('C:/codex_test/output/26210');r=PdfReader(w/'report_before_vector.pdf')
placements=json.loads((w/'placement.json').read_text())['placements']
files=['Figure1_annotation/Figure1_annotation_composite_v1.pdf','Figure2_CNV/Figure2_CNV_composite_v1.pdf','Figure3_cNMF/Figure3_cNMF_candidate_composite_v1.pdf','Figure4_pySCENIC/Figure4_pySCENIC_candidate_composite_v1.pdf','Figure5_protein/Figure5_protein_density_composite_v1.pdf','Figure6_metabolism/Figure6_metabolism_candidate_composite_v1.pdf','Figure7_trajectory/Figure7_trajectory_candidate_composite_v1.pdf']
writer=PdfWriter()
for i,p in enumerate(r.pages):
 for q in placements:
  if q['page']==i+1:
   fp=PdfReader(b/files[q['figure_id']-1]).pages[0];s=q['width']/float(fp.mediabox.width)
   assert abs(float(fp.mediabox.height)*s-q['height'])<.01
   p.merge_transformed_page(fp,Transformation().scale(s).translate(q['x'],q['y']))
 writer.add_page(p)
out=b/'26210_TYMP_report_selected_Fig1_to_Fig7_draft_v6.pdf'
with out.open('wb') as f:writer.write(f)
check=PdfReader(out);assert len(check.pages)==27 and len(placements)==7
info=dict(pages=27,figure_pages=[q['page'] for q in placements],sha256=hashlib.sha256(out.read_bytes()).hexdigest(),bytes=out.stat().st_size)
(w/'FINAL_PDF.json').write_text(json.dumps(info,indent=2));print(json.dumps(info))
