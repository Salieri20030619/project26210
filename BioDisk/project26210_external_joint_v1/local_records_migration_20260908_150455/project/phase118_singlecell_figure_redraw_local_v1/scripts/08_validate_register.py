from pathlib import Path
import json,hashlib,csv,datetime
import pandas as pd
from PIL import Image
r=Path(__file__).resolve().parents[1];d=Path('C:/Users/31937/Downloads/26210/01_单细胞注释总览_v1')
m=json.loads((d/'核查记录/aggregate_validation.json').read_text())
c=pd.read_csv(r.parent/'phase43_part1_client_figures_local_v1/source_remote_v1/celltype_counts_v1.tsv',sep='\t').set_index('cell_type_broad_v1').n_cells.to_dict();assert c==m['counts']
pngs=list(d.rglob('*.png'));pdfs=list(d.rglob('*.pdf'));assert len(pngs)==len(pdfs)==19
for p in pngs:
 with Image.open(p) as im:im.verify()
for p in pdfs:assert p.read_bytes().startswith(b'%PDF')
q={'status':'PASS_FIRST_REVIEW_DELIVERY','png':19,'pdf':19,'counts_match':True,'cells':258120,'visual_review':'Overview includes all 13 markers; both strips and count figures inspected. Coordinate/labels preserved. User style review pending.','delivery':str(d)}
(r/'qc/delivery_validation.json').write_text(json.dumps(q,indent=2),encoding='utf8')
now=datetime.datetime.now().astimezone().isoformat(timespec='seconds')
with (r.parent/'RESULTS_MANIFEST.tsv').open('a',encoding='utf8',newline='') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n')
 for i,p in enumerate(sorted(pngs+pdfs)):
  w.writerow([f'phase118_review_figure_{i+1:02d}_v1','phase118_singlecell_redraw','review_figure',str(p),now,'PASS_DATA_AND_LAYOUT_REVIEW',f'SHA256 {hashlib.sha256(p.read_bytes()).hexdigest()}; counts 14/14 match','First user-review draft; frozen annotation and UMAP'])
text='''\n\n## Phase118 first single-cell figure review delivery (2026-09-08)\n\n- First batch delivered to user-authorized `C:\\Users\\31937\\Downloads\\26210\\01_单细胞注释总览_v1`: 19 PNG and 19 PDF files, Chinese reading note/captions and aggregate validation.\n- Includes overview, horizontal/vertical count plots, cell-identity UMAP, two marker strips and 13 standalone marker UMAPs.\n- Remote rendering retained frozen coordinates, annotations and all 258120 cells. All 14 cell-group counts match the existing local validated table. Only rendered figures and non-identifying aggregates were transferred.\n- Layout review covered overview/all 13 markers, both strips and both bar charts. User style approval remains pending; marker maps currently use expression scatter, not smoothed density contours.\n- Next: collect user feedback on this batch before proceeding to the next analysis figures. Every subsequent completed batch goes to Downloads\\26210 in a new versioned directory.\n'''
for name in ['PROJECT_STATUS.md','CODEX_HANDOFF.md']:
 with (r.parent/name).open('a',encoding='utf8') as f:f.write(text)
(d/'阅读说明与图注.txt').write_text((d/'阅读说明与图注.md').read_text(encoding='utf-8-sig'),encoding='utf-8-sig')
print(json.dumps(q,ensure_ascii=False))
