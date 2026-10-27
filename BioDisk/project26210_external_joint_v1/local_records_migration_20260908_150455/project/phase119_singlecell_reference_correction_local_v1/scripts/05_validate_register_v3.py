from pathlib import Path
import json,csv,hashlib,datetime
from PIL import Image
import pandas as pd
r=Path('D:/project26210_external_joint_v1');d=Path('C:/Users/31937/Downloads/26210/01_单细胞注释总览_v3');meta=json.loads((d/'核查记录/aggregate_validation_v3.json').read_text())
old=pd.read_csv(r/'phase43_part1_client_figures_local_v1/source_remote_v1/celltype_counts_v1.tsv',sep='\t').set_index('cell_type_broad_v1').n_cells.to_dict();assert old==meta['cell_type_counts'];assert len(meta['genes'])==38
files=sorted(list(d.glob('*.png'))+list(d.glob('*.pdf')));assert len(files)==6
for f in d.glob('*.png'):
 with Image.open(f) as im:im.verify()
q={'status':'PASS_DATA_AND_VISUAL_REVIEW','figures':3,'formats':['PNG','PDF'],'cells':258120,'groups':14,'markers':38,'visual':'All three v3 PNGs inspected; vertical labels and brackets separated; tree-to-axis alignment confirmed. User review pending.'}
(r/'phase119_singlecell_reference_correction_local_v1/qc/FINAL_VALIDATION_v3.json').write_text(json.dumps(q,indent=2),encoding='utf8');now=datetime.datetime.now().astimezone().isoformat(timespec='seconds')
with (r/'RESULTS_MANIFEST.tsv').open('a',encoding='utf8',newline='') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n')
 for i,p in enumerate(files):w.writerow([f'phase119_corrected_figure_{i+1}_v3','phase119_reference_format_correction','review_figure',str(p),now,'PASS_DATA_AND_LAYOUT',f'SHA256 {hashlib.sha256(p.read_bytes()).hexdigest()}','38-marker tree dotplot or dashed-major-class UMAP; no changed cell labels'])
s='''\n\n## Phase119 reference-format correction (2026-09-08)\n\n- User clarified that the requested horizontal/vertical figure is a dendrogram-bearing marker dotplot, not a cell-count bar chart. UMAP must show dashed major-class outlines and grouped legends. Reference image was enlarged and inspected.\n- Delivered corrected v3 to `C:\\Users\\31937\\Downloads\\26210\\01_单细胞注释总览_v3`: horizontal and vertical 38-marker/14-group dotplots with computed average-linkage trees, plus UMAP with dashed major-class density outlines. Three PNG and three PDF files; all three PNGs visually inspected.\n- Frozen 258120 cells, UMAP coordinates and cell labels retained; 14/14 source counts match. No per-cell data transferred. Contours are display guides, not new biological classification.\n- v2 is an intermediate layout; v3 fixes axis ordering, vertical group-label crowding and line weights. Use v3 for current review. Phase118 count figures are not fulfillment of the corrected dotplot request.\n- Next breakpoint: user visual review of v3. Continue delivering each versioned batch to Downloads\\26210.\n'''
for name in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md']:
 with (r/name).open('a',encoding='utf8') as f:f.write(s)
print(json.dumps(q))
