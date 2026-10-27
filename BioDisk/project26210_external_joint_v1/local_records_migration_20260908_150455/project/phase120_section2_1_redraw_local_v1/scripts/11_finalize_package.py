from pathlib import Path
import json,hashlib,csv,datetime
import pandas as pd,networkx as nx
from PIL import Image
r=Path(__file__).resolve().parents[2];w=Path(__file__).resolve().parents[1];d=Path('C:/Users/31937/Downloads/26210/02_2.1恶性上皮亚克隆_v1')
png=list(d.glob('*.png'));pdf=list(d.glob('*.pdf'));assert len(png)==len(pdf)==13
for p in png:
 with Image.open(p) as im:im.verify()
for p in pdf:assert p.read_bytes().startswith(b'%PDF')
m=json.loads((d/'核查记录/remote_aggregate_validation_v2.json').read_text());assert sum(m['status_counts'].values())==62652
s=r/'phase32_part1_m03_cnv_clone_trees_local_v1';edges=pd.read_csv(s/'results_local_v1/m03_tree_edges_v1.tsv',sep='\t');nodes=pd.read_csv(s/'results_local_v1/m03_tree_nodes_cnmf_annotation_v1.tsv',sep='\t');trees={}
for pt in ['GSE207493__RCC104','GSE207493__RCC86','GSE156632__RCC4','GSE156632__RCC5']:
 e=edges[edges.infercnv_unit==pt];n=nodes[nodes.infercnv_unit==pt];g=nx.DiGraph();g.add_edges_from(e[['parent_node','child_node']].itertuples(index=False,name=None));assert nx.is_arborescence(g) and set(g)==set(n.node_id);assert g.in_degree('NEUTRAL_ROOT')==0;trees[pt]={'nodes':len(n),'edges':len(e),'cells':int(n.clone_cells.sum())}
q={'status':'PASS_REVIEW_PACKAGE','png':13,'pdf':13,'epithelial_cells':62652,'tumour_epithelial_cells':38670,'conservative_malignant_cells':11950,'four_trees':trees,'patient_effect_reproduced':True,'visual_review':'All main figure families inspected; heatmap v2 and patient contrast v2 supersede initial drafts; individual trees use shared inspected layout.','scope':'Section2.1 first review package, pending user approval. Clinical radar has no matching formal groups and was not fabricated; trajectory belongs to section2.7.'}
(w/'qc/FINAL_VALIDATION_v1.json').write_text(json.dumps(q,indent=2),encoding='utf8');now=datetime.datetime.now().astimezone().isoformat(timespec='seconds')
with (r/'RESULTS_MANIFEST.tsv').open('a',encoding='utf8',newline='') as f:
 writer=csv.writer(f,delimiter='\t',lineterminator='\n')
 for i,p in enumerate(sorted(png+pdf)):
  writer.writerow([f'phase120_section21_fig_{i+1:02d}_v1','phase120_section21_redraw','review_figure',str(p),now,'PASS_REVIEW_CHECKS',f'SHA256 {hashlib.sha256(p.read_bytes()).hexdigest()}','Source-preserving redraw; original patient P/FDR retained'])
s='''\n\n## Phase120 section2.1 figure redraw (2026-09-08)\n\n- User accepted Phase119 v3 single-cell annotation figures and requested the first malignant-epithelial section2.1.\n- Delivered 13 PNG + 13 PDF review figures to `C:\\Users\\31937\\Downloads\\26210\\02_2.1恶性上皮亚克隆_v1`: continuous inferCNV heatmap (final filename v2), six-status epithelial UMAP, light/dark marker maps, CNV group and clone distributions, four-tree composite + four individual trees with event summaries, H-style density maps, and patient-level TYMP clone contrast (final filename v2).\n- Reference images 2/3/4/6/7/8 visually inspected. The CNV heatmap was regenerated from the server checkpoint, not a reused PNG. Frozen cell labels, UMAP coordinates, clone topology and thresholds remain unchanged.\n- Verified 62652 epithelial cells, 38670 tumour epithelial cells and11950 conservative malignant cells; all local clone rows matched source CNV rows. Four tree graphs are valid root-anchored arborescences with exact source node/edge sets.\n- Patient F1-high clone effect median 0.199935300713434 reproduced; n=14, 9/14 positive, P=0.02099609375, BH q=0.073486328125. Nominal association only; no causal TYMP-driven malignant progression claim.\n- Final plot families visually reviewed; initial heatmap/contrast drafts preserved separately. Chinese readme maps figures to source requirements and retains CNV discordance and interpretation boundaries.\n- Clinical exposure/resistance radar lacks a matching formal grouping and was not fabricated; full malignant trajectory is scheduled under original section2.7.\n- Current breakpoint: user reviews section2.1 first version. Keep every new batch/version in Downloads\\26210. No active plot jobs remain.\n'''
for name in ['PROJECT_STATUS.md','CODEX_HANDOFF.md']:
 with (r/name).open('a',encoding='utf8') as f:f.write(s)
with (r/'DECISIONS.md').open('a',encoding='utf8') as f:f.write('\n\n## 2026-09-08 approved single-cell style and section2.1 scope\n\nUser accepted the Phase119 v3 dendrogram dotplots and major-class dashed UMAP. Proceeded to first malignant section2.1, retaining source-based figures and original patient-level evidence. Section2.1 package is now awaiting visual review; all batches must be delivered to Downloads\\26210.\n')
print(json.dumps(q,ensure_ascii=False))
