from pathlib import Path
import json,hashlib,shutil,datetime,re
import numpy as np,pandas as pd
from PIL import Image
from pypdf import PdfReader
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=r/'phase197_top20_network_dynamics_20260917_v4'
Image.MAX_IMAGE_PIXELS=None
tab=pd.read_csv(d/'tables/TOP20_RANKING_ALL85.tsv',sep='\t',index_col=0)
assert len(tab)==85 and tab.selected.sum()==20 and tab.patients.eq(16).all()
sel=json.loads((d/'tables/TOP20_SELECTION.json').read_text())
assert sel['paths']==tab.index[:20].tolist()
# Verify curve values against frozen display arrays.
dy=np.load(r/'phase192_sections28_212_figures_20260917_v1/inputs/DYNAMIC_DISPLAY.npz')
groups=json.loads((d/'tables/DYNAMIC_GENE_GROUPS.json').read_text())
curves=pd.read_csv(d/'tables/DYNAMIC_CURVES.tsv',sep='\t')
for group,genes in groups.items():
 for mode in ['primary','sensitivity']:
  vv=curves[(curves.panel==group)&(curves['mode']==mode)].sort_values('rank')
  expected=dy[mode][[list(dy['genes']).index(g) for g in genes]].mean(0)
  assert len(vv)==40 and np.allclose(vv.mean_logCP10K,expected,atol=1e-12)
# Verify all plotted network edges are actual GMT memberships.
gmt=pd.read_csv(r/'phase177_L2_myeloid_state_integration_20260916_v1/source_tables/program_genes.tsv',sep='\t')
edges=pd.read_csv(d/'tables/GENE_PATHWAY_EDGES.tsv',sep='\t')
for row in edges.itertuples():assert row.gene in gmt.loc[gmt.program==row.pathway,'gene_symbols'].iloc[0].split(';')
checks=[]
for f in (d/'figures').glob('*.pdf'):
 p=PdfReader(f);assert len(p.pages)==1
 im=Image.open(f.with_suffix('.png'));assert min(im.info['dpi'])>599
 checks.append({'file':f.name,'pixels':im.size,'dpi':im.info['dpi']})
norm=lambda s:re.sub(r'\s+','',s).lower()
text=norm(PdfReader(d/'figures/Main_Figure_2_top20_v4.pdf').pages[0].extract_text())
assert all(norm(p) in text for p in sel['paths'])
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
delivery=d/'delivery';delivery.mkdir(exist_ok=True)
for p in (d/'figures').iterdir():shutil.copy2(p,delivery/p.name)
for p in (d/'qc').glob('Main*.jpg'):shutil.copy2(p,delivery/p.name)
for p in (r/'phase196_selected_panels_20260917_v3/delivery').glob('Main_Figure_1*'):
 shutil.copy2(p,delivery/p.name);assert sha(p)==sha(delivery/p.name)
for p in [d/'REVISION_GUIDE_CN.md',d/'tables/TOP20_RANKING_ALL85.tsv',d/'tables/DYNAMIC_GENE_GROUPS.json']:shutil.copy2(p,delivery/p.name)
(delivery/'SHA256.json').write_text(json.dumps({p.name:sha(p) for p in delivery.iterdir() if p.name!='SHA256.json'},indent=2))
(d/'qc/FINAL_VALIDATION.json').write_text(json.dumps({'figures':checks,'top20_labels_present':True,'ranking_patients':16,'all85_retained':True,'network_edges_validated':len(edges),'curve_values_match_frozen':True,'new_tests':False,'visual_review':'Fig2/3 PDF renders reviewed; Fig2 compact tissue redraw rechecked'},indent=2))
inputs=list((d/'scripts').glob('*.py'))+list((d/'tables').glob('*'))
(d/'qc/CODE_AND_TABLE_SHA256.json').write_text(json.dumps({str(p):sha(p) for p in inputs},indent=2))
now=datetime.datetime.now().isoformat();(d/'record_backup').mkdir(exist_ok=True)
for n in ['PROJECT_STATUS.md','DECISIONS.md','CODEX_HANDOFF.md','RESULTS_MANIFEST.tsv']:
 if (r/n).exists():shutil.copy2(r/n,d/'record_backup'/n)
note=f'\n## phase197 | Top20及Fig3非热图重绘 | {now}\n按用户要求：Fig2A从85通路改为Top20，以16名患者内TYMP配对差值（先患者内等权程序，再跨患者中位数）绝对值排序，完整85排名保留。Fig3全部色块热图换为基因-功能归属网络、来源比例带及表达动态曲线，保留等高线和原始轨迹。网络边逐项验证GMT成员，曲线值与冻结数组一致，不新增检验。Fig1完整沿用phase196。2张修改PDF已视觉检查，600dpi PNG和SVG导出。正文未改。交付：{d}/delivery；本地Downloads/26210/27_Top20与Fig3非热图版_v4。\n'
for n in ['PROJECT_STATUS.md','DECISIONS.md','CODEX_HANDOFF.md']:
 with (r/n).open('a',encoding='utf8') as f:f.write(note)
with (r/'RESULTS_MANIFEST.tsv').open('a') as f:f.write(f'\nphase197_top20_network\tphase197\tfigures\t{d}/delivery\t{now}\tPASS\tTop20 ranking; GMT edges; frozen curves; PDF visual QA\tFig2/3 revised; Fig1 unchanged; outside report\n')
(d/'STATUS.json').write_text(json.dumps({'status':'complete','changed_figures':[2,3],'unchanged_figure':1,'report_modified':False},indent=2))
print('READY',len(list(delivery.iterdir())),'files')
