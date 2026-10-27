from pathlib import Path
import json,hashlib,shutil,datetime
import numpy as np,pandas as pd
from PIL import Image
from pypdf import PdfReader
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=r/'phase195_mainfigures_compact_20260917_v2'
Image.MAX_IMAGE_PIXELS=None
z=np.load(r/'part2_rebuild_v1/04_markers_annotation/run_v4_macrophage_only/paired_patient_effects/cluster_4.npz',allow_pickle=True)
tab=pd.read_csv(d/'tables/FUNCTIONAL_GENE_PAIRS.tsv',sep='\t',index_col=0)
assert np.allclose(tab.values,z['effects'][:,[list(z['genes']).index(g) for g in tab.index]].T)
assert tab.columns.tolist()==list(z['patient'])
a=pd.read_csv(d/'tables/SPP1_CONTRIBUTION.tsv',sep='\t').drop(columns='display')
b=pd.read_csv(r/'phase177_L2_myeloid_state_integration_20260916_v1/source_tables/c4_contribution.tsv',sep='\t')
pd.testing.assert_frame_equal(a,b,check_exact=False,atol=1e-14,rtol=1e-14)
qc=[]
for f in sorted((d/'figures').glob('*.pdf')):
 p=PdfReader(f);assert len(p.pages)==1
 im=Image.open(f.with_suffix('.png'));assert min(im.info['dpi'])>599
 qc.append({'file':f.name,'size_mm':[float(p.pages[0].mediabox.width)*25.4/72,float(p.pages[0].mediabox.height)*25.4/72],'pixels':im.size,'dpi':im.info['dpi']})
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
prior=json.loads((r/'phase194_mainfigure_assembly_20260917_v1/tables/SOURCE_SHA256.json').read_text())
assert all(sha(Path(p))==h for p,h in prior.items())
netexact=sha(d/'inputs/user_selected_regulon_network.pdf')==sha(r/'phase130_section2_4_BE_correction_remote_v1/figures/06_E_gene_to_gene_regulon_network_v2.pdf')
(d/'qc/FINAL_VALIDATION.json').write_text(json.dumps({'figures':qc,'gene_pairs_exact':True,'contribution_table_exact':True,'network_user_file_matches_original':netexact,'network_nodes':71,'network_edges':71,'original_report_and_figures_unchanged':True,'visual_review':'All three; final clipped study title corrected and rechecked'},indent=2))
sources=list((d/'scripts').glob('*.py'))+list((d/'inputs').glob('*'))+list((d/'tables').glob('*'))
(d/'qc/SOURCE_AND_CODE_SHA256.json').write_text(json.dumps({str(p):sha(p) for p in sources},indent=2))
delivery=d/'delivery';delivery.mkdir(exist_ok=True)
for f in (d/'figures').iterdir():shutil.copy2(f,delivery/f.name)
for f in (d/'qc').glob('Main*.jpg'):shutil.copy2(f,delivery/f.name)
shutil.copy2(d/'REVISION_GUIDE_CN.md',delivery/'REVISION_GUIDE_CN.md')
hashes={p.name:sha(p) for p in delivery.iterdir() if p.name!='SHA256.json'}
(delivery/'SHA256.json').write_text(json.dumps(hashes,indent=2))
now=datetime.datetime.now().isoformat();(d/'record_backup').mkdir(exist_ok=True)
for n in ['PROJECT_STATUS.md','DECISIONS.md','CODEX_HANDOFF.md','RESULTS_MANIFEST.tsv']:
 if (r/n).exists():shutil.copy2(r/n,d/'record_backup'/n)
note=f'\n## phase195 | 三主图紧凑版及Fig3扩充 | {now}\n用户要求减少留白、复用原脚本、扩充Fig3，并指定加入06_E_gene_to_gene_regulon_network_v2.pdf。已完成3张v2主图：重绘蛋白2x2、原网络、功能与动态矩阵并按内容紧凑拼版。Fig1纳入71节点71边网络，节点位置/原绘图函数保持；Fig3新增26患者24基因配对热图和16患者C4/SPP1RNA贡献图，合计6组面板。无新增模型或推断检验；原正文和原图SHA256均保持。实际PDF渲染逐图复核，600dpi PNG及SVG已导出。产物：{d}/delivery。本地Downloads/26210/25_三张主图紧凑修订_v2。\n'
for n in ['PROJECT_STATUS.md','DECISIONS.md','CODEX_HANDOFF.md']:
 with (r/n).open('a',encoding='utf8') as f:f.write(note)
with (r/'RESULTS_MANIFEST.tsv').open('a') as f:f.write(f'\nphase195_compact_mainfigures\tphase195\tfigures\t{d}/delivery\t{now}\tPASS\t3 PDFs inspected; raw-value and source-hash checks; 600dpi\toutside report; Fig3 expanded; user network included\n')
(d/'STATUS.json').write_text(json.dumps({'status':'complete','main_figures':3,'panel_groups':20,'report_modified':False,'network_original_byte_match':netexact,'new_fig3_panels':2},indent=2))
print('READY',len(hashes),'files','network_original_byte_match',netexact)
