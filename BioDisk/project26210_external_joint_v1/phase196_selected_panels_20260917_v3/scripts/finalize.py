from pathlib import Path
import json,hashlib,shutil,datetime,re
from PIL import Image
from pypdf import PdfReader
import pandas as pd
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=r/'phase196_selected_panels_20260917_v3';prev=r/'phase195_mainfigures_compact_20260917_v2'
Image.MAX_IMAGE_PIXELS=None
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
a=pd.read_csv(d/'tables/DOTPLOT_ORIGINAL_SUMMARY.tsv',sep='\t')
b=pd.read_csv(r/'phase131_sections2_5_2_6_reference_correction_remote_v1/data/PATHWAY_DOTPLOT_SUMMARY.tsv',sep='\t')
pd.testing.assert_frame_equal(a,b,check_exact=False,atol=1e-14,rtol=1e-14)
reg=json.loads((d/'tables/DOTPLOT_VALIDATION.json').read_text())
p2=PdfReader(d/'figures/Main_Figure_2_column_v3.pdf')
norm=lambda s:re.sub(r'\s+','',s).lower()
text=norm(p2.pages[0].extract_text())
missing=[p for p in reg['pathways'] if norm(p) not in text];assert not missing,missing
qc=[]
for p in sorted((d/'figures').glob('*.pdf')):
 rr=PdfReader(p);assert len(rr.pages)==1
 im=Image.open(p.with_suffix('.png'));assert min(im.info['dpi'])>599
 qc.append({'file':p.name,'dpi':im.info['dpi'],'size':im.size})
delivery=d/'delivery';delivery.mkdir(exist_ok=True)
for p in (d/'figures').iterdir():shutil.copy2(p,delivery/p.name)
for p in (d/'qc').glob('Main*.jpg'):shutil.copy2(p,delivery/p.name)
# Include the unchanged third main figure as a convenience; no re-export.
for p in (prev/'delivery').glob('Main_Figure_3*'):
 shutil.copy2(p,delivery/p.name);assert sha(p)==sha(delivery/p.name)
shutil.copy2(d/'LAYOUT_GUIDE_CN.md',delivery/'LAYOUT_GUIDE_CN.md')
(d/'qc/FINAL_VALIDATION.json').write_text(json.dumps({'files':qc,'all_85_labels_in_final_pdf':True,'dotplot_values_unchanged':True,'fig3_byte_identical':True,'visual_review':'Both PDFs reviewed; final full-column layout rechecked','report_modified':False},indent=2))
sources=list((d/'inputs').glob('*'))+list((d/'scripts').glob('*.py'))
(d/'qc/SOURCE_CODE_SHA256.json').write_text(json.dumps({str(p):sha(p) for p in sources},indent=2))
(delivery/'SHA256.json').write_text(json.dumps({p.name:sha(p) for p in delivery.iterdir() if p.name!='SHA256.json'},indent=2))
now=datetime.datetime.now().isoformat();(d/'record_backup').mkdir(exist_ok=True)
for n in ['PROJECT_STATUS.md','DECISIONS.md','CODEX_HANDOFF.md','RESULTS_MANIFEST.tsv']:
 if (r/n).exists():shutil.copy2(r/n,d/'record_backup'/n)
note=f'\n## phase196 | 指定图件与整列逻辑重排 | {now}\nFig1D替换为用户指定H_patient_program_expression_blocks_v2；保留指定调控网络。Fig2A用用户指定85通路点图的原数据及绘图函数适配左侧完整一列，全部85标签在最终PDF检出。右侧按恶性上皮轨迹、MTA生化背景、C4发现、组织配对联系排列；移出重复动态/标志热图，原件保留。Fig3沿用phase195，交付文件SHA256完全一致。2张修改PDF及SVG/600dpiPNG/预览完成检查，无新检验，正文未变。交付：{d}/delivery；本地Downloads/26210/26_三张主图逻辑重排_v3。\n'
for n in ['PROJECT_STATUS.md','DECISIONS.md','CODEX_HANDOFF.md']:
 with (r/n).open('a',encoding='utf8') as f:f.write(note)
with (r/'RESULTS_MANIFEST.tsv').open('a') as f:f.write(f'\nphase196_selected_panels\tphase196\tfigures\t{d}/delivery\t{now}\tPASS\t85 pathway labels; data equality; visual QA\tFig1/2 revised; Fig3 unchanged; outside report\n')
(d/'STATUS.json').write_text(json.dumps({'status':'complete','changed_figures':[1,2],'unchanged_figure':3,'report_modified':False},indent=2))
print('DELIVERY_READY',len(list(delivery.iterdir())))
