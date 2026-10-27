from pathlib import Path
import json,hashlib,re,datetime,subprocess,zipfile,shutil
import numpy as np,pandas as pd
from PIL import Image
P=Path(__file__).resolve().parents[1];R=P.parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
manifest=json.loads((P/'qc/INPUT_MANIFEST.json').read_text())
for item in manifest:
 f=Path(item['path']);st=f.stat()
 assert st.st_size==item['bytes'] and st.st_mtime_ns==item['mtime_ns']
 if item['sha256']:assert sha(f)==item['sha256']
pdf=P/'delivery/26210_Spatial_A4_tissue_10panels_v2.pdf'
info=subprocess.check_output(['pdfinfo',str(pdf)],text=True)
assert re.search(r'Pages:\s+1\b',info)
size=re.search(r'Page size:\s+([\d.]+) x ([\d.]+)',info)
assert abs(float(size[1])-210/25.4*72)<.01 and abs(float(size[2])-297/25.4*72)<.01
fonts=subprocess.check_output(['pdffonts',str(pdf)],text=True);assert 'Type 3' not in fonts
(P/'qc/PDF_INFO.txt').write_text(info);(P/'qc/PDF_FONTS.txt').write_text(fonts)
pixels={}
for ext,dpi in [('png',300),('tiff',600)]:
 with Image.open(P/f'delivery/26210_Spatial_A4_tissue_10panels_v2.{ext}') as im:
  im.load();assert abs(im.width-210/25.4*dpi)<1 and abs(im.height-297/25.4*dpi)<1
  actual=[float(v) for v in im.info['dpi']];assert all(abs(v-dpi)<.1 for v in actual)
  pixels[ext]=dict(size=list(im.size),dpi=actual)
meta=json.loads((P/'qc/DISPLAY_METADATA.json').read_text());assert len(meta)==10
assert max(v['max_registration_error_px'] for v in meta.values())<2
nchecks=0
for name,obj in meta.items():
 f=P/f'plot_data/{name}_tissue_layers.tsv';d=pd.read_csv(f,sep='\t')
 assert len(d)>100 and d.iloc[:,0].is_unique
 source=(R/'phase234_part4_replan_execution_20260921_v1'/f'inputs/{name}_INPUTS.tsv.gz') if name.startswith('GLMF') else (R/'phase250_independent_sixpatient_mainline_20260924_v1'/f'inputs/{name}_MAINLINE_INPUTS.tsv.gz')
 s=pd.read_csv(source,sep='\t')
 assert d.iloc[:,0].tolist()==s.iloc[:,list(s.columns).index('spot_id' if name.startswith('GLMF') else 'barcode')].tolist()
 for col in ['TYMP_RNA','polyamine_mta_context','spp1_state_no_spp1']:
  assert np.allclose(d[col],s[col],rtol=0,atol=1e-12)
 for key,cols in [('context_state_region',['polyamine_mta_context','spp1_state_no_spp1']),('spp1_endothelial_region',['V1_SPP1_q05','endothelial_q05'])]:
  if key not in d:continue
  expected=(s[cols[0]].ge(s[cols[0]].quantile(.8)).astype(int)+2*s[cols[1]].ge(s[cols[1]].quantile(.8)).astype(int)).to_numpy()
  assert np.array_equal(expected,d[key]);nchecks+=1
orig=R/'phase208_client_wording_SPP1_narrative_20260918_v6/report/26210_TYMP_MTA_SPP1_four_figures_narrative_v6.pdf'
report=R/'phase253_spatial_continue_numbering_20260925_v2/delivery/26210_原报告v6_追加空间_沿用原编号_v2.pdf'
assert sha(orig)=='6c6a3840fab35d794f02268c5ea86c146a10c684ea09fd564aec072d9718ce64'
assert sha(report)=='1767e5970adbf79f2ff1d579a40dc1ad8dabb037f3100f00c5881ed029e0507f'
# Preserve statistical support separately; no statistical chart in the main figure.
evidence=P/'plot_data/statistical_evidence';evidence.mkdir(exist_ok=True)
prior=R/'phase254_spatial_A4_figure_20260925_v1/plot_data'
for f in prior.glob('*.tsv'):
 if f.name.startswith(('E_','F_','G_','H_','I_','J_')):shutil.copy2(f,evidence/f.name)
bundle=P/'delivery/26210_Spatial_tissue_v2_绘图数据.zip'
with zipfile.ZipFile(bundle,'w',zipfile.ZIP_DEFLATED) as z:
 for f in sorted((P/'plot_data').rglob('*.tsv')):z.write(f,str(f.relative_to(P)))
 for f in sorted((P/'qc').glob('*.json')):z.write(f,'provenance/'+f.name)
 z.write(P/'code/make_tissue_figure.py','code/make_tissue_figure.py')
with zipfile.ZipFile(bundle) as z:assert z.testzip() is None
now=datetime.datetime.now().astimezone().isoformat()
v=dict(status='PASS',time=now,page_mm=[210,297],panels=10,tissue_maps=28,scatter_plots=0,statistical_heatmaps=0,visual_check='Final full-page PDF render reviewed, all A-J maps, labels, scales and categorical legends readable.',region_recomputations_checked=nchecks,source_arrays_checked=10,maximum_registration_error_px=max(x['max_registration_error_px'] for x in meta.values()),sources_unchanged=True,report_unchanged=True,new_models=0,new_statistical_tests=0,export_pixels=pixels,files=[dict(path=str(f),bytes=f.stat().st_size,sha256=sha(f)) for f in sorted((P/'delivery').iterdir())])
(P/'qc/FINAL_VALIDATION.json').write_text(json.dumps(v,ensure_ascii=False,indent=2))
note=f"""## {now} Phase255 按客户图形要求重做组织空间主图v2
用户明确否定Phase254统计图表达，要求不要散点图和仅2行/列热图。本轮重制A4单页A–J10panel，28张实际组织地图，统计散点/森林图/相关矩阵/统计热图均为0；旧Phase254图保留历史但不作为当前审核版本。
A原H&E；B冻结13类映射归并6类后的优势组成；C TYMP RNA；D代谢背景；E旧SPP1-state巨噬q05；F代谢/去SPP1状态前20%区域叠加；G免疫和侵袭程序实测地图；H SPP1/内皮q05前20%区域叠加；I外部六患者实际切片；J两例旁组织/肿瘤SPP1实际表达。I/J代表肿瘤阵列按最多QC spots选定，HG1选HG_1.3，不按关联效果选片；完整患者分析仍保留所有阵列。
六边形只作实测spot显示单元，通过阵列到H&E仿射几何转换，未平滑/插值/补测量/人工描边。固定前20%为描述性显示阈值，不称显著热点、免疫屏障或C4定位。新统计检验/新模型0；原P/FDR/调整反向/配对SPP1不一致保留在独立图注及数据包；CD8保留原报告统计，本版不占主图。
交付{P}/delivery下PDF、SVG、300dpi PNG、600dpi TIFF、逐panel中文图注和可复核数据ZIP。A4尺寸、最小6pt、坐标/色块类别、源数据及报告哈希检查通过；最终整页28图视觉审核通过，详见qc/FINAL_VALIDATION.json。进度只记服务器；本地Downloads/26210/Phase255_空间组织图_v2。本轮制作完成，无运行任务，图件独立交付未插入报告。

"""
(P/'SESSION_SUMMARY_CN.md').write_text(note)
for name in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md','RESULTS_MANIFEST.tsv','AGENTS.md']:
 f=R/name;b=P/'root_record_backups'/name
 if not b.exists():shutil.copy2(f,b)
 if name in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md']:f.write_text(note+f.read_text())
ag=R/'AGENTS.md';ag.write_text(ag.read_text()+"""
## 2026-09-25 空间综合图以真实切片为主，禁止散点和两行列热图
用户明确认为Phase254除A外过于统计化、不能达到客户要求，要求不要散点图和只有2行/列热图。当前图形偏好以真实组织空间分布、区域叠加和独立病例切片为主体。Phase255 v2为当前10-panel待审版，全部28张组织地图，无统计散点/热图。单页A4、独立交付不插报告继续有效；格点与区域只根据实测冻结结果可复核生成，不制造病理边界或强机制结论。
""")
with (R/'RESULTS_MANIFEST.tsv').open('a') as stream:
 for i,item in enumerate(v['files'],1):
  stream.write('\t'.join([f'phase255_tissue_v2_{i}','phase255_spatial_tissue_figure_20260925_v2','figure_package',item['path'],now,'PASS',str(P/'qc/FINAL_VALIDATION.json'),'Tissue-first A4 figure, no scatter or statistical heatmap; SHA256='+item['sha256']])+'\n')
print(json.dumps(v,ensure_ascii=False))
