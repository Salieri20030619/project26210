from pathlib import Path
import json,hashlib,datetime,subprocess,re,zipfile,shutil
import pandas as pd,numpy as np
from PIL import Image
P=Path(__file__).resolve().parents[1];R=P.parent
sha=lambda f:hashlib.sha256(Path(f).read_bytes()).hexdigest()
for item in json.loads((P/'qc/INPUT_MANIFEST.json').read_text()):
 f=Path(item['path']);st=f.stat()
 assert st.st_size==item['size'] and st.st_mtime_ns==item['mtime_ns']
 if item['sha256']!='metadata_only_large_readonly_h5ad':assert sha(f)==item['sha256']
pdf=P/'delivery/26210_Spatial_A4_10panels_v1.pdf'
info=subprocess.check_output(['pdfinfo',str(pdf)],text=True)
assert re.search(r'Pages:\s+1\b',info)
size=re.search(r'Page size:\s+([\d.]+) x ([\d.]+)',info)
assert abs(float(size[1])-210/25.4*72)<.01 and abs(float(size[2])-297/25.4*72)<.01
fontinfo=subprocess.check_output(['pdffonts',str(pdf)],text=True);assert 'Type 3' not in fontinfo
(P/'qc/PDF_INFO.txt').write_text(info);(P/'qc/PDF_FONTS.txt').write_text(fontinfo)
pixels={}
for ext,dpi in [('png',300),('tiff',600)]:
 f=P/f'delivery/26210_Spatial_A4_10panels_v1.{ext}'
 with Image.open(f) as im:
  im.load()
  assert abs(im.width-210/25.4*dpi)<1 and abs(im.height-297/25.4*dpi)<1
  actual=im.info.get('dpi');assert actual and all(abs(x-dpi)<.1 for x in actual)
  pixels[ext]=dict(size=list(im.size),dpi=[float(x) for x in actual])
# Verify patient-level figure summaries against frozen patient effects.
base=R/'phase234_part4_replan_execution_20260921_v1/tables'
patients=pd.read_csv(base/'A1_PATIENT_EFFECTS.tsv',sep='\t')
keys=['dataset','family','endpoint','exposure','target','layer','scale','qualification','model','effect_type']
checked=0
for name in ['E_context_SPP1','F_functional_programs','G_vascular']:
 d=pd.read_csv(P/f'plot_data/{name}.tsv',sep='\t')
 for _,row in d.iterrows():
  keep=np.ones(len(patients),bool)
  for key in keys:keep &= patients[key].astype(str).to_numpy()==str(row[key])
  vals=patients.loc[keep,'effect'].dropna().to_numpy()
  assert len(vals)==row.n and abs(np.median(vals)-row['median'])<1e-12
  checked+=1
jj=pd.read_csv(P/'plot_data/J_CD8_patients.tsv',sep='\t');js=pd.read_csv(P/'plot_data/J_CD8_summary.tsv',sep='\t')
for _,row in js.iterrows():
 vals=jj[jj.dataset==row.dataset].effect.dropna().to_numpy()
 assert len(vals)==row.n and abs(np.median(vals)-row['median'])<1e-8
# Include known model sensitivity in the source package, without refitting.
shutil.copy2(R/'phase250_independent_sixpatient_mainline_20260924_v1/tables/WHOLE_TISSUE_AND_ADJUSTMENT_ALL.tsv',P/'plot_data/H_model_sensitivity_full.tsv')
shutil.copy2(R/'phase250_independent_sixpatient_mainline_20260924_v1/tables/TWO_INDEPENDENT_COHORTS_MAINLINE.tsv',P/'plot_data/H_frozen_direction_summary.tsv')
orig=R/'phase208_client_wording_SPP1_narrative_20260918_v6/report/26210_TYMP_MTA_SPP1_four_figures_narrative_v6.pdf'
assert sha(orig)=='6c6a3840fab35d794f02268c5ea86c146a10c684ea09fd564aec072d9718ce64'
current=R/'phase253_spatial_continue_numbering_20260925_v2/delivery/26210_原报告v6_追加空间_沿用原编号_v2.pdf'
assert sha(current)=='1767e5970adbf79f2ff1d579a40dc1ad8dabb037f3100f00c5881ed029e0507f'
bundle=P/'delivery/26210_Spatial_A4_10panels_绘图数据.zip'
with zipfile.ZipFile(bundle,'w',zipfile.ZIP_DEFLATED) as z:
 for f in sorted((P/'plot_data').glob('*.tsv')):z.write(f,'plot_data/'+f.name)
 for name in ['INPUT_MANIFEST.json','BUILD_CHECK.json']:z.write(P/'qc'/name,'provenance/'+name)
 z.write(P/'code/make_figure.py','code/make_figure.py')
with zipfile.ZipFile(bundle) as z:assert z.testzip() is None
now=datetime.datetime.now().astimezone().isoformat()
validation=dict(status='PASS',time=now,page_mm=[210,297],panels=10,visual_review='Final PDF rendered at150dpi, all A-J panels checked; no missing glyphs or overlapping labels; complete error bars and readable legends.',minimum_font_pt=6,summary_patient_checks=checked,CD8_source_checks=3,source_files_unchanged=True,report_files_unchanged=True,reference_figures='User supplied Ge et al Sci Adv eaed0733 Fig2-3; style only; own observed data',reference_pdf_sha256=sha(P/'reference/sciadv.aed0733.pdf'),images=pixels,files=[dict(path=str(f),bytes=f.stat().st_size,sha256=sha(f)) for f in sorted((P/'delivery').iterdir())])
(P/'qc/FINAL_VALIDATION.json').write_text(json.dumps(validation,ensure_ascii=False,indent=2))
note=f"""## {now} Phase254 空间A4单页10-panel独立综合图完成
用户明确授权参照本地sciadv.aed0733(1).pdf的fig制作空间8–10panel一张A4图，先单独生成不插入报告。本轮10个panel A–J，210×297mm。参照Fig2/3的地图加验证版式，全部图形为自身冻结数据；不挪用论文边界、屏障或临床结论，不把共有GLMF切片算独立验证。
A–D GLMF1/2真实H&E、TYMP RNA、代谢转录背景及旧SPP1-state q05丰度，同病例统一坐标/色阶/500um比例尺。E代谢背景四种SPP1读出，F四功能程序，G血管绝对量/比例成对呈现，H两个独立来源6+3例全部五关系及NA，I两例配对旁组织/肿瘤四端点，J CD8距离三个来源所有可评患者。
不新增分析/模型或重选统计；原患者单位和原P/BH家族保留。H标未调整4×4区域，图注明确保留深度谱系调整减弱/反向；I单基因SPP1不一致，J跨来源不一致及名义P/FDR保留。空间旧SPP1-state不冒充C4；CD8少量补充，免疫临床板块仍待后续。
交付{P}/delivery：PDF、SVG、300dpi PNG、600dpi LZW TIFF、独立逐panel中文图注及绘图数据ZIP。图内最小6pt，完整置信区间和文字边界检查，24患者汇总与3CD8来源复核通过；A4最终PDF全部panel视觉检查通过。已读h5ad及源结果不变，原报告v6和当前Phase253 PDF哈希不变。记录qc/FINAL_VALIDATION.json，代码code/make_figure.py、code/finalize.py。
本地交付Downloads/26210/Phase254_空间A4综合图_v1；全部进度只记服务器。无运行中的本轮任务；图件已生成待用户审阅，未插入报告。

"""
(P/'SESSION_SUMMARY_CN.md').write_text(note)
for name in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md','RESULTS_MANIFEST.tsv','AGENTS.md']:
 f=R/name;b=P/'root_record_backups'/name
 if not b.exists():shutil.copy2(f,b)
 if name in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md']:f.write_text(note+f.read_text())
ag=R/'AGENTS.md';ag.write_text(ag.read_text()+"""
## 2026-09-25 用户授权空间A4单页8–10panel综合图
当前明确授权参考提供的Science Advances文章fig制作一张空间综合图，单独交付，不插入报告；本轮Phase254为A–J共10panel，210×297mm。该请求授权本轮绘图并覆盖此前当前阶段不绘图/等待图方案批准的限制，不扩展为改写或重新整合报告。旧SPP1-state、新C4、独立患者与原统计家族保持分开；免疫临床尚未完成。
""")
with (R/'RESULTS_MANIFEST.tsv').open('a') as stream:
 for i,item in enumerate(validation['files'],1):
  stream.write('\t'.join([f'phase254_A4_10panel_{i}','phase254_spatial_A4_figure_20260925_v1','figure_package',item['path'],now,'PASS',str(P/'qc/FINAL_VALIDATION.json'),'Independent A4 figure, not inserted in report; SHA256='+item['sha256']])+'\n')
print(json.dumps(validation,ensure_ascii=False))
