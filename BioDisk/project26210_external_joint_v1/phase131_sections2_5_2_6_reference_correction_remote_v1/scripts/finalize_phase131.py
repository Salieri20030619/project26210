from pathlib import Path
import json,hashlib,shutil,datetime,sys
import numpy,pandas,scipy,matplotlib
w=Path(__file__).resolve().parents[1];root=w.parent;q=w/'qc'
for n in ['SECTION25_VALIDATION.json','SECTION26_VALIDATION.json','VISUAL_REVIEW.json']:
 assert json.loads((q/n).read_text())['status']=='PASS'
receipt=json.loads((q/'DELIVERY_VALIDATION.json').read_text());assert receipt['sha256_matches']==20
rows=json.loads((q/'DELIVERY_MANIFEST.json').read_text())
for r in rows:assert hashlib.sha256((w/('figures_'+r['section'])/r['name']).read_bytes()).hexdigest()==r['sha256']
backup=w/'state_backup_before_completion';backup.mkdir(exist_ok=True)
for n in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv']:
 if not (backup/n).exists():shutil.copy2(root/n,backup/n)
note='''## Phase131 sections2.5/2.6 reference-style corrections delivered (2026-09-09)

- 用户指出2.5未体现示例“一圈一圈”的相关图、2.6图件与示例不符。本次重新解析原DOCX图文位置：2.5第1项对应image13/14二维KDE密度等高线；2.6对应image15双分面气泡矩阵、image16横排四联UMAP。image15/16此前被混列为2.5参考，现纠正章节对应。
- 用户进一步明确两种关系都画并区分：MR效应×活性，以及F1×各蛋白活性。两种均已完成独立图件，没有再次用普通散点替代等高线。
- 唯一工作目录：`phase131_sections2_5_2_6_reference_correction_remote_v1`。2.5交付 `C:/Users/31937/Downloads/26210/06_2.5_density_reference_v2`（2 PNG+2 PDF+说明=5文件）；2.6交付 `C:/Users/31937/Downloads/26210/07_2.6_reference_layout_v2`（7 PNG+7 PDF+说明=15文件）。全部20/20 SHA256一致。
- 2.5 MR图仍为CASP9/CHMP1A/PRSS8/PXN四个实际蛋白点，横轴原log OR、纵轴20患者等权平均metaVIPER NES，Pearson r=-0.64833498。Gaussian KDE为4点描述性平滑，实际点和蛋白名全部保留，无4蛋白推断性P；Scott带宽，0.7/1/1.4倍敏感性保存在QC。现图替换旧普通散点作为该参考要求的主图，不增加有效MR蛋白数。
- 2.5 F1图以四蛋白分别绘制密度等高线，11,921细胞、20患者；每患者KDE总权重相等。Weighted cell Pearson r只描述分布，右侧保留原患者内Spearman中位rho/P/FDR（四MR蛋白家族）；复算一致。F1观测支持0–1以虚线明确标记，区间外仅Gaussian平滑尾部，无新增细胞或分数。PXN细胞描述r为正、患者内中位rho为负，两项含义分别标注。缺失TYMP/RRM2B等活性仍不补值。
- 2.6七图：四AUCell横排magma UMAP；四通路TYMP-/+双分面气泡；85通路同式完整气泡；F1+四通路圆点Corrplot；12酶竖排F1四分位气泡；四scFEA代理横排UMAP；四通路两种配对比较的显著性气泡（面积=-log10患者P，颜色=中位delta，黑轮廓=FDR<0.05）。逐图参考对应写入中文说明及registry。
- 双分面描述图按患者×cNMF程序配对，两臂各>=10细胞；F1–F7患者数6/9/0/2/3/2/2，F3以灰叉显示无合格配对。患者等权平均AUCell及AUC>0比例分别映射颜色/面积，不将细胞P伪装为参考P。n<10仅作探索性描述；四重点通路几乎均正分导致面积接近，明确记录真实原因。另以原14患者整体配对统计气泡表达显著性。
- 新圆点Corrplot使用患者内rho中位数，F1四条原效应与原P/FDR保持；通路之间矩阵项仅描述，无新增检验。12酶沿用17位合格患者的原四分位表，转置为基因行/分组列，不按基因缩放颜色。
- 正式细胞、评分、UMAP、模型、主统计结论均保持：三AUCell通路正向支持、甘油磷脂配对阴性；scFEA局部代理及MTA未显式建模范围保留。此次主要为图形匹配修正，新增汇总不取代原患者统计。
- 9/9 PDF单页/边界检查及逐图视觉核查通过；小幅下移短页面标题、扩大等高线视窗使曲线完整，旧初稿保留。原始输入与参考SHA256、密度网格、患者分组表、软件、脚本、复核和交付收据仅存服务器。
- 当前断点：用户审阅2.5密度相关图及2.6参考布局v2；无活动绘图任务。Phase130的2.4 B/E修正、Phase129的2.7完成状态保持，不自动推进下一章节。
'''
for n in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
 p=root/n;old=p.read_text()
 if '## Phase131 sections2.5/2.6' not in old:
  tmp=p.with_name(p.name+'.phase131.tmp');tmp.write_text(note+'\n\n'+old,encoding='utf8');tmp.replace(p)
p=root/'RESULTS_MANIFEST.tsv';old=p.read_text()
with p.open('a') as f:
 for i,r in enumerate(rows):
  key='phase131_reference_correction_'+str(i+1)
  if key+'\t' not in old:f.write(key+'\tphase131_sections25_26_reference_correction\treview_deliverable\t'+str(w/('figures_'+r['section'])/r['name'])+'\t'+datetime.datetime.now().isoformat()+'\tPASS_REVIEW_CHECKS\tSHA256 '+r['sha256']+'\t20/20 delivered; MR and F1 KDE separated; image15/16 section2.6 layout matched; no new inferential tests\n')
manifest=json.loads((q/'SOURCE_MANIFEST.json').read_text());lookup={r['path']:r for r in manifest}
for p in (w/'scripts').glob('*'):
 if p.is_file():lookup[str(p)]={'path':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
(q/'SOURCE_MANIFEST.json').write_text(json.dumps(list(lookup.values()),indent=2))
(q/'SOFTWARE_VERSIONS.json').write_text(json.dumps({'python':sys.version,'executable':sys.executable,'numpy':numpy.__version__,'pandas':pandas.__version__,'scipy':scipy.__version__,'matplotlib':matplotlib.__version__},indent=2))
state={'phase':131,'status':'PASS_REFERENCE_CORRECTIONS_DELIVERED','sections':['2.5','2.6'],'figures':9,'section_figures':{'2.5':2,'2.6':7},'delivered_files':20,'sha256_matches':20,'PDF_checks':9,'visual_review':'PASS','MR_and_F1_density_figures_separated':True,'scientific_models_rerun':False,'active_jobs':False,'next':'User reviews section2.5 KDE and section2.6 reference-layout v2','local_destinations':receipt['destinations']}
for n in ['COMPLETION.json','CURRENT_TASK.json']:(q/n).write_text(json.dumps(state,indent=2))
print(json.dumps(state),flush=True)
