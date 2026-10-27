from pathlib import Path
import json,hashlib,datetime,shutil
w=Path(__file__).resolve().parents[1];root=w.parent;q=w/'qc'
validation=json.loads((q/'INDEPENDENT_VALIDATION.json').read_text())
assert validation['status']=='PASS'
assert json.loads((q/'VISUAL_REVIEW.json').read_text())['status']=='PASS'
receipt=json.loads((q/'DELIVERY_VALIDATION.json').read_text());assert receipt['status']=='PASS' and receipt['sha256_matches']==40
pdf=json.loads((q/'PDF_VALIDATION.json').read_text());assert len(pdf)==19 and all(not r['outside_page_text'] and r['pages']==1 for r in pdf)
rows=json.loads((q/'DELIVERY_MANIFEST.json').read_text())
for r in rows:assert hashlib.sha256((w/('figures_'+r['section'])/r['name']).read_bytes()).hexdigest()==r['sha256']
backup=w/'state_backup_before_completion';backup.mkdir(exist_ok=True)
for name in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md','RESULTS_MANIFEST.tsv']:
 if not (backup/name).exists():shutil.copy2(root/name,backup/name)
note='''## Phase129 sections2.6/2.7 figures delivered (2026-09-09)

- 用户明确要求“继续2.6，2.7绘图”。已完成2.6九套、2.7十套图；每套300 dpi PNG+PDF。2.6交付 `C:/Users/31937/Downloads/26210/07_2.6_metabolism_v1`（19文件），2.7交付 `C:/Users/31937/Downloads/26210/08_2.7_trajectory_v1`（21文件），共40/40 SHA256一致，含两份中文说明与完整汇总统计。
- 唯一工作目录：`phase129_sections2_6_2_7_redraw_remote_v1`。原方案定位为Phase24 checklist_structure_v1.txt P043–P061；参考image20起已属后续髓系章节，未误用为恶性上皮数据或参考分组。
- 2.6：四通路/F1/TYMP UMAP、TYMP及F1高克隆患者配对图、两预测变量相关热图、12指定酶/转运体F1秩四分位气泡图及患者相关图、scFEA四轴UMAP/配对/相关，以及全部85通路TYMP比较补充。11,950细胞、26患者、5队列、原UMAP/cNMF/CNV/代谢模型冻结。
- 14位配对患者的TYMP阳性比较支持嘧啶、Cys/Met、Arg/Pro三通路升高（中位delta=0.009482/0.010186/0.014007，正向13/14、14/14、13/14，预设四通路FDR均0.002262）；F1高主克隆比较也支持三通路升高。甘油磷脂配对差异未获支持（TYMP q=0.801732，克隆q=0.379471）。连续F1关联n=20，三通路q=0.019298/0.002482/0.002482，甘油磷脂q=0.489783。
- 新增85通路TYMP配对比较使用原14合格患者、每臂>=20细胞、患者组中位数，85项独立BH，49/85 q<0.05，完整正负与阴性均保留为探索性补充。不得与原四通路家族混用。12原方案基因F1关联独立BH，TYMP/TK1/ODC1/SRM/PCYT2通过；其余完整保留。每基因可估计n=13–20，非细胞层检验。
- 12基因直接从正式H5AD counts提取，全部保留基因UMI归一化；TYMP原始计数逐值一致，log1p CP10K核对通过。四分位图保留F1同分值，各患者四组均>=5细胞，合格17患者；颜色未按基因缩放。
- scFEA四轴与F1的患者相关均通过原四轴BH。独立重建模块-z聚合误差<5e-14。复核脚本显式将精确恒定模块置零以匹配原R safe_z，修复NumPy浮点标准差残差；原结果不变。M83恒定，胆碱局部代理变化来自M21；MTA不是显式模型输出，未验证匹配代谢组或44倍变化。
- 2.7：五方法总览（Monocle3真实principal graph、Slingshot、scTour、CytoTRACE2、CellRank）、三命运概率、真实转移向量分箱图、四AUCell/四scFEA轴三方法曲线、三组12酶RNA曲线、48项患者相关矩阵及三方法状态构成图。原root、三条Slingshot lineage及主lineage1、全部模型不重跑。
- 三主拟时序不支持稳定跨方法代谢升高；Slingshot嘧啶关联rho=0.116058，14/20正向，P=0.00365448，原12终点BH q=0.04385376，其他方法未重复。CytoTRACE2与F1 rho=-0.203566、q=0.048317，与甘油磷脂rho=-0.131633、q=0.011913；保留其潜能类型。CellRank为0.8 pseudotime+0.2 connectivity，三态概率和为1，非RNA velocity，TYMP-state终态不是独立因果验证。
- 动态显示：20位>=30细胞患者，各自20秩箱；Gaussian h=0.075按箱细胞数加权。曲线为患者拟合中位数、带为患者IQR（非CI），同特征跨方法同纵轴。构成图为患者等权平均且总和为1。RNA动态只作描述。原未平滑患者统计不由显示拟合替换。
- 复核：224条原代谢配对差值、320条代谢相关、956条轨迹患者相关及48汇总全部通过，原P/BH均复算；新增85+12项检验和FDR由R独立重现。19/19 PDF单页/边界检查及逐图渲染视觉审阅通过。所有输入（含1.82GB正式H5AD）SHA256、完整矩阵/患者表/拟合值/软件及脚本保存在服务器。
- 当前断点：用户审阅2.6、2.7；两节绘图和交付完成，无活动绘图任务。下一章节不自动启动。进度仅在服务器记录。
'''
for name in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md']:
 p=root/name;old=p.read_text()
 if '## Phase129 sections2.6/2.7' not in old:
  temp=p.with_name(p.name+'.phase129.tmp');temp.write_text(note+'\n\n'+old,encoding='utf8');temp.replace(p)
mf=root/'RESULTS_MANIFEST.tsv';old=mf.read_text()
with mf.open('a') as out:
 for i,r in enumerate(rows):
  key='phase129_section26_27_figure_'+str(i+1)
  if key+'\t' not in old:
   out.write(key+'\tphase129_sections26_27_redraw\treview_deliverable\t'+str(w/('figures_'+r['section'])/r['name'])+'\t'+datetime.datetime.now().isoformat()+'\tPASS_REVIEW_CHECKS\tSHA256 '+r['sha256']+'\t40/40 delivered; original models frozen; patient-level statistics; scope and negative evidence retained\n')
manifest=json.loads((q/'SOURCE_MANIFEST.json').read_text());lookup={r['path']:r for r in manifest}
for p in (w/'scripts').glob('*'):
 if p.is_file():lookup[str(p)]={'path':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
(q/'SOURCE_MANIFEST.json').write_text(json.dumps(list(lookup.values()),indent=2))
state={'phase':129,'sections':['2.6','2.7'],'status':'PASS_FIGURES_DELIVERED','figure_families':19,'section_figures':{'2.6':9,'2.7':10},'delivered_files':40,'sha256_matches':40,'PDF_checks':19,'visual_review':'PASS','independent_validation':'PASS','active_jobs':False,'next':'User reviews sections2.6 and2.7','local_destinations':receipt['destinations'],'scientific_scope':'RNA activity/flux proxies/cross-sectional ordering distinguished; no matched metabolomics or MTA explicit flux'}
for name in ['COMPLETION.json','CURRENT_TASK.json']:(q/name).write_text(json.dumps(state,indent=2))
print(json.dumps(state),flush=True)
