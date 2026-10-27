from pathlib import Path
import json,hashlib,datetime,shutil
w=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase128_section2_5_metaviper_redraw_remote_v1');root=w.parent;q=w/'qc'
assert json.loads((q/'INDEPENDENT_VALIDATION.json').read_text())['status']=='PASS'
assert json.loads((q/'VISUAL_REVIEW.json').read_text())['status']=='PASS'
assert json.loads((q/'DELIVERY_VALIDATION.json').read_text())['sha256_matches']==21
rows=json.loads((q/'DELIVERY_MANIFEST.json').read_text())
for row in rows:assert hashlib.sha256((w/'figures'/row['name']).read_bytes()).hexdigest()==row['sha256']
backup=w/'state_backup_before_completion';backup.mkdir(exist_ok=True)
for name in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md','RESULTS_MANIFEST.tsv']:
 if not (backup/name).exists():shutil.copy2(root/name,backup/name)
note="## Phase128 section2.5 metaVIPER figures delivered (2026-09-09)\n\n- 用户要求“现在开始做2.5”。按原方案及image13–19参考完成可评估部分与明确标注的补充图，交付至 C:/Users/31937/Downloads/26210/06_2.5_metaVIPER_v1：10 PNG+10 PDF+中文说明，21/21 SHA256一致。\n- 唯一工作目录：`phase128_section2_5_metaviper_redraw_remote_v1`。图件为MR效应散点/患者敏感性、F1主导小提琴/配对中位数、3个F1终点、原始NES UMAP、三方法轨迹曲线、三方法拟合UMAP、末端变化/共同正向比例、TYMP/RRM2B RNA敏感性。\n- 从正式5,566×11,921 RDS直接提取4个MR蛋白与背景HIF1A；与历史导出逐值完全一致。20位患者，原11,950 UMAP保持不变，29个无评分细胞仅灰色背景，不补值。metaVIPER、cNMF、轨迹模型均未重跑。\n- 8个最终MR蛋白只有CASP9/CHMP1A/PRSS8/PXN可推断；KDELC2/POGLUT3、RRM2B、TMEM25、TYMP无合格肾癌网络。TYMP/RRM2B字面蛋白活性端点仍不可评估，RNA补充不替代该端点；历史替代路线的证据类型未改变。\n- 四蛋白MR logOR与患者等权平均活性Pearson r=-0.648335，仅描述，不画四点密度或推断性P。固定四蛋白的患者敏感性20/20负向，中位r=-0.663739，P=1.90735e-6，FDR=3.81470e-6；20患者不增加MR蛋白数。\n- CASP9/CHMP1A/PRSS8与F1连续及患者内四分位终点为正向；PXN负向或未获支持。主F1分组按cNMF主导模块，不按全局中位数：9患者、7,021细胞，F1主导4,069/其他2,952，每臂>=10细胞，按探索性效应与一致性呈现，完整P/FDR保留。\n- 三种主拟时序（Monocle3/Slingshot/scTour）不支持稳定的跨方法末端共激活；Slingshot特异方向、PXN仅Monocle3的正向末端结果均保留。CytoTRACE2/CellRank保留原敏感性类型，不作末端拟时序替代。\n- 新增显示拟合仅用原患者20秩箱均值及箱细胞数，Gaussian带宽固定0.075秩单位，患者内拟合后跨患者等权。带为患者IQR而非CI；拟合UMAP与原始NES明确分开。同蛋白三方法曲线共享纵轴；UMAP均以NES=0居中，拟合图跨方法共享该蛋白完整上下限。\n- HIF1A可在独立metaVIPER肾癌网络中估计且F1关联负向；这不同于2.4 pySCENIC调控子缺失。作为背景目标独立统计，不计入四个MR靶点。\n- 复核300条F1效应、700条轨迹相关、420条末端变化及原分箱均值（最大误差7.99e-15）；所有原BH家族复算通过。10/10最终PDF单页/视觉核查通过；标题、RNA脚注和拟合色条刻度均修正，初稿保留服务器。\n- 当前断点：用户审阅2.5。图件工作完成；缺失的字面TYMP/RRM2B蛋白活性端点仍按不可评估记录。脚本、完整统计/矩阵/拟合值、源SHA256、软件版本及进度仅保存在服务器；无活动绘图任务。\n"
for name in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md']:
 p=root/name;old=p.read_text()
 if '## Phase128 section2.5' not in old:p.write_text(note+'\n\n'+old,encoding='utf8')
mf=root/'RESULTS_MANIFEST.tsv';old=mf.read_text()
with mf.open('a') as f:
 for i,row in enumerate(rows):
  key='phase128_metaviper_figure_'+str(i+1)
  if key+'\t' not in old:f.write(key+'\tphase128_section25_metaviper_redraw\treview_deliverable\t'+str(w/'figures'/row['name'])+'\t'+datetime.datetime.now().isoformat()+'\tPASS_REVIEW_CHECKS\tSHA256 '+row['sha256']+'\t21/21 delivered; unavailable activity endpoints retained; fitted and RNA views separated\n')
manifest=json.loads((q/'SOURCE_MANIFEST.json').read_text());known={r['path'] for r in manifest}
for p in (w/'scripts').iterdir():
 if p.is_file() and str(p) not in known:manifest.append({'path':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
(q/'SOURCE_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
state={'status':'PASS_FIGURES_DELIVERED_WITH_SCOPE_LIMITATION','section':'2.5','figure_families':10,'delivered_files':21,'sha256_matches':21,'PDF_checks':10,'independent_validation':'PASS','unavailable_literal_protein_activity':['TYMP','RRM2B','KDELC2/POGLUT3','TMEM25'],'RNA_not_protein_activity':True,'active_jobs':False,'next':'User reviews section2.5','local_delivery':'C:/Users/31937/Downloads/26210/06_2.5_metaVIPER_v1'}
for name in ['COMPLETION.json','CURRENT_TASK.json']:(q/name).write_text(json.dumps(state,indent=2))
print(json.dumps(state))
