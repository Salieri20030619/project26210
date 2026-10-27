from pathlib import Path
import json,hashlib,datetime,shutil
w=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase127_section2_4_pyscenic_redraw_remote_v1');root=w.parent;q=w/'qc'
assert json.loads((q/'INDEPENDENT_VALIDATION.json').read_text())['status']=='PASS'
assert json.loads((q/'VISUAL_REVIEW.json').read_text())['status']=='PASS'
assert json.loads((q/'DELIVERY_VALIDATION.json').read_text())['sha256_matches']==15
rows=json.loads((q/'DELIVERY_MANIFEST.json').read_text())
for row in rows:assert hashlib.sha256((w/'figures'/row['name']).read_bytes()).hexdigest()==row['sha256']
backup=w/'state_backup_before_completion';backup.mkdir(exist_ok=True)
for name in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md','RESULTS_MANIFEST.tsv','AGENTS.md']:
 if not (backup/name).exists():shutil.copy2(root/name,backup/name)
note="## Phase127 section2.4 pySCENIC redraw completed (2026-09-09)\n\n- 用户已明确授权“所有已知数据库与本地的互相传输”，并要求进入2.4。本次7 PNG+7 PDF+中文说明已交付 C:/Users/31937/Downloads/26210/05_2.4_pySCENIC_v1，15/15 SHA256一致；此前2.3也已补交9/9一致。\n- 唯一工作目录：`phase127_section2_4_pyscenic_redraw_remote_v1`。按原方案image12形式提供两算法F1–F7 RSS排名、cNMF患者平衡热图、候选连续关联、探索性分组效应、环形网络及AUCell UMAP。\n- 保留原11,950细胞、26患者、5队列；GENIE3 299、GRNBoost2 314调控子，共有197。模型、AUCell、cNMF和UMAP未重跑；仅新增RSS与绘图描述性汇总。\n- ILF2/HOXD3/HNF4A在20位患者的连续F1关联均双算法FDR<0.05；F1主导组比较n=9，以方向和效应量呈现。ILF2在GRNBoost2、HOXD3在GENIE3的分组结果仅名义显著；完整P/FDR保留，不将历史3/4分级改写为因果驱动。无4/4高置信候选。\n- HIF1A在两矩阵缺失，Phase28 8/8靶向敏感性仍未恢复；不画虚构的HIF1A活性。未检出不等于生物学无作用。\n- 热图24调控子，21患者提供37份患者×模块谱；每谱>=10细胞，患者内中位数再跨患者中位数，两算法分别行z标准化。F1–F7患者n分别9/12/2/4/4/2/4，不当作独立重复相加；不新增检验。\n- 环形图按两种算法真实边及方法内权重排序，支持类别各最多10靶点，排除自身。ILF2/HOXD3/HNF4A各显示14/24/30个靶点，共86边；完整非自身靶点数分别GENIE3 4/14/87、GRNBoost2 1133/53/518，共同2/6/39。ILF2→TYMP仅GRNBoost2预测；其他两TF均无TYMP边，不建立直接调控共识。\n- 校验：输入SHA256一致；156候选患者记录资格/效应复算，原家族BH复算，RSS独立误差2.22e-16；7/7 PDF最终结构及渲染视觉核查通过。D轴标题和E长基因名间距已修正。\n- 脚本、完整患者表、逐细胞矩阵、网络及进度仅保存在服务器。当前断点：用户审阅2.4；无活动绘图任务。\n"
for name in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md']:
 p=root/name;old=p.read_text()
 if '## Phase127 section2.4' not in old:p.write_text(note+'\n\n'+old,encoding='utf8')
p=root/'AGENTS.md';old=p.read_text()
auth="\n## 2026-09-09 用户明确传输授权\n\n用户指令：“授权所有已知数据库与本地的互相传输，现在进入2.4”。此授权持续适用于已知项目数据源、服务器与本地之间的传输；无需重复请求已授权的传输。分析和进度仍只在服务器记录，图件交付Downloads/26210。当前阶段以CODEX_HANDOFF.md顶部最新记录为准。\n"
if '## 2026-09-09 用户明确传输授权' not in old:p.write_text(old+auth,encoding='utf8')
mf=root/'RESULTS_MANIFEST.tsv';old=mf.read_text()
with mf.open('a') as f:
 for i,row in enumerate(rows):
  key='phase127_pyscenic_figure_'+str(i+1)
  if key+'\t' not in old:f.write(key+'\tphase127_section24_pyscenic_redraw\treview_deliverable\t'+str(w/'figures'/row['name'])+'\t'+datetime.datetime.now().isoformat()+'\tPASS_REVIEW_CHECKS\tSHA256 '+row['sha256']+'\tPatient-unit evidence; n9 group exploratory; 15/15 delivered\n')
state={'status':'PASS_COMPLETE_PENDING_USER_REVIEW','section':'2.4','figure_families':7,'delivered_files':15,'delivery_sha256_matches':15,'independent_validation':'PASS','pdf_final_checks':7,'active_jobs':False,'next':'User reviews section2.4','local_delivery':'C:/Users/31937/Downloads/26210/05_2.4_pySCENIC_v1'}
for name in ['COMPLETION.json','CURRENT_TASK.json']:(q/name).write_text(json.dumps(state,indent=2))
print(json.dumps(state))
