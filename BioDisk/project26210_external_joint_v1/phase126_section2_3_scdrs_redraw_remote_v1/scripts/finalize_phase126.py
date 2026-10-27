from pathlib import Path
import json,hashlib,datetime,shutil
w=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase126_section2_3_scdrs_redraw_remote_v1');root=w.parent;q=w/'qc'
assert json.loads((q/'INDEPENDENT_VALIDATION.json').read_text())['status']=='PASS'
assert json.loads((q/'VISUAL_REVIEW.json').read_text())['status']=='PASS'
assert len(json.loads((q/'PDF_VALIDATION.json').read_text()))==4
rows=json.loads((q/'DELIVERY_MANIFEST.json').read_text())
for row in rows:
 f=w/'figures'/row['name'];assert hashlib.sha256(f.read_bytes()).hexdigest()==row['sha256']
backup=w/'state_backup_before_completion';backup.mkdir(exist_ok=True)
for name in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md','RESULTS_MANIFEST.tsv']:
 if not (backup/name).exists():shutil.copy2(root/name,backup/name)
note="## Phase126 section2.3 scDRS redraw completed; local export pending (2026-09-09)\n\n- 用户最新指令：“进入2.3绘图”。已按原方案2.3及image11参考完成绘图。\n- 唯一工作目录：`phase126_section2_3_scdrs_redraw_remote_v1`。正式图件4 PNG + 4 PDF + 中文说明：三面板总图、A恶性上皮scDRS UMAP、B F1-scDRS密度着色散点、C患者层效应图。\n- 复用Phase25 compute_v2/analysis_v1正式结果与冻结F1/UMAP；11,950细胞、26患者、5队列。原20基因GWAS位点集过滤后17基因。无重跑scDRS/cNMF或细胞标签变更。\n- 当前数据不支持稳定F1-scDRS关联：20合格患者8正12负，中位rho=-0.0490078，患者bootstrap 95%区间[-0.1291551,0.0635501]，双侧符号翻转P=0.4833。单个预设主要终点无需额外FDR。细胞Pearson r=0.009127仅描述，不计算细胞P；scDRS细胞描述性BH FDR<0.10/0.20均0个。\n- 全细胞身份、F1、scDRS评分、坐标核对通过；患者效应以rank-Pearson独立复算，bootstrap/P复算通过。4/4 PDF单页且渲染图视觉核查通过。完整表、源哈希、脚本、软件版本均在服务器。\n- 交付清单：qc/DELIVERY_MANIFEST.json，9个文件，约3.9MB。拟交付路径：C:/Users/31937/Downloads/26210/04_2.3_scDRS_v1。\n- 本地正式导出尚未执行：自动审批两次拒绝该批下载，理由为项目派生图件属敏感数据，可信用户消息未明确授权本批载荷及目的地，服务器交接记录不能替代用户授权。已向用户发出明确授权问题。不得绕过；得到授权后可重试原交付操作并核对9/9 SHA256。\n- 当前断点：图件完成；等待用户选择“授权下载此批图件”或“仅保留服务器版本”，随后审阅2.3。无活动绘图任务。\n"
for name in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md']:
 f=root/name;old=f.read_text()
 if '## Phase126 section2.3' not in old:f.write_text(note+'\n\n'+old,encoding='utf8')
mf=root/'RESULTS_MANIFEST.tsv';old=mf.read_text()
with mf.open('a') as out:
 for i,row in enumerate(rows):
  key='phase126_scdrs_figure_'+str(i+1)
  if key+'\t' not in old:out.write(key+'\tphase126_section23_scdrs_redraw\tserver_review_artifact\t'+str(w/'figures'/row['name'])+'\t'+datetime.datetime.now().isoformat()+'\tPASS_REVIEW_CHECKS\tSHA256 '+row['sha256']+'\tFrozen scDRS/F1/UMAP; local export pending explicit authorization\n')
state={'status':'FIGURES_COMPLETE_EXPORT_PENDING_AUTHORIZATION','figures':4,'formal_files':9,'independent_validation':'PASS','pdf_visual_checks':4,'local_export_completed':False,'active_jobs':False,'next':'User authorizes local export or chooses server-only; review section2.3'}
(q/'COMPLETION.json').write_text(json.dumps(state,indent=2))
(q/'CURRENT_TASK.json').write_text(json.dumps(state,indent=2))
print(json.dumps(state))
