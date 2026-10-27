from pathlib import Path
import json,hashlib,shutil,datetime,sys
import numpy,pandas,scipy,matplotlib
w=Path(__file__).resolve().parents[1];root=w.parent;q=w/'qc'
assert (q/'R_INDEPENDENT_VALIDATION.txt').read_text().startswith('PASS')
assert json.loads((q/'PRIOR_AXIS_RECHECK.json').read_text())['status']=='PASS'
assert json.loads((q/'VISUAL_REVIEW.json').read_text())['status']=='PASS'
r=json.loads((q/'DELIVERY_VALIDATION.json').read_text());assert r['sha256_matches']==9
rows=json.loads((q/'DELIVERY_MANIFEST.json').read_text())
for row in rows:assert hashlib.sha256((w/'figures'/row['name']).read_bytes()).hexdigest()==row['sha256']
backup=w/'state_backup_before_completion';backup.mkdir(exist_ok=True)
for n in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv','AGENTS.md']:
 if not (backup/n).exists():shutil.copy2(root/n,backup/n)
note='''## Phase133 sections2.5/2.7 supportive supplementation complete (2026-09-09)

- 用户明确调整定位：湿实验已验证，生信只需方向相关或名义显著的对应线索，不要求独立强证据闭环。用户同意2.3舍弃，并要求2.5、2.7补充。2.3已从正文主图清单移出，原阴性scDRS数据/图件保留档案；不保留其支持遗传风险富集的主张。当前正文状态见本阶段qc/MAIN_TEXT_SECTION_STATUS.tsv。
- 唯一工作目录：`phase133_sections2_5_2_7_targeted_support_remote_v1`。本轮新端点先冻结SUPPLEMENT_ANALYSIS_PLAN.md再计算，全部双侧、固定目标/阈值/家族，无挑方法或患者。交付 `C:/Users/31937/Downloads/26210/09_2.5_2.7_targeted_support_v1`：4 PNG+4 PDF+含全部统计的中文说明，9/9 SHA256一致。
- 新2.5补充：20患者TYMP RNA与CASP9/CHMP1A/PRSS8活性的患者内Spearman中位rho=0.076646/0.117349/0.110424，正向18/20、16/20、17/20，P=0.000636/0.001825/0.000834，四蛋白BH FDR=0.001668/0.002434/0.001668。PXN未调整rho=-0.026829，P=FDR=0.562821。
- 技术协变量敏感性：调整原文库量和检出基因数的秩后，CASP9/CHMP1A减弱；PRSS8 rho=0.054238，P=0.045793/FDR=0.091586（名义）；PXN rho=0.046029，16/20正向，P=0.004823/FDR=0.019292。主/敏感性各为预定义四蛋白BH，全部结果展示，不择优模型。RRM2B为独立RNA端点：未调整rho=0.024532、13/20正向、P=0.045793；调整后接近零P=0.695063，不称蛋白活性。
- 新2.7补充：固定13酶基因×3种原拟时序，以患者早20%/末20%细胞构建raw-UMI pseudobulk，按全基因UMI归一化log2(CPM+1)，20独立患者、每臂>=5细胞。39项统一BH。scTour TYMP中位delta=+0.312972、13/20正向、P=0.034919/FDR=0.340461；Slingshot SRM delta=+0.635434、13正6负1零、P=0.025520/FDR=0.331758，均明确名义/方法特异。逐患者留出中位数均为正，范围0.253763–0.372180及0.540414–0.730454。
- 主要反向端点保留：Slingshot PCYT1A delta=-0.689216、P=0.002091/FDR=0.081553；scTour SMS delta=-0.226539、P=0.016042/FDR=0.312820。完整39项图和全部4个名义端点配对图一并交付；不把不同方法的不同基因拼成共激活结论。新pseudobulk端点补充而非替换旧cell-median RNA端点。
- 恢复旧Phase30名义线索：TYMP严格功能代理的F1高主克隆比较n=14，delta=0.057178，11正2负1零，P=0.027708/FDR=0.055416；患者内F1四分位n=20，delta=0.042944，13正5负2零，P=0.031104/FDR=0.062208。
- TYMP广义代理末端：Slingshot delta=0.162888，15/20正向；scTour delta=0.300574，14/20正向；两者P=0.044054/FDR=0.088108。按当前定位可作两方法名义方向支持，原P/FDR及旧等级记录未改。Monocle3对应为负且不显著，完整保留。
- 旧RRM2B广义代理与F1的反向结果完整展示。严格代理虽在Slingshot/scTour有名义P，但各仅5患者正向、14或15零，总体中位数0，不作为普遍共激活。旧字段F1_dominant_group经源代码确认实际为F1高主克隆比较，显示已改正为F1-high clone，避免与原2.5的9患者F1主导细胞比较混淆。
- 复核：200条患者相关（含偏相关）与780条计数推导末段效应由R独立复算，相关最大误差2.51e-15、P最大误差2.23e-16，所有本轮BH通过。旧72项汇总效应/正负零计数及全部原BH复核通过。4/4 PDF单页、边界及逐图视觉核查完成。所有原始/汇总/LOO数据、脚本、方案、源SHA256和记录只在服务器。
- 已向用户请求简要湿实验干预和方向，尚未收到；不要求重新上传完整原始数据。当前完成的是方向无偏的支持性补充，具体实验对应待方向清单，不自行假定湿实验结果。现有Slingshot all_lineages以na=FALSE导出，不代表分支归属；未开展无依据的分支检验。
- 当前断点：用户审阅补充结果并对齐湿实验。2.3正文退出；2.5/2.7新增线索和历史名义证据已准备好，无活动计算。原Phase130–132图件保持。
'''
for n in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
 p=root/n;old=p.read_text()
 if '## Phase133 sections2.5/2.7 supportive supplementation complete' not in old:
  tmp=p.with_name(p.name+'.phase133.tmp');tmp.write_text(note+'\n\n'+old,encoding='utf8');tmp.replace(p)
p=root/'AGENTS.md';old=p.read_text();rule='''\n## 2026-09-09 用户明确的生信证据定位与正文取舍\n\n用户明确：湿实验承担主要验证，生信以方向性关联、患者一致性或名义显著提供对应支持；不要把“独立完成强因果闭环”作为每个生信结果可用的门槛。P<0.05但FDR未通过时标为名义显著，完整P/FDR保留，关键反向结果不能隐藏。\n2.3 scDRS从正文主线/主图清单移出，原数据和阴性记录完整保留；重点补充2.5与2.7。具体湿实验干预和变化方向按用户说明或已提供材料对应，未提供的内容不自行假定。\n'''
if '用户明确的生信证据定位与正文取舍' not in old:p.write_text(old+rule,encoding='utf8')
p=root/'RESULTS_MANIFEST.tsv';old=p.read_text()
with p.open('a') as f:
 for i,row in enumerate(rows):
  key='phase133_targeted_support_'+str(i+1)
  if key+'\t' not in old:f.write(key+'\tphase133_sections25_27_targeted_support\treview_deliverable\t'+str(w/'figures'/row['name'])+'\t'+datetime.datetime.now().isoformat()+'\tPASS_REVIEW_CHECKS\tSHA256 '+row['sha256']+'\t9/9 delivered; 2.3 omitted from main narrative; nominal and FDR evidence distinguished; all opposing results retained\n')
manifest=json.loads((q/'SOURCE_MANIFEST.json').read_text());lookup={r['path']:r for r in manifest}
for p in (w/'scripts').glob('*'):
 if p.is_file():lookup[str(p)]={'path':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
(q/'SOURCE_MANIFEST.json').write_text(json.dumps(list(lookup.values()),indent=2))
(q/'SOFTWARE_VERSIONS.json').write_text(json.dumps({'python':sys.version,'executable':sys.executable,'numpy':numpy.__version__,'pandas':pandas.__version__,'scipy':scipy.__version__,'matplotlib':matplotlib.__version__},indent=2))
state={'phase':133,'status':'PASS_EXPLORATORY_SUPPLEMENT_COMPLETE','section23_main_text':'OMITTED_BY_USER','new_25_patient_correlations':200,'new_27_patient_pseudobulk_effects':780,'restored_prior_summary_rows':72,'figures':4,'delivered_files':9,'sha256_matches':9,'independent_validation':'PASS','visual_review':'PASS','wetlab_direction_details':'AWAITING_USER','active_jobs':False,'next':'Review supportive findings and align exact wet-lab intervention directions','local_delivery':r['local_delivery']}
for n in ['CURRENT_TASK.json','COMPLETION.json']:(q/n).write_text(json.dumps(state,indent=2))
print(json.dumps(state),flush=True)
