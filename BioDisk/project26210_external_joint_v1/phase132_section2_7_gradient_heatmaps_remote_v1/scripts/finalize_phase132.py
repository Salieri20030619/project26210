from pathlib import Path
import json,hashlib,shutil,datetime
w=Path(__file__).resolve().parents[1];root=w.parent;q=w/'qc'
v=json.loads((q/'VALIDATION.json').read_text());assert v['status']=='PASS'
assert json.loads((q/'VISUAL_REVIEW.json').read_text())['status']=='PASS'
receipt=json.loads((q/'DELIVERY_VALIDATION.json').read_text());assert receipt['sha256_matches']==7
rows=json.loads((q/'DELIVERY_MANIFEST.json').read_text())
for r in rows:assert hashlib.sha256((w/'figures'/r['name']).read_bytes()).hexdigest()==r['sha256']
backup=w/'state_backup_before_completion';backup.mkdir(exist_ok=True)
for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv']:
 if not (backup/name).exists():shutil.copy2(root/name,backup/name)
note='''## Phase132 section2.7 gradient heatmaps delivered (2026-09-09)

- 用户指出2.7缺原方案底部image19的渐变热图。本次补齐Monocle3主图及Slingshot/scTour同基因对照，三张均为蓝紫—橙黄渐变表达热图，含左侧树/四簇色条与数量、顶部真实状态构成/拟时序色带、右侧完整图例及12酶基因引线。
- 唯一工作目录：`phase132_section2_7_gradient_heatmaps_remote_v1`。交付 `C:/Users/31937/Downloads/26210/08_2.7_gradient_heatmaps_v1`：3 PNG+3 PDF+中文说明，7/7 SHA256一致；编号11/12/13承接原2.7十套图。每图360 dpi。
- 基因为冻结轨迹2,500 HVG，加11个不在HVG中的原方案酶（TK1已在HVG），共2,511行全部保留，未按显著性或外观选基因。原始HVG counts以全保留基因文库UMI归一化log1p CP10K；新增酶沿用Phase129正式同源提取。20位>=30细胞患者共11,921细胞，另6小样本患者29细胞不参与显示拟合，源对象不改。
- 每位患者内20拟时序秩箱，Gaussian h=0.075按箱细胞数加权拟合，跨20位患者等权均值；200个秩位置0.025–0.975，基因内沿秩z标准化。仅水平插值至1,400显示位置，不跨基因行平滑。三方法原轨迹和方向不改，不重跑模型。
- Monocle3的z曲线Ward/Euclidean树切为4显示簇C1/C2/C3/C4，基因数408/741/340/1022。仅旋转原树分支以改善阅读，树距离不改；三方法完全固定同基因行顺序/簇，便于比较。色阶统一[-2.5,2.5]，越界饱和并标注；C簇非新细胞亚型、非已确认显著动态基因。
- 顶部采用原Base/Hypoxia/Cycle/TYMP-state的患者等权构成比例，每秩点和为1；不把共同秩位置虚构为单一细胞类型。12酶标签全部对应实际行，文字避让不改数据位置。
- 独立加权和公式复核拟合最大误差2.7e-15；全部基因/患者纳入、行序与状态比例核查通过。3/3 PDF单页、边界与逐图视觉检查通过。图为表达动态描述，无新增推断性检验；原2.7患者层相关P/FDR和跨方法不一致结论保留。
- 完整患者曲线/均值/z矩阵、基因登记、行簇与标签位置、来源SHA256、软件和脚本均在服务器。当前断点：用户审阅2.7新增渐变热图；无活动绘图任务。Phase130/131修正及原其他章节完成状态保持。
'''
for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
 p=root/name;old=p.read_text()
 if '## Phase132 section2.7' not in old:
  tmp=p.with_name(p.name+'.phase132.tmp');tmp.write_text(note+'\n\n'+old,encoding='utf8');tmp.replace(p)
p=root/'RESULTS_MANIFEST.tsv';old=p.read_text()
with p.open('a') as f:
 for i,r in enumerate(rows):
  key='phase132_gradient_heatmap_'+str(i+1)
  if key+'\t' not in old:f.write(key+'\tphase132_section27_gradient_heatmap\treview_deliverable\t'+str(w/'figures'/r['name'])+'\t'+datetime.datetime.now().isoformat()+'\tPASS_REVIEW_CHECKS\tSHA256 '+r['sha256']+'\t7/7 delivered; 2511 fixed genes; 20 patients; 3 methods; no new inferential tests\n')
manifest=json.loads((q/'SOURCE_MANIFEST.json').read_text());lookup={r['path']:r for r in manifest}
for p in (w/'scripts').glob('*'):
 if p.is_file():lookup[str(p)]={'path':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
(q/'SOURCE_MANIFEST.json').write_text(json.dumps(list(lookup.values()),indent=2))
state={'phase':132,'section':'2.7','status':'PASS_GRADIENT_HEATMAPS_DELIVERED','figures':3,'methods':['Monocle3','Slingshot','scTour'],'genes':2511,'patients':20,'delivered_files':7,'sha256_matches':7,'PDF_checks':3,'visual_review':'PASS','active_jobs':False,'next':'User reviews section2.7 gradient heatmaps','local_delivery':receipt['local_delivery']}
for name in ['COMPLETION.json','CURRENT_TASK.json']:(q/name).write_text(json.dumps(state,indent=2))
print(json.dumps(state),flush=True)
