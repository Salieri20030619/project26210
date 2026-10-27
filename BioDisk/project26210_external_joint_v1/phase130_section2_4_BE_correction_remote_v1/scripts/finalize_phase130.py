from pathlib import Path
import json,hashlib,shutil,datetime
w=Path(__file__).resolve().parents[1];root=w.parent;q=w/'qc'
assert json.loads((q/'VALIDATION.json').read_text())['status']=='PASS'
assert json.loads((q/'E_NETWORK_VALIDATION.json').read_text())['status']=='PASS'
assert json.loads((q/'VISUAL_REVIEW.json').read_text())['status']=='PASS'
receipt=json.loads((q/'DELIVERY_VALIDATION.json').read_text());assert receipt['sha256_matches']==11
rows=json.loads((q/'DELIVERY_MANIFEST.json').read_text())
for r in rows:assert hashlib.sha256((w/'figures'/r['name']).read_bytes()).hexdigest()==r['sha256']
backup=w/'state_backup_before_completion';backup.mkdir(exist_ok=True)
for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv']:
 if not (backup/name).exists():shutil.copy2(root/name,backup/name)
note='''## Phase130 section2.4 B/E corrections delivered (2026-09-09)

- 用户先指出2.4B遗漏左侧标记、顶部细胞类型信息和图例，随后明确纠正2.4E：参考C/D/E是基因间连接，旧图错误画成基因与方法之间的连接。本阶段同时处理两项反馈。
- 唯一工作目录：`phase130_section2_4_BE_correction_remote_v1`。交付 `C:/Users/31937/Downloads/26210/05_2.4_B_E_correction_v2`：B热图、E基因网络总图和ILF2/HOXD3/HNF4A三个单图，每图360 dpi PNG+PDF，加中文说明共11文件；11/11 SHA256一致。
- B补齐原树k=4显示聚类的C1(9)、C2(8)、C3(4)、C4(3)、括号/色条；顶部明确Cell type=Malignant epithelial、cNMF program=F1–F7、Algorithm=G3/GB2，并提供完整独立图例。Cluster标题与Algorithm标签间距已修正。聚类为调控子显示分组，不是新细胞亚型。
- B原24调控子、14列、行列顺序及热图数值不变（患者平衡汇总复算误差<1e-15）；21位患者/37份患者×程序谱，F1–F7 n=9/12/2/4/4/2/4。未重跑pySCENIC/AUCell，未新增统计。
- **旧E图形结构撤回，由v2替代。** 旧图将GENIE3/GRNBoost2设为连线端点，错误表达了真实TF→靶基因关系。v2全部节点为基因，连线直接从ILF2/HOXD3/HNF4A指向靶基因；算法只在边颜色与图例中表示支持来源。
- E总图使用原展示68个靶基因+3个TF共71个节点，保留原68条唯一基因边，并在同一节点集合内补齐正式源表存在的3条边：ILF2→PTER、HNF4A→PTN、HNF4A→TRPV4。总计71条唯一有向基因边/89条方法支持记录；双方法18、GENIE3特有20、GRNBoost2特有33。无算法节点，无虚构靶基因间边。
- E单图ILF2/HOXD3/HNF4A分别15/24/32靶基因；共享靶点PTER/PTN/TRPV4在总图各仅一个节点。固定原展示节点选择，增加边仅因补齐该节点集合内真实连接，不扩大或重新筛选基因。图为预测调控子子网络，不是实验物理互作或完整调控网络；未新增TYMP直接调控主张。
- B数值/标签复核、E全部边源表核查、5/5 PDF单页与边界检查、5/5逐图视觉审核均通过。初稿/旧版本保留服务器；完整矩阵、节点/边清单、原始方法权重、哈希、脚本与软件记录齐全。
- 当前断点：用户审阅2.4 B/E v2。2.6/2.7 Phase129完成状态保持；本轮图件无活动任务，不自动推进其他章节。后续进度继续仅写服务器。
'''
for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
 p=root/name;old=p.read_text()
 if '## Phase130 section2.4 B/E' not in old:
  tmp=p.with_name(p.name+'.phase130.tmp');tmp.write_text(note+'\n\n'+old,encoding='utf8');tmp.replace(p)
p=root/'RESULTS_MANIFEST.tsv';old=p.read_text()
with p.open('a') as out:
 for i,r in enumerate(rows):
  key='phase130_section24_BE_'+str(i+1)
  if key+'\t' not in old:out.write(key+'\tphase130_section24_BE_correction\treview_deliverable\t'+str(w/'figures'/r['name'])+'\t'+datetime.datetime.now().isoformat()+'\tPASS_REVIEW_CHECKS\tSHA256 '+r['sha256']+'\t11/11 delivered; B annotations restored; E gene-to-gene topology supersedes old method-hub graphic\n')
manifest=[]
for name in ['SOURCE_MANIFEST.json','E_SOURCE_MANIFEST.json']:manifest+=json.loads((q/name).read_text())
lookup={r['path']:r for r in manifest}
for p in (w/'scripts').glob('*'):
 if p.is_file():lookup[str(p)]={'path':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
(q/'SOURCE_MANIFEST.json').write_text(json.dumps(list(lookup.values()),indent=2))
state={'phase':130,'sections':['2.4B','2.4E'],'status':'PASS_CORRECTIONS_DELIVERED','PDF_figures':5,'delivered_files':11,'sha256_matches':11,'B_matrix_unchanged':True,'E_gene_nodes':71,'E_gene_edges':71,'E_method_nodes':0,'old_E_figure_superseded':True,'active_jobs':False,'next':'User reviews section2.4 B/E v2','local_delivery':receipt['local_delivery']}
for name in ['CURRENT_TASK.json','COMPLETION.json']:(q/name).write_text(json.dumps(state,indent=2))
print(json.dumps(state),flush=True)
