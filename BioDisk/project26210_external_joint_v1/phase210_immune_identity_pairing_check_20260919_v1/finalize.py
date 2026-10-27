from pathlib import Path
import pandas as pd,json,shutil,hashlib,datetime,csv
H=Path.home();R=H/'project26210_external_joint_v1';O=R/'phase210_immune_identity_pairing_check_20260919_v1'
assert json.loads((O/'VALIDATION.json').read_text())['passed']
q=pd.read_csv(O/'SAME_TUMOR_LIBRARY_READINESS.tsv',sep='\t');e=pd.read_csv(O/'EXTERNAL_PATIENT_SAMPLE_READINESS.tsv',sep='\t');extra=[]
for (ds,tp),g in e.groupby(['dataset','timepoint']):
 for ep in ['CD8','NK','Treg_CD4']:
  a=g[g[ep+'_count_ready']];extra.append(dict(dataset=ds,timepoint=tp,endpoint=ep,patients=a.patient.nunique(),samples=len(a)))
pd.DataFrame(extra).to_csv(O/'EXTERNAL_TIMEPOINT_READINESS.tsv',sep='\t',index=False)
report='''# 26210 免疫接收端身份与配对资格核对

Phase210，2026-09-19。承接Phase209；前两部分与冻结C4保持原定义。用户本轮授权核对身份与配对；本轮未重聚类、未建立新免疫关联模型、未绘图、未更改源注释。

## 结论：精确映射与数据配对成立，正式细分身份仍需整理

86444个当前T/NK细胞全部精确匹配到旧CellTypist逐细胞预测，研究与文库字段一致；原始counts抽取的文库总量全部吻合。冻结C4的2674个ID完整匹配。数据源及配对基础可用。

但旧标签不能直接升为当前正式注释：存在T/NK交叉标签、Treg逐细胞检出不足和外部功能障碍T群身份未定。已交付候选标签、标志证据、冲突/未定细胞及逐患者文库资格，不把严格标志阳性子集当作全体T细胞。

## 1. 数据身份与成熟巨噬分母已核实

- 当前发现对象：/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/03_malignancy_cnv/cnv_annotated_h5ad_v1/part1_scvi_integrated_cell_annotated_cnv_v1.h5ad。
- 旧预测：/mnt/BioDisk/ludingyu/project26210/object/01_annotation_adata_annotated.h5ad。当前T/NK全部有精确ID，不做条码截短或模糊匹配。旧预测提供线索，不是独立验证。
- 成熟巨噬继承Phase154口径：part2_rebuild_v1/08_regulons/input_v1/CELL_METADATA.tsv.gz排除candidate 1和Monocyte；共22444细胞（全组织），不是原Macrophage大类42546细胞。C4成员2674不变，与T/NK无ID交集。
- 75文库、患者×研究身份、组织来源保留。只在同肿瘤文库检查巨噬与接收端；同患者多个文库不得增加统计n。

## 2. 当前T/NK需要处理的具体问题

肿瘤T大类中，旧预测有CD8候选30700、CD4候选13525、Treg候选5640、NK候选3238、非经典T1911、增殖未定T35、其余未定1667。旧CD4候选包括明确的CD4细胞毒标签，未按“cytotoxic”字符串全部归入CD8。γδT、MAIT、NKT独立保留。

肿瘤NK大类中，旧预测有4440个CD8候选，其中1517满足本轮T骨架和CD8标志检出规则；当前T大类中又有3238个旧NK候选。两类交叉不能直接决定谁是正确标签，须结合表达/群体结构核实；正式资格筛查暂要求旧候选与当前大类相容，冲突细胞完整列出而不重命名源对象。

T骨架采用CD3D/CD3E/TRAC至少2个检出，CD8额外要求CD8A或CD8B。CD4/IL7R只用于已属helper候选的支持，不把IL7R单独当作CD4身份。NK采用KLRD1/NCR1/KLRF1/NCAM1/FCGR3A面板，并检查T骨架。

Treg候选5640个中，857同时有T骨架、FOXP3及IL2RA或CTLA4检出；再排除本轮混合谱系复核标记后856个。该差异不能解释为其余均不是Treg：单细胞掉零、状态和测序深度会影响检出。必须看群体标志与一致性，不能仅保留FOXP3检出的细胞后宣布获得无偏Treg比例。

原多数投票与原逐细胞预测不是独立算法验证。严格标志支持子集用于检验注释可信度与缺口，不作为正式免疫统计的默认入选集合；尤其不能按细胞毒或功能障碍分数筛选CD8后再研究该功能端点。

## 3. 同肿瘤文库的候选资格与基因覆盖

以下规则在读取结果前登记为技术筛查：成熟巨噬≥50；CD8/NK候选≥20；Treg/CD4的CD4候选分母≥50。另完整保留10/20/50细胞门槛对照，不按统计显著性挑阈值。CD4分母为相容的CD4候选+Treg候选；不混入NK、γδT、MAIT或未定细胞。比例允许零C4；C4内表达另需≥20 C4。

| 拟定免疫端点 | 计数资格：文库/患者 | 加程序表达覆盖后：文库/患者 | 当前含义 |
|---|---|---|---|
| CD8细胞毒 | 42 / 31 | 42 / 31 | 5研究，候选身份与数据覆盖可衔接 |
| CD8功能障碍 | 42 / 31 | 32 / 25 | 5研究，部分文库功能障碍基因覆盖不足 |
| NK细胞毒 | 37 / 26 | 37 / 26 | 4研究，候选独立次要端点 |
| Treg/CD4 | 37 / 27 | 未将逐细胞FOXP3检出当作比例资格 | 4研究，主要缺口为可靠分子分母身份 |

程序表达覆盖仅要求至少3个程序基因各在至少5个候选细胞中检出；本轮未计算功能评分或关联。细胞毒面板NKG7/GNLY/GZMB/PRF1/CTSW；功能障碍面板PDCD1/CTLA4/LAG3/HAVCR2/TIGIT/TOX/CXCL13/ENTPD1。它们作为继承旧程序的覆盖探针，不意味着已冻结最终功能障碍评分；需进一步与活化区分。

严格标志支持子集对应的计数患者数分别为CD8 25、NK 24、CD4分母17。差异说明身份处理会改变资格，不可只挑得到较强关联的版本。以上全部是候选技术资格，不能作为已完成正式注释后的最终n。

目前达到上述成熟巨噬及接收端门槛的文库均有C4，但算法没有先筛掉零C4；完整零值及未达标文库留在表中。C4比例、整体成熟巨噬SPP1与C4内SPP1的资格分开，不共用一个固定n。

患者汇总建议：先固定同肿瘤文库集合，再患者内对合格文库等权汇总；不让某端额外的非配对文库混入；本轮只登记该规则，尚未计算患者功能效应。

## 4. 外部注释来源已核对，GSE264586存在重要端点缺口

已读取Phase7 local_v2最终cluster映射、原marker/规范面板和人工复核说明；这些是分队列独立marker加规范面板注释，不是原作者标签直接导入。本轮复用原群级表达表，未重新加载外部全量RDS或重审原始临床材料。

**GSE264586：CD8功能障碍覆盖不完整。** 5902个T细胞处于cluster4，原记录明确TOX富集、Exhaustion_candidate、低置信、CD4/CD8/Treg身份不足且患者集中。现有CD8标签只有cluster15（1384细胞，增殖）和cluster21（954细胞，干扰素响应）。因此不能把旧CD8标签集合视为完整CD8功能障碍人群，也不能把TOX群直接改成CD8。无冻结Treg，不能给Treg/CD4做零填充。

**GSE294109：更适合优先承接。** CD8 cluster1保留16636个带功能障碍线索的细胞，有CD8A/B、IFNG及LAG3/PDCD1依据；Treg cluster5为6341个高置信细胞。CD8 cluster7的5370个细胞虽有CD8证据，却低置信且接近单患者集中，必须单列复核。增殖、干扰素标签是状态字段，不能自动删除或与谱系混同。

按既有外部标签、同sample/timepoint及成熟巨噬C1QC+SPP1≥50筛查：GSE264586 CD8计数资格8患者/13个样本时间点、NK9/14，Treg不可估计；GSE294109三个候选端点均13患者/22个样本时间点。它们不是新增独立患者，且治疗/取材资格仍需按预设研究问题分开；各时间点人数另表，不混为单一未治疗样本集。外部没有本项目冻结C4，不能将Macrophage_SPP1重命名为C4。

GSE207493同时存在Phase64旧细胞标签与Phase7群级注释两套历史入口；本轮读取的Phase7细分与Phase72所用旧细分不完全相同，不把二者混成同一注释版本。其与发现队列的患者重叠也使它不能直接计作新增独立外部验证。

## 5. 执行结论与下一步

核对工作完成，数据和配对可追溯。正式分析前需要完成接收端细分注释整理：优先处理发现队列T/NK交叉、CD4/Treg群体证据，以及GSE264586未定TOX群和GSE294109低置信患者集中群。采用身份证据决定标签，不依据免疫关联方向决定取舍。

CD8细胞毒、NK细胞毒已有较好的候选技术覆盖；CD8功能障碍需按覆盖保留不同资格；Treg/CD4必须先确认分子分母。可优先准备GSE294109外部衔接；GSE264586保持有缺口的补充队列定位。当前不启动关联、通信、治疗或空间模型。

## 6. 验证与文件

11项一致性检查通过：唯一ID、T骨架、Treg规则、严格集合的谱系相容、共享NK适配蛋白未误作髓系混合证据、75文库、C4总数、成熟巨噬总数、两端门槛及患者去重。原始counts块读取通过整数/非负/有限值及总量核对。

初版把FCER1G/TYROBP纳入混合谱系提示，这两者与NK共享，已纠正为不据它们判断NK混合；初版及所有输出保存在attempt_v1。正式脚本check_v2复用已核对的原始基因缓存，不重跑归一化或聚类。候选配对同时修正为要求与当前T/NK大类相容。首轮校验仅在写JSON时遇到numpy布尔类型序列化错误，external_and_verify_v2转换原生布尔后成功；两轮日志均保留。

主文件：
- CELL_IDENTITY_AUDIT.tsv.gz：每细胞旧预测、原大类、候选、标志检出、冲突与严格支持。
- CELL_PANEL_RAW_COUNTS.tsv.gz、MARKER_EVIDENCE_BY_CANDIDATE.tsv：原始计数与群体表达证据。
- ALL_LIBRARY_IDENTITY_COUNTS.tsv、SAME_TUMOR_LIBRARY_READINESS.tsv、PATIENT_PAIRING_REGISTER.tsv：完整患者/文库和不合格项。
- PROGRAM_EXPRESSION_COVERAGE.tsv、PROGRAM_READINESS_SUMMARY.tsv：程序基因覆盖及候选n。
- EXTERNAL_TNK_CLUSTER_REVIEW.tsv、EXTERNAL_PATIENT_SAMPLE_READINESS.tsv、EXTERNAL_TIMEPOINT_READINESS.tsv：外部来源、资格和时间点。
- PLAN_BEFORE_CHECK.json、INPUT_MANIFEST.json、SOFTWARE.json、VALIDATION.json及执行日志：方法、身份和验收。

局限：候选注释不等于最终生物学身份；逐细胞marker不检出不等于阴性身份；计数/覆盖资格不等于统计支持。本轮不计算P或FDR，不产生免疫机制结论。
'''
(O/'IMMUNE_IDENTITY_PAIRING_REVIEW_CN.md').write_text(report)
# Preserve key source annotation evidence hashes.
base=R/'local_records_migration_20260908_150455/project';paths=[base/'phase7_annotation_local_v2/ANNOTATION_SUMMARY.md',base/'phase7_annotation_local_v2/runs/20260810_234538_376/tables/final_cluster_annotations.tsv']
(O/'EXTERNAL_SOURCE_SHA256.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},indent=2))
status=json.loads((O/'STATUS.json').read_text());status['final_validation_passed']=True;status['completed_at']=datetime.datetime.now().astimezone().isoformat();status['authoritative_scripts']=['check_v2.py','external_and_verify_v2.py','finalize.py'];status['formal_annotation_frozen']=False;(O/'STATUS.json').write_text(json.dumps(status,indent=2))
backup=O/'root_record_backups';backup.mkdir(exist_ok=False)
for name in ['PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv','CODEX_HANDOFF.md']:shutil.copy2(R/name,backup/name)
note='''\n## 2026-09-19 Phase210：免疫身份与配对核对完成\n86444 T/NK细胞精确对应旧CellTypist，研究/文库/原始count总量一致；成熟巨噬按Phase154净化口径22444全组织细胞，C4 2674不变。存在T/NK交叉预测、Treg掉零，旧标签只作候选，未冻结正式细分。相容候选同肿瘤文库计数+覆盖：CD8细胞毒31患者/42库、CD8功能障碍25/32、NK细胞毒26/37；Treg/CD4候选计数27/37，身份分母待整理。数字为技术资格，不是最终n。外部Phase7群级证据核对：GSE2645865902个TOX富集T细胞身份未定，CD8标签主要增殖/干扰素，且无Treg；GSE294109保留功能障碍CD8/Treg，但有低置信单患者集中群。外部计数资格CD8 8/13患者（GSE264586/GSE294109），不能把样本时间点当患者。11项检查通过，初版共享NK标志复核修正及JSON写入错误完整留档。\n权威：phase210_immune_identity_pairing_check_20260919_v1/IMMUNE_IDENTITY_PAIRING_REVIEW_CN.md；正式check_v2、external_and_verify_v2。核对任务完成，尚未开展新细分注释/免疫关联/通信/空间。下一步应整理接收端身份与冲突群，保持C4及前两部分冻结。\n'''
for name in ['PROJECT_STATUS.md','DECISIONS.md']:
 with (R/name).open('a') as f:f.write(note)
p=R/'CODEX_HANDOFF.md';p.write_text(note+'\n'+p.read_text())
manifest={str(p.relative_to(O)):hashlib.sha256(p.read_bytes()).hexdigest() for p in O.rglob('*') if p.is_file() and 'root_record_backups' not in p.parts and 'attempt_v1' not in p.parts}
(O/'OUTPUT_SHA256.json').write_text(json.dumps(manifest,indent=2))
with (R/'RESULTS_MANIFEST.tsv').open('a') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n')
 for p in sorted(O.iterdir()):
  if p.is_file() and p.name not in ['check.py','check.log','external_and_verify.py','external_and_verify.log']:
   w.writerow(['phase210_'+p.stem,'phase210','identity_pairing_audit',str(p),status['completed_at'],'PASS_AUDIT_NOT_FORMAL_ANNOTATION',str(O/'VALIDATION.json'),'No new association; provisional identities distinguished'])
print('PHASE210_FINALIZED',len(manifest))
print(pd.DataFrame(extra).to_string(index=False))
