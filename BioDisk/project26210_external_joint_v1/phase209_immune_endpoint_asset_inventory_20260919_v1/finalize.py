from pathlib import Path
import csv,json,hashlib,shutil,datetime
h=Path.home();r=h/'project26210_external_joint_v1';o=r/'phase209_immune_endpoint_asset_inventory_20260919_v1'
report='''# 26210 下游免疫终点现有资产盘点

日期：2026-09-19；Phase209。前两部分已交付，以Phase208正文v6及冻结C4为当前入口。本次仅盘点已有文件、细胞元数据、组织/患者覆盖及历史统计，不启动新评分、聚类、通信模型或图件；不推进空间分支。

## 结论：数据基础齐备，主要缺口在接收免疫细胞身份和端点定义

现有75文库主对象包含T/NK及同文库巨噬，冻结C4的2674个细胞ID全部匹配。具备从当前主线进入患者级免疫终点的原始数据基础。主对象仅有T cell与NK cell大类，尚无本主对象内冻结的CD8/CD4/Treg细分身份，因此CD8专属功能、Treg/CD4比例及其有效样本量尚未成立。

旧Phase72已经完成三队列患者级免疫生态分析，数据、脚本、评分、原始P/FDR及验证记录均在。其预测量是旧髓系程序，结局是全T功能障碍、T+NK混合细胞毒和Treg/全T比例；不是新C4到CD8/Treg的正式分析。旧结果可完整复用为背景和补充，不能通过改名成为新终点。

## 1. 当前发现队列：75文库、40患者、5研究

主对象为 /mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/03_malignancy_cnv/cnv_annotated_h5ad_v1/part1_scvi_integrated_cell_annotated_cnv_v1.h5ad 。本轮只读obs、var及层结构，未加载表达矩阵；对象258120细胞、含counts层。历史SHA256登记于文件清单，本轮未重新全量hash。

| 当前标签 | 全组织细胞/患者/文库 | 肿瘤细胞/患者/文库 | 用途 |
|---|---|---|---|
| T cell | 63071 / 40 / 74 | 56716 / 37 / 51 | 接收端母集合，需CD8/CD4/Treg身份核查 |
| NK cell | 23373 / 40 / 71 | 20174 / 36 / 49 | 次要细胞毒端点候选，需核查与细胞数资格 |
| Macrophage大类 | 42546 / 40 / 75 | 37637 / 38 / 52 | 原大类覆盖，不等同净化后成熟巨噬分母 |
| 冻结C4 | 2674 / 37 / 55 | 2437 / 35 / 47 | 沿用成员，不重分群 |
| Monocyte | 8215 / 39 / 70 | 6727 / 36 / 50 | 单核背景，不能混入成熟巨噬预测量 |
| B cell | 2571 / 36 / 64 | 1803 / 33 / 46 | 补充资产，暂不扩展主要端点 |
| Dendritic cell | 1070 / 36 / 55 | 905 / 34 / 44 | 补充资产，暂不扩展主要端点 |

同肿瘤文库有巨噬与T：51文库/37患者；巨噬与NK：49文库/36患者。C4与T、C4与NK均47文库/35患者。上述只要求至少一个细胞存在，不是检验资格，不是最后n。C4缺失但巨噬合格的患者，后续C4比例可按规则保留真实零；不能把35例存在患者自动用作比例模型入选名单。C4内SPP1表达则需要足够C4细胞，分母不同。

25个身份/功能候选基因均存在于特征索引，包括CD3D/E、TRAC、CD4、CD8A/B、FOXP3、IL2RA、NKG7、GNLY、PRF1、PDCD1、TOX等。索引存在不等于在目标细胞中充分表达；本轮没有验证表达覆盖或完成注释。

发现队列GSE207493也出现在历史Phase72的19例队列中，不能将其当作发现队列之外的独立验证。两者样本范围不同，后续须按患者ID明确重叠，不可简单相加。

## 2. 外部单细胞：对象与冻结标签可用，需重新核对目标区室

三份Phase7 Seurat对象在 project26210_phase71_immune_advanced_v1/inputs_v1/seurat/，大小约1.92/0.94/1.46 GB；已有历史传输hash和ID映射核验。本轮核对存在及文件信息，未重新加载RDS。

| 外部队列 | 冻结细分标签存在覆盖 | 既有患者级分析范围 | 缺口 |
|---|---|---|---|
| GSE264586 | CD8_cytotoxic 2338、CD4_naive_memory 15211、NK 13518，均12患者；另5902未解决T细胞 | 12患者、19患者时间点；5例重复取样 | 没有冻结Treg标签；不得把无标签当作生物学零。CD8_cytotoxic是否覆盖非细胞毒/功能障碍CD8须核查 |
| GSE294109 | CD8_cytotoxic 46039、CD4_naive_memory 13138、Treg 6341、NK 12312，均15患者 | Phase72细胞数门槛后13患者 | 15例存在覆盖不是13例分析资格；CD4分母及CD8标签选择性需核查 |
| GSE207493 | 旧版有CD4_T/CD8_T/Treg/Exhausted_T/Proliferating_T/NK | Phase72 19患者 | 与当前发现样本重叠；不能直接把旧标签迁到75库主对象 |

外部细胞级映射路径：project26210_external_joint_v1_phase67_direct_program_replication_v1/inputs_v1/metadata_v2/ 。这里的metadata可用，不代表Phase67目录下三次失败的主程序已成功；失败标记与.part结果保留，不能引用为正式结果。

GSE264586存在Baseline/Cycle2/End_of_treatment，治疗字段为sitravatinib+nivolumab+ipilimumab，response均Unknown。GSE294109同时有未治疗与治疗后手术，多个方案混合。字段来自冻结元数据，本轮未重新审核原始临床来源；不能将全队列解释为单一治疗暴露，也不能把时间点扩充成患者n。

## 3. 已有免疫统计：完整保留，不能升级为C4结论

权威目录：project26210_phase72_immune_ecosystem_v1 。既有结果229326免疫细胞，患者单位19/12/13。脚本门槛为每患者时间点至少50 Myeloid/DC和50 T细胞，跨时间点取患者中位数。该门槛不自动适用于新的CD8/CD4端点。

- 旧SPP1程序→T功能障碍：发现rho=-0.418，P=0.0753，FDR=0.3388；外部12例rho=0.510，P=0.0899，FDR=0.1903；外部13例rho=0.242，P=0.4262，FDR=0.6347。不能称稳定外部支持。
- 旧SPP1程序→细胞毒、SPP1→Treg等未获预设支持；9条主关系无一达到双外部规则。既有关键反向和无支持结果保留。
- 髓系免疫调节程序→T功能障碍有合并方向证据：两外部固定效应rho=0.614、q=0.0109；单队列q=0.1903/0.1488，其中13例P=0.01654为名义显著。复合生态端点合并rho=0.522、q=0.0347。合并结果不能替代单队列验证，也不能转称C4效应。
- GSE264586纵向仅5例，GSE294109疗效探索8例CR/PR对3例No_effect，历史无支持记录保留。

Phase71是髓系内部代理程序分析，不是T/NK接收端结果；其MTA_RESPONSE不等于实测MTA或本次C4定义。冻结程序、基因成员及评分代码可复用作历史复现，不作为新增免疫终点的预设答案。

本轮核对历史验证摘要8/8通过及完整统计表可读，没有重新计算原P/FDR。既有统计原样保存在existing_evidence。

## 4. 通信、临床和空间资产的分工

通信：phase35_part1_c01_c04_cell_communication_v1有原CellChat/NicheNet脚本及患者结果；Phase109 analysis_v3及Phase144审计有旧SPP1-state细分通信。已确认的重点是内皮联系，并非新C4→CD8机制。先复用配体受体资源、患者映射和原始结果；新C4及目标T/NK身份确认后再判断是否需要定向追加。旧失败版本不使用。

临床：Phase64、69及Phase81/82保留CheckMate、JAVELIN等表达/临床资产。Phase81记录CheckMate旧Phase64 PFS/OS事件编码倒置，相关Cox/治疗交互已被修正版替代；JAVELIN原事件规则仍正确。修正版RNA资产CheckMate n=311、JAVELIN n=726。不能把旧生存显著结果拿来补强免疫功能，也不能把bulk表达指定给C4。当前免疫终点盘点不重做临床模型。

空间：旧SPP1-state的血管或FFPE T功能障碍结果属于另一分支，不能填补免疫端点阴性或身份缺口；本轮不重跑空间模型。

## 5. 下一步可执行顺序（未在本轮启动）

1. 基于当前75文库冻结对象核对T/NK谱系，建立CD8、CD4、Treg及未定/混合标签；优先检查可追溯旧注释的映射与marker证据，不按关联强弱挑选细胞。保留C4成员和前两部分。
2. 用净化后的成熟巨噬建立整体SPP1、C4比例和C4内SPP1三种独立预测量；按同患者同肿瘤文库列出CD8/CD4/NK端点的细胞数、表达覆盖和可估计性，固定多文库规则。
3. 固定CD8细胞毒、CD8功能障碍、Treg/CD4三个优先读出；NK独立为次要读出。先区分活化与功能障碍，Treg分母明确包括哪些CD4状态。实际检验家族按同类假设定义，不能把比例与功能等全部无关检验混成一个池。
4. 完成资格后才进行患者层关联、有限研究调整及外部验证；通信作为候选机制对应。n<10按探索性口径，完整P/FDR及关键反向结果保留。

本次交付止于盘点；上述是后续建议，不代表已完成或已启动。

## 6. 文件与核查边界

- ASSET_FILE_INVENTORY.tsv：638项文件信息，包含主对象；小于2 MB文件记录本轮SHA256，大文件仅元数据，历史hash明确区分。
- PATIENT_LIBRARY_LINEAGE_COUNTS.tsv：患者×文库×组织×大类/C4数量。
- DISCOVERY_LINEAGE_COVERAGE.tsv、DISCOVERY_TUMOR_CO_PRESENCE.tsv：覆盖及同文库存在性。
- EXTERNAL_FINE_LABEL_PRESENCE.tsv及两份ANNOTATION_COUNTS.tsv：外部现有标签与时间/治疗字段。
- IDENTITY_ENDPOINT_GENE_INDEX_COVERAGE.tsv：25基因索引存在性，不是表达有效性。
- existing_evidence/：Phase72原始统计及证据账本、Phase81纠错说明。
- inventory.py、inventory.log、INVENTORY_STATUS.json及FINAL_VALIDATION.json：代码、执行退出状态和完整性检查。

本盘点不等同对所有表达对象、细分注释、临床来源或历史全流程的重新验证。已有资产和未完成的资格审查分开登记。
'''
(o/'IMMUNE_ASSET_INVENTORY_CN.md').write_text(report)
e=o/'existing_evidence';e.mkdir(exist_ok=True)
for p in [h/'project26210_phase72_immune_ecosystem_v1/results_v1/phase72_primary_relations_all_v1.tsv',h/'project26210_phase72_immune_ecosystem_v1/PHASE72_IMMUNE_EVIDENCE_REGISTER_v1.tsv',h/'project26210_phase72_immune_ecosystem_v1/results_v1/phase72_cross_cohort_decisions_v1.tsv',h/'project26210_phase81_bulk_ici_transcriptome_rebuild_v1/PROJECT_STATUS_PHASE81_82_v1.md']:
 shutil.copy2(p,e/p.name)
status=json.loads((o/'INVENTORY_STATUS.json').read_text())
with (o/'DISCOVERY_LINEAGE_COVERAGE.tsv').open() as f:cov=list(csv.DictReader(f,delimiter='\t'))
checks={'source_cells':status['source_cells']==258120,'C4_full_match':status['C4_ids_matched']==2674,'metadata_cell_total':sum(int(x['cells']) for x in cov if x['scope']=='All' and x['lineage']!='C4_frozen')==258120,'report_exists':(o/'IMMUNE_ASSET_INVENTORY_CN.md').stat().st_size>0,'no_models':status['no_new_models']}
assert all(checks.values()),checks
(o/'FINAL_VALIDATION.json').write_text(json.dumps(checks,indent=2))
backup=o/'root_record_backups';backup.mkdir(exist_ok=True)
for name in ['PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv','CODEX_HANDOFF.md']:
 shutil.copy2(r/name,backup/name)
note='''\n## 2026-09-19 Phase209：免疫终点资产盘点完成\n用户确认前两部分已交付，后续进入免疫终点，本轮仅盘点。主对象75库40患者258120细胞，肿瘤T56716/37患者、NK20174/36患者；C4肿瘤2437/35患者，2674全体ID完全匹配。存在性不是端点资格。当前主对象缺CD8/CD4/Treg细分标签；旧Phase72是全T、T+NK及Treg/全T口径，9关系无双外部支持，保留合并方向及关键阴性；不能改名为C4结果。外部GSE264586无冻结Treg，GSE294109存在15患者但旧模型仅13；GSE207493与发现重叠。已标注Phase81替代旧CheckMate事件错误模型。638文件盘点、元数据覆盖表、完整旧统计与中文报告已保存；未启动新模型/新图/空间。\n权威：phase209_immune_endpoint_asset_inventory_20260919_v1/IMMUNE_ASSET_INVENTORY_CN.md。下一步建议接收端身份及同肿瘤文库资格核查，待后续任务；前两部分继续冻结。\n'''
for name in ['PROJECT_STATUS.md','DECISIONS.md']:
 with (r/name).open('a') as f:f.write(note)
p=r/'CODEX_HANDOFF.md';p.write_text(note+'\n'+p.read_text())
now=datetime.datetime.now().astimezone().isoformat()
with (r/'RESULTS_MANIFEST.tsv').open('a') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n')
 for p in sorted(o.iterdir()):
  if p.is_file():w.writerow(['phase209_'+p.stem,'phase209','asset_inventory',str(p),now,'PASS_METADATA_SCOPE',str(o/'FINAL_VALIDATION.json'),'No new biological model; existing result definitions retained'])
manifest={str(p.relative_to(o)):hashlib.sha256(p.read_bytes()).hexdigest() for p in o.rglob('*') if p.is_file() and 'root_record_backups' not in p.parts}
(o/'OUTPUT_SHA256.json').write_text(json.dumps(manifest,indent=2))
print('FINALIZED',len(manifest),'artifacts',checks)
