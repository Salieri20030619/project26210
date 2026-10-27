from pathlib import Path
import json,datetime,hashlib,shutil,csv
import pandas as pd
H=Path.home();R=H/'project26210_external_joint_v1';O=R/'phase211_immune_receiver_refreeze_20260919_v1';F=O/'freeze_v2';validation=json.loads((O/'VALIDATION.json').read_text());assert validation['status']=='PASS';assert json.loads((F/'STATUS.json').read_text())['status']=='FROZEN_VALIDATED'
d=pd.read_csv(F/'DISCOVERY_RECEIVER_CELL_FREEZE.tsv.gz',sep='\t',index_col='cell_id');e=pd.read_csv(F/'DISCOVERY_ENDPOINT_LIBRARY_FREEZE.tsv',sep='\t');ee=pd.read_csv(F/'EXTERNAL_ENDPOINT_SAMPLE_FREEZE.tsv',sep='\t')
comp=e[e.eligible].groupby(['endpoint','study']).agg(patients=('patient','nunique'),libraries=('library','size')).reset_index();comp.to_csv(F/'DISCOVERY_ENDPOINT_STUDY_COMPOSITION.tsv',sep='\t',index=False)
summary=d.groupby(['source_tissue','receiver_label_v1']).agg(cells=('receiver_label_v1','size'),patients=('curated_patient_id','nunique'),libraries=('sample_id','nunique')).reset_index();summary.to_csv(F/'DISCOVERY_FROZEN_LABEL_COVERAGE.tsv',sep='\t',index=False)
transition=pd.crosstab(d.cell_type_broad_v1,d.receiver_label_v1);transition.to_csv(F/'OLD_BROAD_TO_FROZEN_IDENTITY.tsv',sep='\t')
report='''# 26210 免疫接收端整理与重新冻结

Phase211，2026-09-19。正式版本：freeze_v2。前两部分继续以已交付Phase208为准，C4保持原2674个成员及原文件SHA256。本轮完成接收端身份整理、逐细胞标签、同肿瘤文库配对、患者端点资格及外部分层冻结；未计算免疫关联、P/FDR、通信或空间模型，未绘图。

## 一、正式冻结已完成，接收端与C4可按固定名单衔接

全部86444个发现队列T/NK候选细胞都进入新标签表，没有静默丢弃。新的常规CD8、CD4、Treg与NK以群体身份及多基因证据确定，不按逐细胞FOXP3/CD8检出或细胞毒/功能障碍评分挑选细胞。

| 冻结身份（全组织） | 细胞数 | 处理 |
|---|---:|---|
| CD8_T | 31480 | 常规CD8功能接收端 |
| CD4_Tcon | 13816 | CD4常规T细胞，进入CD4分母 |
| CD4_Treg | 3956 | 完整Treg群进入分子与CD4分母 |
| NK | 12520 | 常规NK独立次要接收端 |
| gdT_enriched | 6562 | γδT富集群单列，不混入常规CD8/NK |
| T_unresolved | 10838 | T身份可见、CD4/CD8等未定，保留完整记录 |
| T_NK_unresolved | 7031 | T/NK混合或低深度不足以定类 |
| ILC_candidate | 140 | ILC候选，暂不计入常规T/NK端点 |
| non_TNK_unresolved | 101 | 非T/NK线索，保留未定，不强行命名pDC |

Treg群有连贯的FOXP3、IKZF2、IL2RA/CTLA4与T骨架证据。保留全群3956细胞，包括逐细胞未检出FOXP3的成员，避免把掉零变成比例筛选条件。CD4分母固定为CD4_Tcon+CD4_Treg。NK群不再吸收具有TCR/CD3和γδT特征的群。

## 二、方法与身份决策可追溯

使用原75文库对象的30维scVI表征，仅在冻结T/NK候选子集上建立k=30近邻图；Leiden分辨率1.0、主种子262211形成20群，种子42形成19群。不重训scVI、不修改源对象或C4。群体标志来自原始counts的检出率及log1pCP10K均值，旧CellTypist预测只作支持证据。

两套分群全局ARI=0.7357，不把所有小群称为稳定亚型。主群3在另一分区中分成两个具有清晰相反身份标志的部分：3.0有2518细胞，IL7R检出76.4%、CD8A/B仅3.4%/3.3%，归CD4；3.8有2629细胞，CD8A/B检出48.3%/50.2%、CD4仅2.3%，归CD8。余下成员保留未定。主群16虽有不同种子片段，但标志未形成足够清晰的身份，因此不据片段强制归类。种子42在此辅助混合群拆分，不冒充独立验证。

群体身份的正式映射见DISCOVERY_UNIT_ANNOTATION.tsv；每细胞保留父群、子单元、置信度、判断理由、患者集中标记及分区重叠情况。患者集中与身份是不同问题，不单凭一个患者占比高就否定已有CD8身份；后续可按冻结标记做排除该群的单独检查。

## 三、同肿瘤文库及患者端点名单已固定

主要技术资格固定如下：
- 成熟巨噬≥50；继承Phase154净化口径，排除candidate1与Monocyte，全组织成熟巨噬22444细胞。
- CD8或NK≥20；Treg/CD4要求已分类CD4分母≥50。分子为0可保留，不要求每个患者都检出Treg。
- 每文库未解决T/NK身份≤20%；10%与30%对照完整保留，不用于替代主名单。
- 功能面板至少3个基因各在至少5个接收细胞中检出。此项只核对可测性，尚未计算评分。
- C4比例与整体巨噬SPP1允许真实零值；C4内SPP1单独要求≥20 C4，不将该门槛施加于另外两个预测量。

| 发现队列端点 | 合格文库 | 独立患者 | 研究数 | 可用于C4内SPP1的患者 |
|---|---:|---:|---:|---:|
| CD8细胞毒 | 19 | 16 | 4 | 11 |
| CD8功能障碍 | 19 | 16 | 4 | 11 |
| NK细胞毒 | 16 | 13 | 2 | 10 |
| Treg/CD4 | 16 | 13 | 2 | 10 |

这是本次冻结后的正式技术资格，取代Phase209的存在人数和Phase210的候选人数。它们不是已经完成的免疫效应或支持证据。研究构成另表；NK与Treg只有2个研究，不能套用CD8的4研究覆盖。

未知身份门槛10%/20%/30%下，CD8两端点患者数15/16/16、NK13/13/13、Treg13/13/14。主规则20%在标签审核前登记，未依据关联或P值选择。

患者内必须使用同一端点合格文库集合，预测量和接收端不能各自掺入额外文库；患者内文库等权汇总，一患者一个端点值。具体名单、排除原因及C4内表达单独名单已逐行冻结。

## 四、外部旧NK标签已纠正为未定，不能沿用旧比例

已重新载入两份Seurat对象，精确核实Phase7群号、细胞ID、患者、样本、治疗和原始计数。GSE264586冻结36969个T/NK候选、GSE294109冻结77830个；所有细胞都保留在新表。

GSE264586旧NK群0/23中，分别85.1%/81.0%细胞满足至少2个T骨架基因检出，T骨架与CD8B共同检出42.8%/40.5%，并有γδT信号。GSE294109旧NK群3/13的T骨架比例80.9%/99.1%，群13同时有较强γδT信号。故这些群改为T_NK_unresolved或gdT_NKlike_unresolved，不作为常规NK，也不因细胞毒或CD8单个信号直接全部改作常规CD8。这里的处理是注释纠正与明确不确定性，不是断言外部样本没有NK。

GSE264586的TOX富集群4继续T_unresolved；TOX不能解决CD4/CD8身份，且仍无已解决的Treg群。其源对象tissue_status为Unknown，tissue_confirmed记录的是肝、肺、淋巴结等取材部位，不能直接等价于已确认肿瘤性质。本版严格同肿瘤外部验证不纳入该队列；原细胞、标签、时间点及全部排除原因保留，可在原始取材证据补齐后以新版本重新登记，不覆盖本冻结版。

GSE294109保留具有明确T骨架/CD8证据的低置信单患者集中群7，同时记录患者集中标记；其功能障碍、增殖、干扰素状态不作为移除CD8的理由。Treg群5保留。通过本轮身份、同样本、计数及覆盖门槛的CD8两端点和Treg/CD4，在各取材/方案层均为：

| 外部分层 | 合格患者 | 合格样本 |
|---|---:|---:|
| 未治疗手术 | 1 | 2 |
| 治疗后Axitinib/Pembrolizumab | 2 | 2 |
| 治疗后Ipilimumab/Nivolumab | 4 | 9 |
| 治疗后Lenvatinib/Pembrolizumab | 2 | 2 |

这些分层均不足10例，只作为探索性对应；单一未治疗患者不进行推断性检验。不得把多个取材区域当独立患者，也不自动把不同治疗层合并成同质验证队列。外部没有本项目冻结C4，Macrophage_SPP1不重命名为C4。本版无合格的常规NK外部端点。

## 五、下游接口已经固定，关联分析尚未启动

主要预测量为同文库成熟巨噬SPP1原始计数伪汇总log2(CPM+1)；C4比例及C4内SPP1作为独立次要预测量。功能读出预设为对应谱系内等权基因平均log1pCP10K；CD8细胞毒五基因与功能障碍八基因保存于FREEZE_DEFINITION.json，活化面板独立记录，不用单个PDCD1命名耗竭。

基因集、细胞资格与同患者汇总规则冻结；具体研究调整模型在后续统计方案中明确。功能、组成、外部验证和临床不混入同一个FDR池。CD8两个主要功能假设为同类家族，Treg比例独立；次要预测量及NK另按预设同类家族报告，不重写历史P/FDR。

本轮完成标准为可复核且固定的标签与名单，而不是保证预期免疫方向成立。未定细胞保留；后续发现问题须新版本记录，不能依据阳性结果悄悄改标签或门槛。

## 六、验证与权威文件

VALIDATION.json记录独立复核：全部发现/外部ID及来源字段、C4原SHA256、原始来源计数抽查、完整Treg群而非FOXP3阳性筛选、逐文库与逐患者资格重新统计、外部治疗/取材分层及无C4替代。计算通过不等同全部身份获得实验验证；当前身份是多基因与群体结构支持的生信注释。

正式目录freeze_v2；freeze_v1因外部元数据按错误列读取索引而在精确ID匹配处停止，FAILED_STATUS明确其非权威，未使用模糊匹配或产生正式结果。首次准备的字符串ID查找较慢，已改为索引查找后正常完成；中断日志及脚本保留。

- DISCOVERY_RECEIVER_CELL_FREEZE.tsv.gz：每细胞正式标签、来源、置信度和完整未定集合。
- DISCOVERY_UNIT_ANNOTATION.tsv：群/子单元的身份依据。
- DISCOVERY_ENDPOINT_LIBRARY_FREEZE.tsv、DISCOVERY_PATIENT_ENDPOINT_FREEZE.tsv：主要及C4内表达资格名单与全部排除原因。
- DISCOVERY_ENDPOINT_STUDY_COMPOSITION.tsv：端点实际研究构成。
- 两份外部RECEIVER_CELL_FREEZE、EXTERNAL_CLUSTER_ANNOTATION_FREEZE.tsv及外部PATIENT_ENDPOINT_FREEZE：原群、新标签、取材及治疗分层。
- FREEZE_DEFINITION.json、STATUS.json、FREEZE_SHA256.json：规则、完成状态及不可静默改变的文件身份。

项目权威进度仍只记服务器。前两部分、C4原数据和历史免疫统计均未覆盖。
'''
(O/'IMMUNE_RECEIVER_REFREEZE_CN.md').write_text(report)
now=datetime.datetime.now().astimezone().isoformat();source_paths=[R/'part2_rebuild_v1/09_downstream_crosswalk/C4_CANDIDATE_CELL_IDS.tsv.gz',R/'part2_rebuild_v1/08_regulons/input_v1/CELL_METADATA.tsv.gz',R/'phase210_immune_identity_pairing_check_20260919_v1/CELL_IDENTITY_AUDIT.tsv.gz']
(O/'SOURCE_SHA256.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths},indent=2))
# Seal formal artifacts after all final tables have been written.
seal={str(p.relative_to(F)):hashlib.sha256(p.read_bytes()).hexdigest() for p in F.rglob('*') if p.is_file()};(F/'FREEZE_SHA256.json').write_text(json.dumps(seal,indent=2))
status={'status':'COMPLETE_FROZEN_VALIDATED','completed_at':now,'authoritative_freeze':str(F),'report':str(O/'IMMUNE_RECEIVER_REFREEZE_CN.md'),'validation':str(O/'VALIDATION.json'),'new_association_tests':False,'C4_modified':False,'source_reports_modified':False};(O/'STATUS.json').write_text(json.dumps(status,indent=2))
backup=O/'root_record_backups';backup.mkdir(exist_ok=False)
for name in ['PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv','CODEX_HANDOFF.md','AGENTS.md']:shutil.copy2(R/name,backup/name)
note='''\n## 2026-09-19 Phase211：免疫接收端重新冻结完成\n用户授权先整理再重新冻结。权威phase211_immune_receiver_refreeze_20260919_v1/freeze_v2，报告IMMUNE_RECEIVER_REFREEZE_CN.md。86444接收细胞全部保留：CD8 31480、CD4常规13816、Treg3956、NK12520、gdT富集6562，其余未定/候选明确登记；没有按逐细胞FOXP3或功能评分挑细胞。原scVI表征内T/NK子图20群；混合主群3两个有多标志依据的分区分别归CD4/CD8，其他未解决群保留。成熟巨噬继承Phase154的22444全组织口径，C4原2674成员及SHA256不变。\n主要规则50巨噬、20CD8/NK、50CD4分母、未定身份<=20%、功能至少3基因各检出5细胞。发现队列CD8两端点各16患者/19库/4研究，NK与Treg各13患者/16库/2研究；C4内SPP1另有11/10患者。逐文库/患者名单及排除原因冻结。取代Phase209/210的存在/候选人数，不冒充已得关联。\n外部原NK群经原始计数确认具有明显TCR/CD3及gamma证据，单列未定，不冒充常规NK或CD8。GSE264586还缺肿瘤性质字段确认及Treg/TOX群身份，暂不进入严格肿瘤外部主验证；GSE294109通过资格的未治疗1例、三治疗层2/4/2例，分别探索性，不以多样本/治疗层扩充n。外部无C4替代。完整冻结及独立验证通过，无新免疫关联/通信/空间/图。后续从freeze_v2固定名单进入统计方案；任何标签变更须新版本，不按显著性改名单。\n'''
for name in ['PROJECT_STATUS.md','DECISIONS.md']:
 with (R/name).open('a') as f:f.write(note)
p=R/'CODEX_HANDOFF.md';p.write_text(note+'\n'+p.read_text())
with (R/'AGENTS.md').open('a') as f:f.write('\n## 2026-09-19 用户授权免疫接收端整理并重新冻结\n本轮已完成Phase211，正式接收端与患者资格以phase211_immune_receiver_refreeze_20260919_v1/freeze_v2为准。前两部分与C4成员不变。未定细胞保留；外部混合T/NK群不能继续直接当常规NK。后续不得按关联显著性改变标签或纳入名单；必要更改须新版本。当前仅冻结完成，免疫关联/通信/空间尚未启动。\n')
with (R/'RESULTS_MANIFEST.tsv').open('a') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n')
 for p in list(sorted(F.iterdir()))+[O/'IMMUNE_RECEIVER_REFREEZE_CN.md',O/'VALIDATION.json',O/'STATUS.json']:
  if p.is_file():w.writerow(['phase211_'+p.stem,'phase211','immune_receiver_freeze',str(p),now,'FROZEN_VALIDATED',str(O/'VALIDATION.json'),'Exact cell IDs and endpoint eligibility; no immune association tested'])
# Check seal after registration, without changing sealed contents.
assert all(hashlib.sha256((F/k).read_bytes()).hexdigest()==v for k,v in seal.items())
print('FINALIZED_PHASE211',len(seal),'sealed files',validation['n_checks'],'independent checks')
print(comp.to_string(index=False));print(summary[summary.source_tissue=='Tumor'].to_string(index=False))
