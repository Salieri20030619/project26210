from pathlib import Path
import json,csv,hashlib,datetime,shutil
import pandas as pd
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');O=R/'phase173_L1_MTA_evidence_integration_20260916_v1'
status=json.loads((O/'EVIDENCE_BUILD_STATUS.json').read_text());assert status['checks']==status['passed']
serial=pd.read_csv(O/'L1_SERIAL_MODEL_CONTEXT.tsv',sep='\t')
assert len(serial)==36 and ((serial.serial_low<=0)&(serial.serial_high>=0)).all()
biochem="""# L1生化衔接说明：MTA连接多胺合成与蛋氨酸补救循环

本文件是整合文稿的机制解释素材，不是新的完整报告。五条反应均已在2026-09-16核查Reactome并保存JSON。

| 环节 | 生化依据 | 本项目使用方式 |
|---|---|---|
| TYMP核苷入口 | TYMP催化胸苷/脱氧尿苷与磷酸参与的反应，生成相应碱基和脱氧核糖-1-磷酸 | 支持核苷代谢入口；该反应不直接产生MTA |
| AMD1供给 | SAM脱羧形成dcSAM | 解释SAM与多胺反应的物质联系，不声称本项目测定其酶活或通量 |
| SRM支路 | dcSAM与腐胺参与亚精胺合成，同时生成MTA | 为原报告SRM/多胺程序引出MTA提供具体反应依据 |
| SMS支路 | dcSAM与亚精胺参与精胺合成，同时生成MTA | 解释公共SMS—实测MTA关联的生化相容性；不要求第一部分SMS必须转录上升 |
| MTAP清除 | MTA磷酸解生成腺嘌呤和甲硫核糖-1-磷酸，进入蛋氨酸补救循环 | 清除背景独立描述，不预设缺失或普遍失活 |

关键推理是：第一部分观察到TYMP相关肿瘤状态伴随蛋氨酸/多胺相关程序，并有ODC1/SRM对应；这些程序所涉及的多胺反应产生MTA，因此MTA值得作为候选代谢介质进一步考察。该推理不等于TYMP直接催化MTA合成，也不证明转录变化已经造成MTA净积累或胞外暴露。

SRM与SMS并非必须在全部患者、细胞状态和轨迹上同时升高。将MTA列为候选无需抹去SMS轨迹反向；但也不能用这一生化关系把任何阴性解释成阳性。合成、清除、运输和区室均影响实际丰度，本轮没有直接估计这些通量。

来源：
- [TYMP反应](https://reactome.org/content/detail/R-HSA-112265)
- [AMD1反应](https://reactome.org/content/detail/R-HSA-351222)
- [SRM反应](https://reactome.org/content/detail/R-HSA-351215)
- [SMS反应](https://reactome.org/content/detail/R-HSA-351210)
- [MTAP反应](https://reactome.org/content/detail/R-HSA-1237160)
"""
(O/'L1_MTA_BIOCHEMICAL_LINKS.md').write_text(biochem,encoding='utf-8')
transition="""# L1可并入统一文稿的正文与衔接素材

## 使用说明
第一部分PDF是已完成结果的来源容器，最终可拆开重排。以下材料按结果、机制解释和衔接分别提供，不另起一份独立Methods—Results—Discussion报告。原结果、统计家族和关键反向证据保持；标题与编号待整篇编排时统一。

## 模块A：承接原代谢结果（与原2.6节合并，避免重复）
### TYMP相关肿瘤状态伴随蛋氨酸及多胺相关程序增强
在14位可配对患者中，TYMP阳性恶性细胞的半胱氨酸/蛋氨酸与精氨酸/脯氨酸通路分别有14位和13位同向升高，AUCell差值中位数分别为0.01019和0.01401，P分别为0.00110和0.00137，FDR均为0.00226。酶基因层面，ODC1与F1在18/20位患者中正向，中位rho=0.08108，P=0.000273、FDR=0.00327；SRM在14/19位患者中正向，中位rho=0.08945，P=0.01491、FDR=0.04262。这些结果将TYMP相关状态的代谢特征定位到蛋氨酸及部分多胺合成相关环节。

## 模块B：从上述发现引出MTA（接在代谢结果后的研究依据）
多胺合成提供了连接上述转录状态与MTA的具体生化依据：SRM和SMS催化的氨丙基转移反应分别参与亚精胺和精胺合成，同时生成MTA；MTA又通过MTAP参与的分解进入蛋氨酸补救循环。结合肿瘤细胞中已观察到的蛋氨酸/多胺相关程序，我们将MTA作为后续考察肿瘤代谢与髓系微环境联系的候选代谢物。
[SRM反应](https://reactome.org/content/detail/R-HSA-351215)、[SMS反应](https://reactome.org/content/detail/R-HSA-351210)、[MTAP反应](https://reactome.org/content/detail/R-HSA-1237160)。

## 模块C：已有公共实测代谢组的分段支持（可接候选引出后）
### 公共实测代谢组显示MTA与多胺合成及清除背景相联系
四个配对组织RNA—代谢组队列共255位肿瘤患者的既有分析显示，控制TYMP、MTAP及髓系RNA代理后，SMS与实测MTA仍在四队列同向，随机效应合并rho=0.24993，P=0.000578、原六基因家族FDR=0.00289。MTAP与MTA的未调整关联四队列均负，随机效应合并rho=-0.25044，P=0.00957、对应六基因家族FDR=0.02870；加入髓系代理后该关联减弱，rho=-0.21583，P=0.06011、FDR=0.12023。上述结果分别支持MTA与合成、清除背景的联系，效应随模型而异。
组织层TYMP与MTA的直接未调整关联较弱，随机效应合并rho=0.05843、P=0.58509；因此，公共数据对候选选择的支持主要来自代谢网络及分段关联，尚未确认TYMP对组织MTA的完整作用路径。该原始检验的历史固定效应FDR为0.72422，Phase160相应随机效应背景单项FDR为0.58509，两者统计口径不同，完整并列于源表，不能互相替换。

## 模块D：进入髓系部分前的研究问题（研究转折，不作为完成结果）
基于上述代谢特征和MTA的候选定位，下一步将从单核—巨噬细胞的组成、状态及功能出发，检视MTA相关作用是否与SPP1富集巨噬表型相衔接。在此基础上，再分别评估其免疫功能联系和空间组织联系。

## 模块E：并入统一Discussion的边界与关键反向结果
现有单细胞证据表征转录状态，scFEA未显式输出MTA；TYMP的经典核苷反应也不直接生成MTA。公共组织MTA不能直接代表肿瘤细胞胞外暴露，整体组织SPP1不能归为巨噬来源。既有MTA—组织SPP1关联为名义正向（rho=0.13894、P=0.02926、历史FDR=0.09755），而在进一步条件模型及串联模型中支持减弱；36个串联模型的间接乘积区间均跨零。Phase170空间RNA代理桥接同样未获得稳健支持。这些结果限定了完整线路的解释，但不抹去已建立的肿瘤代谢状态及SMS—MTA分段联系。

原第一部分SMS、PCYT1A的轨迹反向结果、MTAP未一致下降以及不同轨迹不一致均继续保留。不能把MTA候选表述成已证明“爆发性积累”“跨界释放”或已诱导SPP1。

## 模块F：实验验证位置
客户现有shTYMP代谢图在最终验证段使用，重述其报告的MTA及相关代谢物变化，不作为模块A的发现起点。不重复计算客户实验，也不让该图承担胞外释放、SPP1响应、免疫抑制和空间屏障等未展示的终点。原第一部分3.2节及表7应随整篇结构移至实验对应/综合讨论，而不是在每个章节重复放置。

## 证据编号
模块A对应L1_TYMP_METABOLIC_EVIDENCE.tsv的auc_pair/enzyme_f1；模块C对应L1_PUBLIC_MTA_EVIDENCE.tsv的Phase162/SMS/b_adjust_MTAP/myeloid、MTAP/XY/raw与myeloid及Phase160/MTA/a/raw。正文未使用新检验或新的更小FDR家族。
"""
(O/'L1_TRANSITION_TEXT_CN.md').write_text(transition,encoding='utf-8')
assembly="""# 统一文稿拆分编排清单

| 原报告素材 | 后续整合位置 | 处理方式 |
|---|---|---|
| 原标题、摘要、结尾 | 全文最终标题/摘要/结论 | 待所有模块齐备后统一改写，不把第一部分摘要当全文结论 |
| 1.1—1.8设计、数据及方法 | 统一Methods | 与髓系、免疫、空间方法分节合并；原参数和统计家族保留，不把全部检验混校正 |
| 2.1恶性身份及CNV背景 | 前部肿瘤图谱与状态发现 | 保留必要证据；QC完整明细可入补充，主要冲突/界限不隐藏 |
| 2.2 F1/TYMP及2.4/2.5候选背景 | 肿瘤状态与调控模块 | 根据整篇主线重排；不强制全部候选进入MTA或SPP1路径 |
| 2.6代谢结果 | 肿瘤代谢结果模块 | 作为本轮模块A落点，与MTA选择依据顺接 |
| 2.7轨迹结果 | 状态动态支持与补充 | 保留方法差异、SMS/PCYT1A反向，不强写统一末端 |
| 原3.1讨论 | 对应结果的综合解释 | 拆入全文Discussion，不复制成独立第一部分Discussion |
| 原3.2及表7客户实验对应 | 最终实验验证/综合讨论 | 不作为生信发现起点，不替代未展示的下游实验 |
| 原3.3/3.4统计与测量边界 | 统一Discussion及技术附件 | 保留关键边界，避免各结果后重复防御性说明 |
| 原统计附表 | 全文补充表及溯源档案 | 保留原精度、校正家族和端点身份，最后统一编号 |
| 原Fig1—Fig7 | 全文图件素材库 | 原图号暂作溯源ID，最终整体编号；本轮不修改PDF、重绘或拼图 |
| 本轮L1材料 | 代谢结果到MTA候选的衔接 | 模块化插入，不另起独立报告 |
| 后续免疫分支 | 独立免疫功能结果模块 | 状态、组成与通信结局 |
| 后续空间分支 | 独立空间组织结果模块 | 位置、邻近、边界、CD8分布结局 |

已完成数据与统计作为固定依据，不等于固定文稿结构。最终文章只保留一套协调的研究问题、Methods、Results、Discussion和图表编号。
"""
(O/'MANUSCRIPT_ASSEMBLY_MAP_CN.md').write_text(assembly,encoding='utf-8')
summary="""# 第一步完成：从既有肿瘤代谢结果合理引出MTA

## 已完成
- L1-01：提取原报告83条相关统计，核对十二酶报告数值与源表、AUCell原统计。
- L1-02：核查并归档5条官方生化反应，区分TYMP核苷入口、AMD1供给、SRM/SMS生成MTA和MTAP清除。
- L1-03：整理公共实测代谢组91条既有结果，核对RC18/RC20/CPTAC/CPTAC_val分别97/52/50/56患者及MTA特征；保留原P/FDR、历史家族和主要阴性。
- L1-04：形成可并入统一文稿的承接段、MTA选择依据、公共支持、下一层转折及Discussion素材。
- 另形成第一部分拆分编排表、详细计划v2，明确统一文稿而非独立报告串联。

## 结论
MTA可由第一部分已建立的蛋氨酸/多胺代谢状态、ODC1/SRM对应和明确生化反应引出；公共SMS—实测MTA与MTAP—MTA结果提供合成/清除背景的分段支持。当前不将其表述为已证明TYMP直接生成MTA、MTA胞外积累或完整中介。
下一层需要回答的是MTA相关作用与本项目单核/巨噬及SPP1状态之间的具体对应，而不是再次检验第一部分是否成立。

## 完整性
本轮复用原统计，不生成新的P/FDR，不重跑第一部分；所有主要阴性、条件关联衰减及Phase170空间代理不足均保留。第二层及免疫/空间新分析尚未启动，图件未修改。

## 阅读顺序
1. L1_TRANSITION_TEXT_CN.md：可直接拆分并入全文的文字模块。
2. MANUSCRIPT_ASSEMBLY_MAP_CN.md：第一部分与后续统一编排。
3. L1_MTA_BIOCHEMICAL_LINKS.md：生化理由及官方来源。
4. L1_TYMP_METABOLIC_EVIDENCE.tsv / L1_PUBLIC_MTA_EVIDENCE.tsv：完整证据账本。
5. DETAILED_ANALYSIS_PLAN_v2_CN.md：更新后的执行与编排口径。
"""
(O/'L1_COMPLETION_CN.md').write_text(summary,encoding='utf-8')
plan=O/'DETAILED_ANALYSIS_PLAN_v2_CN.md'
s=plan.read_text().replace('第一步证据整理与衔接已启动','第一步证据整理与衔接已完成').replace('本次启动第一步L1-01至L1-04','本次完成第一步L1-01至L1-04').replace('**本次状态：详细计划已形成，新增分析未启动。**','**本次状态：L1-01至L1-04已完成；无新统计检验，第二层及免疫/空间分支新增分析未启动。**')
plan.write_text(s,encoding='utf-8')
now=datetime.datetime.now().astimezone().isoformat()
status.update(status='PASS_L1_EVIDENCE_INTEGRATION_COMPLETE',completed_at=now,manuscript_structure='single_integrated_manuscript_with_reassembled_existing_part1',serial_models_retained=36,serial_CI_all_cross_zero=True)
(O/'STATUS.json').write_text(json.dumps(status,ensure_ascii=False,indent=2))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
pd.DataFrame([dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(O.rglob('*')) if p.is_file() and 'authority_backups' not in p.parts and p.name not in ['ARTIFACT_HASHES.tsv','FINAL_STATE_VALIDATION.json']]).to_csv(O/'ARTIFACT_HASHES.tsv',sep='\t',index=False)
B=O/'authority_backups';B.mkdir(exist_ok=False)
names=['PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv','CODEX_HANDOFF.md','AGENTS.md']
for n in names:shutil.copy2(R/n,B/n)
notice=f"""
## Phase173 第一步L1完成，文稿按整体拆分整合（{now}）
用户明确第一部分现为完整报告结构，但最终拆开加入后续统一文稿；固定结果而非章节/图号，不是多个独立报告串联。已完成L1-01至L1-04：原报告83条统计、公共MTA91条既有结果、5条官方反应；255患者实测MTA身份核对。141检查PASS，无新P/FDR或模型。
MTA由蛋氨酸/多胺背景及SRM/SMS生成反应合理引出，公共SMS—MTA与MTAP—MTA提供分段支持；TYMP—组织MTA弱、原串联36模型区间跨零与Phase170空间不足保留。不称直接生成/胞外积累/完整中介。
模块正文：{O/'L1_TRANSITION_TEXT_CN.md'}
编排清单：{O/'MANUSCRIPT_ASSEMBLY_MAP_CN.md'}
计划v2：{plan}
后续入口为L2直接作用资料—单核/巨噬/SPP1状态对照；免疫功能和空间组织继续两条独立分支。原PDF未改、无新图、第二层未执行。
"""
for n in ['PROJECT_STATUS.md','DECISIONS.md','CODEX_HANDOFF.md']:
 p=R/n;t=p.with_name(n+'.phase173.tmp');t.write_text(notice+p.read_text(),encoding='utf-8');t.replace(p)
p=R/'AGENTS.md'
rule="""

## 2026-09-16 用户确认全文整合方式与第一步执行
第一部分现有PDF的完整报告结构仅为阶段性交付。最终需要拆分、重排并入后续统一文稿，共享Methods、Discussion及最终图表编号；固定的是已完成数据/结论，不是独立报告结构。不得再按“完整第一部分报告+独立后续报告”机械串联。用户已授权第一步L1，Phase173完成证据整理及模块衔接；后续沿计划v2执行，免疫与空间终点分开。
"""
t=p.with_name('AGENTS.md.phase173.tmp');t.write_text(p.read_text()+rule,encoding='utf-8');t.replace(p)
p=R/'RESULTS_MANIFEST.tsv'
with p.open(newline='') as f:rows=list(csv.reader(f,delimiter='\t'))
assert len(rows[0])==8 and not any(x and x[0].startswith('phase173_') for x in rows[1:])
reg=['L1_COMPLETION_CN.md','L1_TRANSITION_TEXT_CN.md','L1_MTA_BIOCHEMICAL_LINKS.md','L1_TYMP_METABOLIC_EVIDENCE.tsv','L1_PUBLIC_MTA_EVIDENCE.tsv','L1_PUBLIC_SAMPLE_MEASUREMENT_AUDIT.tsv','L1_BIOCHEMICAL_LINKS.tsv','MANUSCRIPT_ASSEMBLY_MAP_CN.md','DETAILED_ANALYSIS_PLAN_v2_CN.md','SOURCE_MANIFEST.tsv','VALIDATION.tsv','STATUS.json','build_L1_evidence.py','finalize_L1.py','ARTIFACT_HASHES.tsv']
for n in reg:rows.append(['phase173_'+n,'Phase173','L1_evidence_integration',str(O/n),now,'PASS','141 source checks; no new statistical tests','sha256='+sha(O/n)])
t=p.with_name(p.name+'.phase173.tmp')
with t.open('w',newline='') as f:csv.writer(f,delimiter='\t',lineterminator='\n').writerows(rows)
t.replace(p)
gate=dict(stage='Phase173',registered=len(reg),all_paths_exist=all((O/n).exists() for n in reg),four_records_updated=all((R/n).read_text().startswith(notice) for n in ['PROJECT_STATUS.md','DECISIONS.md','CODEX_HANDOFF.md']),AGENTS_updated=True)
assert gate['all_paths_exist'] and gate['four_records_updated']
(O/'FINAL_STATE_VALIDATION.json').write_text(json.dumps(gate,indent=2));print(json.dumps(status,ensure_ascii=False));print(json.dumps(gate))
