from pathlib import Path
import json,datetime,hashlib,shutil
import pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');R=P/'part2_rebuild_v1';O=P/'phase152_TYMP_SPP1_direct_bridge_plan_20260915_v1'
z=pd.read_csv(O/'SAME_LIBRARY_FEASIBILITY.tsv',sep='\t');e=z[z.eligible_same_library]
assert len(e)==16 and e.curated_patient_id.nunique()==16 and e.sample_id.is_unique
text="""# 下一阶段直接补强肿瘤TYMP—巨噬SPP1联系

## 主任务调整

三TF的当前证据已完成阶段性定位：保留为髓系状态相关调控背景，不将JUN/JUNB/CEBPD设为连接TYMP和SPP1的必经中介。下一阶段围绕肿瘤细胞TYMP与巨噬SPP1的直接对应、细胞组成解释、空间邻近及独立复现组织。

已有17患者结果为候选线索：肿瘤TYMP RNA—巨噬SPP1表达rho=0.507353、P=0.037633、FDR=0.188164；研究调整后P=0.4324。保持该结果原统计家族，不因本轮缩小关注范围重算旧FDR或升级证据。

## 第一优先：同患者、同肿瘤文库的直接对应

已从冻结H5AD读取元数据，与当前成熟巨噬细胞表精确对齐。使用原恶性细胞资格规则、Tumor来源，单文库至少20恶性细胞和50成熟巨噬细胞，共16文库、16名独立患者满足门槛：
GSE156632为4名、GSE159115为1名、GSE178481为1名、GSE207493为10名。细胞数不是样本量。其余细胞/患者保留资格表，不能计入该主比较。这里确认共享sample_id文库；其他跨文库同患者信息仅作单独敏感性。

拟执行端点在计算前固定：
1. 主要目标：同文库恶性细胞TYMP伪汇总表达，与同文库全部成熟巨噬SPP1伪汇总表达的关联。以log2(CPM+1)描述连续表达，避免二分阈值。
2. 分解端点：同文库巨噬中SPP1检出比例，用于区分整体表达与表达细胞比例。
3. 状态端点：冻结C4占成熟巨噬比例；不将C4等同全部SPP1阳性细胞，不重新调整C4成员。

三项关联同属直接桥接家族，完整P和BH FDR共同保存；正、负方向都保留。原始跨患者Spearman为整体描述，研究内秩调整/置换用于评估研究构成影响，逐研究和留一患者/留一研究列明方向。GSE159115及GSE178481在此各n=1，不能做研究内推断；GSE207493属于发现资料，不能另称独立验证集。预先固定3项，避免从大量指标中只挑最显著的一个。

该步骤预期解决：此前正关联是否在同一肿瘤样本的两个细胞区室间存在，是否由研究构成、SPP1检出比例或C4比例主导。它仍为观察性关联，不自动产生TYMP诱导SPP1的因果结论。

## 第二优先：具有细胞来源约束的空间邻近检验

先审计可复用的TJ-RCC与GSE175540现有表达、映射、坐标和切片身份。原GSE175540组成调整TYMP—SPP1程序关系两种制备均不支持，原结果保留；不能重复同端点换展示方式当作新验证。

有新增判别力的目标是“恶性区域的TYMP暴露，与相邻巨噬SPP1表达/程序的关系”，同时考虑恶性/巨噬组成和测序深度，采用切片内空间块或空间保持置换，最后以独立患者/已核实切片身份汇总。若现有数据不能区分TYMP的细胞来源，必须保留为区域暴露而非肿瘤细胞特异表达。先固定邻域、读出和检验再计算。

复用既有cell2location结果；旧SPP1-state丰度保持原名。C4若作为空间读出，需要先通过新参考可识别性审查，不能直接把旧后验改名。本计划不启动新的cell2location训练。

## 第三优先：独立资料复现同一个直接联系

优先核对现有外部单细胞资料是否同时包含合格恶性细胞和巨噬，以及可匹配的患者/组织样本。免疫富集或只有巨噬数据的队列不能直接检验肿瘤细胞TYMP→巨噬SPP1。GSE207493与主发现数据重叠，不作独立复现。

若具备两区室输入，复现第一步冻结的同样本端点与方向，队列分开报告，不混入发现队列提高n。若仅有bulk，可支持组织层共变，但不能据此指定到细胞类型；bulk或细胞系MTA结果也不能代替巨噬SPP1证据。

## 机制背景与报告组织

肿瘤TYMP—代谢背景—SPP1状态可以继续作为目标框架，但每一段使用其对应证据。客户湿实验由客户负责，继续按已提供方向衔接，不要求重新分析原始湿实验，也不把计算结果包装成干预验证。

先以直接跨细胞区室关联作为核心，SPP1群身份和功能作为承接，空间/独立资料作为复现；三TF不再是该主线的必需节点。只有在明确主关系后，才对支持的候选机制作定向延伸。

## 本轮实际完成

已完成方向调整、已有证据核对及同文库配对元数据审计；没有新算上述TYMP—SPP1三端点，没有启动空间或独立队列模型，没有新图。全部记录在服务器版本目录。输入与资格表为SAME_LIBRARY_FEASIBILITY.tsv和FEASIBILITY_STATUS.json。下一执行入口为第一优先的16患者同文库直接桥接。
"""
(O/'DIRECT_TYMP_SPP1_NEXT_STAGE_PLAN_CN.md').write_text(text,encoding='utf-8')
plan={'status':'FEASIBILITY_AND_PLAN_COMPLETE','next_step':'same_library_direct_TYMP_SPP1_bridge','eligible_patients':16,'eligible_libraries':16,'primary_cell_scope':'frozen malignant cells and all mature macrophage candidates; tumour only; no C1/Monocyte as macrophages','endpoints':['malignant_TYMP_pseudobulk_vs_macrophage_SPP1_pseudobulk','malignant_TYMP_pseudobulk_vs_macrophage_SPP1_detection_fraction','malignant_TYMP_pseudobulk_vs_frozen_C4_fraction'],'fdr_family':'three related direct bridge endpoints; raw and study-adjusted separate; retain prior Phase151 FDR unchanged','threeTF_mandatory_mediator':False,'new_hypothesis_results_computed':False,'figures_created':False,'model_rerun':False,'time_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(O/'STATUS.json').write_text(json.dumps(plan,indent=2))
shutil.copy2(R/'scripts/39_direct_bridge_feasibility_20260915_v1.py',O/'01_feasibility.py')
B=O/'state_before_update';B.mkdir(exist_ok=False)
for q in [P/'CODEX_HANDOFF.md',P/'PROJECT_STATUS.md',P/'DECISIONS.md',P/'RESULTS_MANIFEST.tsv',R/'STATUS.json']:shutil.copy2(q,B/q.name)
note="""## 2026-09-15 Phase152 主任务转向直接补强TYMP—SPP1联系

用户明确下一步不再围绕原三个TF，转向补强TYMP与SPP1关系。三TF阶段性收束，不作必经中介。权威计划phase152_TYMP_SPP1_direct_bridge_plan_20260915_v1/DIRECT_TYMP_SPP1_NEXT_STAGE_PLAN_CN.md。
已完成同患者同Tumor文库元数据盘点：原冻结恶性资格≥20细胞、当前成熟巨噬≥50细胞，16文库/16患者（GSE156632 4、GSE159115 1、GSE178481 1、GSE207493 10）可做直接配对。后续从该配对集合计算恶性TYMP与巨噬SPP1伪汇总表达、SPP1检出比例、冻结C4比例三端点，三项家族校正，分研究与留一敏感性，保留Phase151名义阳性及原FDR。
顺序：同样本直接桥接→现有空间来源/邻近性端点审查→外部双细胞区室队列复现。旧GSE175540 TYMP-SPP1程序阴性保留；旧SPP1后验不改名C4，不启动新cell2location。用户“先审内容后图”仍生效。当前仅可行性与计划完成，尚未产生本轮三端点新统计；下一入口是16患者配对分析。
"""
for fn in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
 q=P/fn;q.write_text(note+'\n'+q.read_text(),encoding='utf-8')
q=R/'STATUS.json';d=json.loads(q.read_text());d['remaining_analysis'].update(current_step='direct_TYMP_SPP1_bridge_prepared',next_plan=str(O/'DIRECT_TYMP_SPP1_NEXT_STAGE_PLAN_CN.md'),threeTF_mandatory_mediator=False);q.write_text(json.dumps(d,ensure_ascii=False,indent=2))
with (P/'RESULTS_MANIFEST.tsv').open('a') as f:f.write('\t'.join(['phase152_direct_bridge_plan','phase152','plan_and_feasibility',str(O/'DIRECT_TYMP_SPP1_NEXT_STAGE_PLAN_CN.md'),plan['time_utc'],'PLAN_READY','16 same-library patients; no new effect computed','No figures'])+'\n')
print(json.dumps(plan));print(hashlib.sha256((O/'DIRECT_TYMP_SPP1_NEXT_STAGE_PLAN_CN.md').read_bytes()).hexdigest())

