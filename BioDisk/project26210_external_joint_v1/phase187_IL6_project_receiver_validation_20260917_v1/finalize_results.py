from pathlib import Path
import pandas as pd,json,shutil,hashlib,datetime
root=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=root/'phase187_IL6_project_receiver_validation_20260917_v1'
def md(name,cols):
 a=pd.read_csv(d/name,sep='\t')[cols]
 lines=['| '+' | '.join(cols)+' |','|'+'|'.join(['---']*len(cols))+'|']
 for row in a.itertuples(index=False,name=None):
  lines.append('| '+' | '.join(f'{x:.6g}' if isinstance(x,float) else str(x) for x in row)+' |')
 return '\n'.join(lines)
report="""# Phase187：IL6接收端与实测MTA对应检验

## 结论

**IL6与STAT3转录代理的正向联系得到本项目及公共队列支持；当前数据不支持IL6已连接MTA与SPP1。**

本轮已连续完成来源、受体、整体巨噬患者内关联和公共实测MTA对应核查，无需另行选择后才计算。成纤维中转仍保留为有外部依据的候选；不将样本不足写成机制不存在，也不把已成立的IL6—STAT3联系抹去。

## Results 1：成纤维来源检验不可估计，巨噬受体具备RNA测量基础

冻结75文库258120细胞中，肿瘤Fibroblast标签仅39细胞、1患者，IL6原始计数为0。没有足够独立生物学重复，成纤维IL6—巨噬SPP1/STAT3跨患者关联不可估计。未重分群、未改标签，也未把Pericyte/SMC或内皮细胞替代成纤维来源。

整体巨噬受体覆盖：34个达到至少20巨噬细胞的患者，IL6R及IL6ST各30个患者至少10个细胞检出；共同检出19个患者达到同样门槛。此为转录测量资格，不代表已证明受体蛋白、IL6暴露或信号传递。

IL6存在多类来源：按每类至少20细胞，其他成熟巨噬19/34、内皮16/34、保守恶性5/20患者达到至少10个IL6检出细胞。来源概览全表保留；不同来源的分母不同，不用这些比例声称全组织分泌贡献。

## Results 2：整体巨噬IL6—SPP1没有稳定正向关联

主资格为每患者至少50个成熟肿瘤巨噬，IL6和SPP1分别至少10个检出细胞。19患者通过；主深度调整仅10/19正向，中位rho接近0。加入C4组成或检出基因数调整后方向仍不稳定。

"""+md('MACRO_SUMMARY.tsv',['adjustment','n_patients','positive','median_rho','CI_low','CI_high','P'])+"""

这里检验的是同一巨噬细胞群中的IL6转录与SPP1关系，不是胞外IL6蛋白的作用；阴性不能否定外部IL6干预实验，但不能据此声称本项目IL6—SPP1成立。

## Results 3：IL6—STAT3代理跨患者和跨资源同向

IL6不在两套冻结STAT3靶集中，因此可以使用原评分，未为追求阳性调整靶基因。IL6ST在两套靶集中，未将IL6ST—评分相关作为受体功能证据。

19患者IL6—STAT3代理在主CollecTRI和DoRothEA敏感性中均19/19正向；控制C4或检出基因数方向保留。

"""+md('IL6_STAT3_MACRO_SUMMARY.tsv',['resource','adjustment','n','positive','median_rho','P'])+"""

公共4队列255例中，调整实测MTA、MTAP及髓系背景后，IL6—STAT3代理同样4/4正向，增加缺氧调整后仍保留。

"""+md('IL6_STAT3_PUBLIC_META.tsv',['resource','adjustment','n','rho','CI_low','CI_high','P','positive_cohorts'])+"""

这支持IL6与STAT3相关转录状态联系，不是pSTAT3直接测量，也不证明IL6位于MTA下游或STAT3是必经中介。已有STAT3—SPP1患者内联系保留，不能凭两段共同相关自动推断IL6—SPP1因果关系。

## Results 4：公共实测MTA—IL6不支持稳定正向联系

预设两条边：A为MTA—IL6调整MTAP及髓系背景；B为IL6—SPP1同时调整MTA、MTAP及髓系背景。主家族A/B两项一起BH校正；缺氧敏感性另作同类家族。

"""+md('PUBLIC_META.tsv',['adjustment','edge','n','rho','CI_low','CI_high','P','FDR','positive_cohorts'])+"""

主MTA—IL6仅1/4队列正向；IL6—SPP1主合并为弱正向，但P=0.0652、FDR=0.1303，不能称名义显著。增加缺氧调整后IL6—SPP1接近0。完整逐队列、逐研究剔除和患者数据保留，不只展示合并结果。未计算或宣称正向中介乘积。

## Methods

数据与注释：继承Phase182冻结全细胞10基因计数缓存及细胞索引；成熟巨噬来自冻结精细髓系标签，排除单核candidate1/Monocyte，C4保持原成员。缓存行与258120细胞索引完全一致。保留患者和研究身份，多文库在患者内汇总进行患者层面的效应推断。

单细胞：基因表达log1p(CP10K)，秩相关残差调整深度，另加检出基因数或C4指示变量作敏感性；数值12位小数后排名，防止浮点并列秩不一致。统计单位为独立患者，不是细胞。患者rho采用双侧Wilcoxon近似检验加连续性校正；主IL6—SPP1的中位效应区间为5000次患者重抽样百分位区间。每个独立预设主要终点不作FDR；资源/协变量敏感性明确标记，不挑其中最小P替代主结果。

STAT3：复用Phase182现有两套zmean评分及其原定义，IL6不参与评分；38个患者×资源的行数与SPP1缓存值核验通过。未将评分变化解释为STAT3磷酸化变化。

公共队列：RC18、RC20、CPTAC、CPTAC_val分别97、52、50、56例，合计255例。RNA采用原log2(TPM+1)，依据MasterMapping进行患者匹配；四队列TYMP/SPP1逐样本均复现既有冻结值。未将跨队列患者与单细胞患者拼成同患者中介。各队列偏Spearman相关、t近似双侧P；Fisher z固定效应汇总权重n-k-3，k为协变量数。完整P及家族内FDR保留。IL6—STAT3为独立接收端终点，主CollecTRI、DoRothEA敏感性，未并入MTA两边的校正池。

验证：32项公共关联采用QR独立实现复算全部通过；另核对来源资格、输入SHA256、255例映射和38项评分对齐。初版脚本因公共首列名称不同未完成，修订版analyze_v2.py读取实际首列并成功复现既有患者值。旧脚本保留作为过程记录。

## Discussion与全文处理

1. 可以在对应Results保留“IL6与STAT3相关转录状态正向关联”，包括患者一致性、效应、P和敏感性。
2. Phase186“MTA使成纤维IL6 RNA小幅上升”是外部细胞系探索性观察；外部IL6中和降低髓系SPP1为机制依据。二者均保留，不能写成在本项目同一体系已证明。
3. 本项目成纤维缺少重复、巨噬IL6—SPP1阴性、公共MTA—IL6阴性必须共同呈现。现阶段不把IL6插为全文必经中介。
4. 第一部分TYMP相关代谢发现继续作为起点，客户已有代谢图承担后续实验对应；不能为圆链条倒置起点或补写客户未做的实验。
5. 本轮已到该候选在现有数据中的可判断边界：保留IL6—STAT3阳性与MTA中转假说，停止用当前队列反复更换评分/细胞标签追求闭环。只有新的匹配来源数据或直接MTA上清—IL6阻断—髓系SPP1读出才会改变该中转结论。免疫与空间仍分别保留为下游终点，未启动被暂缓的模型。

## 主要外部依据

- Phase186完整结果与原始数据资格：../phase186_MTA_fibroblast_panel_20260916_v1/MTA_FIBROBLAST_PANEL_RESULTS_CN.md
- 基质IL6清除与髓系SPP1：[Nature Communications 2024](https://doi.org/10.1038/s41467-024-50438-2)。
- IL6时序及OPN反馈：[Oncogenesis 2026](https://doi.org/10.1038/s41389-026-00596-3)。
"""
(d/'IL6_PROJECT_VALIDATION_RESULTS_CN.md').write_text(report)
now=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()
backup=d/'project_record_backups';backup.mkdir(exist_ok=True)
note=f"""## Phase187 IL6接收端与实测MTA检验完成（{now}）
本轮来源/受体/患者内/公共队列/STAT3无循环检验均已完成，无新图。
冻75库Fibroblast肿瘤仅39细胞1患者、IL6计数0，跨患者中转不可估计；未换成Pericyte/SMC。巨噬IL6R/IL6ST各30患者>=10细胞，共同19患者。
整体巨噬19患者IL6—SPP1深度10正，中位rho0.000944/P0.7323，C4/基因数敏感性无稳定联系。IL6不在两套STAT3靶集中：IL6—代理19/19正，两资源及协变量同向，主rho0.077535/P0.000143。
公共255例：MTA—IL6主rho-0.09149/P0.1596/FDR0.1596；IL6—SPP1调整MTA等rho0.12074/P0.06517/FDR0.13035，缺氧后接近0。IL6—STAT3主rho0.52782/P4.52e-19，4/4正，缺氧后0.34403/P6.53e-8。
32项QR、4队列映射复现、38评分对齐通过。IL6ST参与原评分，禁止其与评分相关证明受体功能。旧结果全保留。
报告：{d}/IL6_PROJECT_VALIDATION_RESULTS_CN.md
结论：IL6—STAT3状态证据可保留；当前数据不支持IL6打通MTA—SPP1，不将跨体系两段拼成确认中介。现有数据此分支检验结束，不反复改评分找阳性。新来源数据或直接MTA上清/阻断才改变中转结论。免疫/空间仍分别暂缓，无活跃任务。
"""
for n in ['PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv','CODEX_HANDOFF.md']:
 p=root/n;shutil.copy2(p,backup/n)
 if n=='RESULTS_MANIFEST.tsv':
  with open(p,'a') as f:f.write('\n'+'\t'.join(['phase187_IL6_receiver','phase187','report',str(d/'IL6_PROJECT_VALIDATION_RESULTS_CN.md'),now,'PASS','32 QR; 38 score alignment; 255 mapping','fibroblast source not estimable; IL6-STAT3 positive; MTA-IL6 unsupported'])+'\n')
 else:p.write_text(note+'\n'+p.read_text())
files=['IL6_PROJECT_VALIDATION_RESULTS_CN.md','MACRO_SUMMARY.tsv','MACRO_PATIENT_EFFECTS.tsv','MACRO_ELIGIBILITY.tsv','PUBLIC_META.tsv','PUBLIC_COHORT_EDGES.tsv','PUBLIC_LEAVE_COHORT.tsv','IL6_STAT3_MACRO_SUMMARY.tsv','IL6_STAT3_PUBLIC_META.tsv','IL6_SOURCE_SUMMARY.tsv','IL6_RECEIVER_COVERAGE.tsv','FIBROBLAST_SOURCE_GATE.json','INDEPENDENT_QR_VALIDATION.tsv','STAT3_ALIGNMENT_VALIDATION.json']
dest=d/'delivery';dest.mkdir(exist_ok=True)
for n in files:shutil.copy2(d/n,dest/n)
ha={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in dest.iterdir()}
(dest/'SHA256.json').write_text(json.dumps(ha,indent=2))
print('REPORT',len(report),'DELIVERY',len(files))
