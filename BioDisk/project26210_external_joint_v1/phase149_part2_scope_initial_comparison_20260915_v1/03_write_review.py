from pathlib import Path
import json,hashlib,datetime,shutil,csv
import pandas as pd,numpy as np
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');R=P/'part2_rebuild_v1';O=P/'phase149_part2_scope_initial_comparison_20260915_v1';E=O/'evidence'
def read(name):return pd.read_csv(E/name,sep='\t')
def table(df):
 def fmt(v):
  if isinstance(v,float):return 'NA' if np.isnan(v) else f'{v:.6g}'
  return str(v).replace('|','/').replace('\n',' ')
 return '| '+' | '.join(map(str,df.columns))+' |\n| '+' | '.join(['---']*len(df.columns))+' |\n'+'\n'.join('| '+' | '.join(fmt(v) for v in row)+' |' for row in df.itertuples(index=False,name=None))
rows=[
 ['P2-1','提取Monocyte/Macrophage并细分巨噬','计算完成；发现定位为探索性','33,098轨迹细胞/39患者/70有效文库；主数据范围仍75文库。12巨噬候选、22,444细胞；C4冻结2,674细胞。','已完成7张发现阶段图；不能把17候选称17个稳定生物亚型。'],
 ['P2-2','Monocle3轨迹分析','计算完成；支持范围有限','5研究×主/去应激增殖两种模型，共10模型；根敏感性、DPT、组织限定均完成。','主模型仅442/2674 C4有有效拟时序；C4与单核可比患者n=8。'],
 ['P2-3','pySCENIC网络及JUN/JUNB/CEBPD活性','计算及患者统计完成','GRNBoost2、双motif数据库剪枝、AUCell完成；235 regulons、67,949预测边、33,098细胞。','三TF均得到可评估regulon；并未出现三者一致级联上扬。'],
 ['P2-4','Monocyte到TYMP+/SPP1+TAM轨迹流线图','交付未完成；按后续指令暂缓','模型、坐标与伪时间已保存，可用于后续结果导向呈现。','9月14日用户明确禁止新图；未生成本轮正式轨迹流线图。'],
 ['P2-5','Regulon瀑布热图','交付未完成；按后续指令暂缓','全细胞AUCell和动态统计已保存。','新瀑布热图及其最终平滑/分箱展示矩阵未交付；不能画成预设一致上扬。'],
 ['P2-6','最终分化为TYMP+/SPP1+促癌状态的结论','原假设未被完整支持；不是缺跑模型','SPP1富集和C4晚于单核的部分排序得到支持；TYMP相对富集、普遍终末位置和因果转化未建立。','当前可写SPP1富集的C4转录状态；不能将表达检出等同TYMP驱动。'],
 ['P2-7','JUN/JUNB/CEBPD级联上扬','当前数据不支持这一整体表述','三TF在C4相对其他巨噬均降低；轨迹JUN正、JUNB依赖方法、CEBPD无主模型支持。','保留实际方向；不把达到方案预想方向作为计算完成标准。']
]
coverage=pd.DataFrame(rows,columns=['编号','原方案要求','完成判定','证据','保留事项']);coverage.to_csv(O/'ORIGINAL_PLAN_COVERAGE.tsv',sep='\t',index=False)
extra=[
 ['功能定位','完成','患者程序/ORA、排除SPP1/TYMP版本、研究方向及原始pseudobulk核验'],
 ['单核输入复核','完成','533新增边界细胞排除；根资格与冻结C4核对'],
 ['Monocle3与动态程序','完成，覆盖限制保留','10模型、基因/程序动态、同分区比较、根与方法敏感性'],
 ['pySCENIC','完成','235 regulons；全部33,098细胞AUCell'],
 ['调控整合','已完成主要计算，正文需据实组织','患者状态差异、轨迹关联、深度调整、靶基因与marker交集'],
 ['稳健性','完成已计划计算；不等同跨队列独立验证','肿瘤限定、根/方法敏感性、留一患者；本次新增同患者2×2网络/标签对照'],
 ['下游交接','依赖登记完成','旧状态不能改名为C4；第三四部分新模型未重跑，属于后续范围'],
 ['图件和客户报告','尚未最终交付','先审内容后绘图的用户指令仍生效']
]
pd.DataFrame(extra,columns=['后续工作项','判定','依据']).to_csv(O/'EXPANDED_PLAN_COVERAGE.tsv',sep='\t',index=False)
matched=read('MATCHED_UNIVERSE_2X2_SUMMARY.tsv');net=read('OLD_NEW_REGULON_TARGET_COMPARISON.tsv')
newtf=read('NEW_ALL_REGULON_CONTRASTS.tsv');traj=read('SELECTED_TF_ALL_TRAJECTORY_SENSITIVITIES.tsv')
major=newtf[(newtf.scope=='all_tissues')&(newtf.endpoint=='C4_minus_other_macrophages')&(newtf.TF.isin(['CEBPD','JUN','JUNB','PPARG','TFEB','THRA','CEBPA']))]
major=major[['TF','n_patients','median_effect','positive','negative','P','FDR_all_regulons']]
keytraj=traj[(traj.scope=='all_tissues')&(traj['mode']=='primary')&(traj.pt_version=='pt_primary')&(traj.adjustment=='rho')&(traj.TF.isin(['JUN','JUNB','CEBPD','PPARG','TFEB','THRA']))][['TF','n_patients','median_effect','positive','negative','P','FDR_all_regulons','FDR_prespecified_3TF']]
# Verify quoted marker numbers directly.
markers=pd.read_csv(R/'04_markers_annotation/run_v4_macrophage_only/full_marker_tables/cluster_4_all_genes.tsv.gz',sep='\t')
markers[markers.gene.isin(['TYMP','SPP1','GPNMB','CCL18'])].to_csv(E/'C4_KEY_MARKERS_DIRECT.tsv',sep='\t',index=False)
tymp=markers[markers.gene=='TYMP'].iloc[0];spp1=markers[markers.gene=='SPP1'].iloc[0]
assert np.isclose(tymp.median_patient_log2CPM_difference,-.110926,atol=1e-5)
report="""# 第二部分原方案完成度与初版结果对照审阅

日期：2026-09-15。审阅范围：第二部分巨噬细分、功能、轨迹和调控；不将第一部分2.1–2.7误作第二部分，不重算第三/第四部分。权威目录：phase149_part2_scope_initial_comparison_20260915_v1。

## 核心计算已完成，但尚未完成原方案的全部交付

第二部分所需的巨噬细分、Monocle3、pySCENIC及本轮扩展的功能/患者统计已完成。上游验收1,500/1,500、后续验收1,364/1,364通过。原方案要求的轨迹流线图和regulon瀑布热图尚未交付，原因是9月14日用户明确要求先审内容、不生成新图；图件和最终客户正文仍应列为待办。

原方案期待的“普通单核逐步成为TYMP+/SPP1+促癌TAM，JUN/JUNB/CEBPD级联上扬”未被完整支持。这是已完成分析得出的结果，不应记成算法没跑完，也不能按预设方向改写结果。

与最早初版相比，新版补齐了从无监督分群到完整marker、患者/研究核查的可追溯链条；与8月规则状态版相比，重新发现了不同的C4群，并采用保留分区的轨迹和新网络。它不是旧SPP1群的直接延续。

## 原方案逐项核对

原始《分析待做（新）》与《分析待做》的第二部分文字一致；已核对原DOCX SHA256与既有提取记录。第二部分明确要求巨噬细分、Monocle3、pySCENIC、轨迹流线图、regulon瀑布热图。四种轨迹算法、CytoTRACE/CellRank、metaVIPER和恶性细胞四通路通量属于第一部分；不能列为第二部分缺项。

{coverage}

后续获授权的扩展工作（功能定位、单核复核、轨迹、动态程序、SCENIC、整合稳健性及下游交接）覆盖表另见EXPANDED_PLAN_COVERAGE.tsv。下游交接完成指已说明新旧标签依赖，不代表第三/第四部分已按C4重做。

## 三套SPP1定义不同，新版不能直接继承初版标签

这里同时比较两套历史结果，避免“初版”指代混淆。7月七亚群版是原巨噬Seurat对象；8月规则状态版是此前Phase33/34轨迹和SCENIC的直接前版。9月14日当日已完成的C4发现版与本次SCENIC版使用相同冻结C4成员，SCENIC没有再筛选或修改C4。

| 项目 | 7月最早七亚群版 | 8月规则状态版（Phase33/34） | 当前重构版 |
| --- | --- | --- | --- |
| 基础对象/标签 | 41,105巨噬、7亚型；最终分辨率0.1 | 50,761髓系、40患者；SPP1/C1QC患者内程序阈值 | 33,098轨迹/SCENIC细胞、39患者；当前基础仍75文库 |
| SPP1群 | SPP1+Macro 11,290细胞 | SPP1-state 8,892细胞 | 无监督C4 2,674细胞、出现于37患者 |
| 主要研究构成 | 旧SPP1群93.92%来自GSE242299 | 丰度受患者内阈值约束 | C4最大来源GSE207493占51.05%，五研究均出现 |
| 与C4重叠 | 36细胞，占C4的1.35% | 1,681细胞，占C4的62.86% | C4身份冻结 |
| 分群依据 | 有聚类对象；原完整marker命名链与4,204细胞排除依据未恢复 | 预先定义程序高分状态，不是无监督SPP1群发现 | 无监督后完整marker、患者/研究组成、稳定性检查 |
| 稳定性结论 | 存在于五研究不能证明跨研究稳健 | 条件状态可分析，不能冒称自然群体比例 | 17候选全局稳定门槛未通过；C4核心有稳定性支持，整体仍探索性 |

对旧七亚群的新增精确对照：C4中1,430细胞原标为CX3CR1+Macro，349为HLA-II+Macro，325不在旧最终巨噬对象；只有36原标SPP1+Macro。完整八类去向见NEW_C4_VS_JULY_SEVEN_SUBTYPES.tsv。36例交集已用独立csv读取再次核对。重叠率不是以新版为真值计算的旧标签准确率。

## C4具有明确的SPP1富集与功能差异，TYMP没有同步富集

相对其他成熟巨噬候选，C4的SPP1患者配对log2(CPM+1)差中位+1.448645，26名有效患者中24名正向，P=6.67781e-5、FDR=0.000544714；GPNMB、CCL18等多基因支持该群身份。C4并非每个细胞均SPP1阳性，SPP1检出率约47.79%；群命名应表达相对富集。

TYMP在C4有表达，但患者配对差中位−0.110926，P=0.0308635、FDR=0.0586967，为名义降低。不能把“TYMP可检出”替换成“TYMP相对富集”，也不能由此声称TYMP诱导C4。

去除SPP1/TYMP后，26名患者的KEGG胆固醇代谢和溶酶体程序均26/26升高，中位差分别+0.624684、+0.349194，P均8.80367e-6、FDR均2.00098e-5；抗原呈递25/26降低。吞噬体程序14/26升高、P=0.542159，当前数据不支持统一升高。缺氧、炎症、干扰素及巨噬内血管生成程序整体降低；Hallmark脂肪酸代谢轻度升高。不同数据库的胆固醇程序含义和基因集不同，不能混为一项。

这些结果补充了初版缺少的当前C4功能证据。它们不能直接反证或验证旧细胞通信中的内皮响应，因为“巨噬内血管生成程序”与“巨噬到内皮通信”不是同一端点。

## C4晚于单核的排序方向保留，但普遍终末分化不成立

初版Phase33利用冻结X_scVI作表示，learn_graph设use_partition=FALSE、close_loop=TRUE；新版使用研究内HVG/PCA，去除SPP1/TYMP拟合特征，use_partition=TRUE、close_loop=FALSE，缺少可信单核起点的分区保留无有效拟时序。输入、坐标、图结构和群定义均改变。

初版五研究SPP1-minus-Monocyte患者中位方向均正，但部分研究仅n=5–7；不能写五研究均显著。新版主模型C4-minus-Monocyte在8名有效患者中8/8正向，中位差+0.555638，属于小样本探索性排序。原表P=0.0142662、FDR=0.0285324完整保留在技术附件。去应激/增殖敏感性模型19/19正向，中位差+0.623473，P=0.000143020、FDR=0.000286040。

关键覆盖事实：主模型只有442/2,674个C4细胞（16.53%）具有有效拟时序；敏感性模型为1,743/2,674（65.18%）。主模型中GSE207493的1,361个C4位于缺少足够可信单核根的分区，敏感性重建后该研究1,365个C4进入有根图。应保留这一结构敏感性，不能把小部分可估计C4的排序外推为全部C4。

C4相对其他巨噬的主模型效应仅4/8正向，中位+0.028087；去应激/增殖模型7/19正向，中位−0.069368、P=0.324166。当前数据不支持C4普遍比其他成熟巨噬更晚。DPT与Monocle3在这一端点不一致，须保留，不能只选择支持晚期的DPT结果。

## 原定三个TF并非C4特异上扬，新发现候选需分别解释

当前C4相对其他巨噬的主要regulon差异如下。n为独立患者；effect为同一regulon的患者内平均AUCell差中位数，不用于跨TF比较生物学驱动力大小。

{major}

JUN/JUNB/CEBPD分别25/26、26/26、26/26降低；PPARG、TFEB、THRA、CEBPA均26/26升高。以上方向在五研究的患者中位数均一致，但各研究有效n=2/3/4/5/12，不能据此称五个独立验证队列均显著；留一患者的方向均保留。

在全组织主Monocle3轨迹中的变化：

{traj}

JUN活性有正向轨迹关联；JUNB主轨迹P=0.0186352，完整235-regulon家族FDR=0.0568737，仅名义显著，但原预设三TF家族FDR=0.0279528通过，两种口径均须保留。JUNB在DPT方向转负、去应激/增殖Monocle3中未获支持，因此不能称跨方法级联上扬。CEBPD主轨迹无支持。群间降低和分区内轨迹相关升高可以同时成立：后者涵盖同分区其他分支，不是C4专属路径，也不代表C4达到最高活性。

PPARG在C4中升高，在主Monocle3、DPT、去应激/增殖模型中均呈正向关联，可作为当前优先阐述的调控候选。其25个可评估预测靶基因中19个与C4上调marker重叠（ORA P=6.25852e-15、FDR=6.84071e-14），为同数据的调控对应证据。TFEB的状态差异与上调marker重叠明确，但主Monocle3轨迹无升高证据，不能表述为持续上升驱动。THRA网络仅5个靶基因，虽状态和部分动态结果支持，暂列小靶集候选；不因AUCell差较大就称其为核心驱动。

## 与8月初版比较，CEBPD和JUNB动态结论发生变化

旧版CEBPD随Monocle3正向（40患者，rho中位0.272055，P=1.04050e-5，三TF FDR=3.12149e-5）；旧SPP1-state相对单核升高（35患者，差+0.003551）。新版CEBPD在C4相对单核降低（24患者，差−0.011555），主轨迹rho中位0.032499、P=0.660282。旧版“CEBPD活化伴随单核向SPP1状态排序”的结论不能直接迁移给C4。

旧JUNB主轨迹下降（rho−0.282797，P=0.00116718），新版主Monocle3方向转正但方法间不一致。JUN从旧版主轨迹无支持转为新版主Monocle3正向。旧SPP1-minus-C1QC和新C4-minus-other-macrophages使用不同比较群，不能把两个效应直接相减或声称重复验证。

新旧网络也不同：210→235 regulons。JUN靶集Jaccard0.291、JUNB0.444、CEBPD0.106；THRA从157靶基因变为5个，TFEB旧网络未检测到。相同TF名称不代表相同计分基因集合。

### 同一细胞与患者集合的2×2对照确认定义和网络均有影响

为避免仅凭不同样本量解释方向变化，本审阅固定33,098个共同可评估细胞，并限定两种目标群与单核均达到20细胞的同一24名患者。比较旧SPP1-state/新C4两种标签与旧/新网络两种评分，均以同一患者单核为对照；不重新训练网络。本次为事后敏感性，每个网络×标签的三TF分别BH，不替代原主要统计。

{matched}

以CEBPD为例：旧网络+旧标签为正（+0.003672，18/24正）；旧网络+新C4为弱正且无支持（+0.000751，14/24正）；新网络+旧标签为负（−0.007130，21/24负）；新网络+新C4为负（−0.011555，22/24负）。这表明可见变化同时随状态定义和网络计分而改变；不能只归因于患者数减少，也不能据此断言某个网络是实验真值。JUNB在这四个条件下相对单核均负向，群间降低比其动态方向更一致。

## 当前审核稿需要怎样修订

1. 将“全部完成”限定为核心计算完成；正式图件和最终正文仍待内容审核后制作。
2. 在主结果保留主模型16.53% C4伪时间覆盖、敏感性65.18%及“未普遍晚于其他巨噬”的事实；不是仅放入无关QC角落。
3. 主线改为“C4身份与功能差异→相对单核的可估计排序→差异调控程序”；不沿用“TYMP规训终末SPP1状态、三TF级联上扬”。
4. PPARG可优先进入候选叙述；TFEB为状态关联，THRA为小靶集候选。保留JUN/JUNB/CEBPD降低及动态方法差异。
5. 原技术审核稿未充分呈现同定义比较与关键敏感性，应补入本次对照；完整统计留在附表，正文围绕结论组织。
6. 第三/第四部分旧SPP1-state的通信、空间后验不能改名为C4。新C4下游需另做明确端点的分析；本次仅登记依赖，不算第二部分核心计算缺项。

## 方法与审计边界

本轮读取原方案提取文本并核对DOCX哈希，审阅初版对象追溯、Phase33/34原脚本和结果，核对新版完整统计表、验证状态及冻结C4身份。新增计算仅为标签重叠、覆盖率/研究方向汇总、靶集重叠及固定患者的现成AUCell评分对照；没有重做聚类、Monocle3、SCENIC或第三/第四部分模型，没有绘图。统计以独立患者为单位；n<10按探索性，n=2研究方向仅作描述，主要检验和审计性三TF家族分开。

当前证据是同一发现数据中的条件性关联。横断面轨迹不等于谱系转化，推断regulon不等于TF结合实验，靶基因/marker重叠不是独立验证。不同网络AUCell绝对值和不同细胞群的效应不作直接定量差值推断。上述边界不否定已确认的C4身份、表达和程序差异。

## 可复核材料

原方案覆盖：ORIGINAL_PLAN_COVERAGE.tsv；后续授权覆盖：EXPANDED_PLAN_COVERAGE.tsv。
证据：evidence/ 下保留初版/新版完整核心统计、原方案第二部分文本、全部轨迹敏感性、七亚群与规则状态两套重叠、TF研究方向/LOPO、靶集变化以及2×2对照。
来源哈希：SOURCE_MANIFEST.tsv；本轮脚本：01_collect_evidence.py、02_matched_comparison.py、03_write_review.py。
本报告为完成度与结果对照审阅，不是已完成排版的客户最终报告。
"""
report=report.replace('{coverage}',table(coverage)).replace('{major}',table(major)).replace('{traj}',table(keytraj)).replace('{matched}',table(matched))
assert '???' not in report and '{coverage}' not in report
(O/'PART2_COMPLETENESS_AND_INITIAL_COMPARISON_CN.md').write_text(report,encoding='utf-8')
# Copy the executed collection script into the audit folder for reproducibility.
shutil.copy2(R/'scripts/36_part2_scope_initial_audit_20260915_v1.py',O/'01_collect_evidence.py')
summary="""# 第二部分审阅结论

第二部分的核心计算已经完成，原方案的全部交付尚未完成。待交付的是轨迹流线图、regulon瀑布热图和最终结果正文，均承接“先审内容、后绘图”的安排。原方案的三TF一致级联上扬不受当前结果支持，不能记为因算法未运行而待补。

新版C4与最早SPP1+Macro只重叠36细胞（C4的1.35%），与8月规则SPP1-state重叠1681细胞（62.86%）；它不是旧群体的直接延续。新版在无监督分群、完整marker、患者和研究核查方面有更完整的证据链，17候选整体仍为探索性分类。

C4的SPP1、胆固醇代谢和溶酶体程序升高；TYMP相对其他巨噬没有同步富集。主轨迹中8名有效患者均显示C4晚于单核，但只有442/2674（16.53%）C4细胞有有效拟时序；去应激/增殖后覆盖65.18%、19/19患者方向一致。当前数据不支持C4普遍晚于其他成熟巨噬。

新SCENIC中JUN/JUNB/CEBPD均在C4相对其他巨噬降低。JUN有主轨迹正向关联，JUNB方向依赖方法，CEBPD主轨迹无支持。PPARG在群间差异及多种排序中方向更一致，可作为调控候选；TFEB支持状态差异，THRA仅5靶基因，需保持小靶集候选定位。

固定同一24名患者的旧/新标签×旧/新网络对照表明，CEBPD变化同时随群定义和网络改变；不能简单解释为样本数变化。旧结果与新结果保留各自定义，不混成同一群体的连续证据。

建议第二部分正文围绕“C4身份与功能特征、可估计的单核—巨噬排序、差异调控程序”组织；图件按真实结果呈现。第三/第四部分新C4通信和空间验证是下一阶段的定义衔接，不作为本部分额外缺跑算法。
"""
(O/'REVIEW_SUMMARY_CN.md').write_text(summary,encoding='utf-8')
print('REPORT_READY',len(report),'characters')
print('TYMP',tymp.median_patient_log2CPM_difference,tymp.patient_wilcoxon_p,tymp.BH_FDR_cluster_gene_family)
print('SPP1',spp1.pct_in)

