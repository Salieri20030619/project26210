from pathlib import Path
import json,hashlib,datetime,shutil
import pandas as pd,numpy as np
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/part2_rebuild_v1');O=R/'11_reports/remaining_analysis_review_v1';O.mkdir(exist_ok=False)
F=R/'06_functional_programs/run_v1';T=R/'07_trajectory/monocle3_v3';D=R/'07_trajectory/dynamics_v2';N=R/'08_regulons/network_v1';S=R/'08_regulons/statistics_v1';X=R/'09_downstream_crosswalk/remaining_analysis_v1';V=R/'provenance/remaining_analysis_validation_v1';assert json.loads((V/'STATUS.json').read_text())['status']=='PASS'
def table(df):
 def val(v):
  if isinstance(v,float):return 'NA' if np.isnan(v) else f'{v:.6g}'
  return str(v).replace('|','/').replace('\n',' ')
 return '| '+' | '.join(map(str,df.columns))+' |\n| '+' | '.join(['---']*len(df.columns))+' |\n'+'\n'.join('| '+' | '.join(val(v) for v in row)+' |' for row in df.itertuples(index=False,name=None))
f=pd.read_csv(F/'C4_CONTENT_REVIEW_KEY_RESULTS.tsv',sep='\t');pq=pd.read_csv(T/'PARTITION_QC.tsv',sep='\t');pt=pd.read_csv(D/'PSEUDOTIME_SUMMARY.tsv',sep='\t');dt=pd.read_csv(D/'DYNAMIC_TREND_SUMMARY.tsv',sep='\t');rc=pd.read_csv(S/'AUC_CONTRAST_SUMMARY.tsv',sep='\t');rt=pd.read_csv(S/'AUC_TRAJECTORY_SUMMARY.tsv',sep='\t');ov=pd.read_csv(X/'C4_LEGACY_OVERLAP_METRICS.tsv',sep='\t');validation=json.loads((V/'STATUS.json').read_text());network=json.loads((N/'STATUS.json').read_text())
primary_pt=pt[(pt['mode']=='primary')&(pt.scope=='all_tissues')&(pt.pt_version=='pt_primary')]
primary_dt=dt[(dt['mode']=='primary')&(dt.scope=='all_tissues')&(dt.pt_version=='pt_primary')&(dt.adjustment=='rho')]
tf3=rc[rc.TF.isin(['JUN','JUNB','CEBPD'])];tf3t=rt[rt.TF.isin(['JUN','JUNB','CEBPD'])&(rt['mode']=='primary')&(rt.pt_version=='pt_primary')&(rt.adjustment=='rho')]
rank=rc[(rc.scope=='all_tissues')&(rc.endpoint=='C4_minus_other_macrophages')&(rc.n_patients>=10)&(rc.P<.05)].copy();rank['abs_effect']=rank.median_effect.abs();selected=pd.concat([rank[rank.median_effect>0].sort_values('abs_effect',ascending=False).head(10),rank[rank.median_effect<0].sort_values('abs_effect',ascending=False).head(10)]);selected.to_csv(O/'DISCOVERY_TF_CONTENT_SHORTLIST.tsv',sep='\t',index=False)
cols=['n_patients','median_effect','positive','negative','P'];ar=cols+['FDR_all_regulons','reporting_tier']
text='''# 第二部分剩余分析内容审核稿（无图）

本稿供内容审核，尚未获准绘图。原第二部分和新分群发现结果保留，当前工作独立写入part2_rebuild_v1。用户确认名义显著可作为探索性支持；本稿同时保留效应方向、患者一致性、独立n、P与FDR。所有主要阴性及反向结果保留。

## 1. C4突出胆固醇处理及溶酶体程序，整体功能并非全面增强

冻结C4共2,674细胞，出现于37患者，其中2,437肿瘤来源、237正常组织来源；在26名有效患者内相对其他巨噬候选的SPP1差中位1.448645，24/26正向，P=6.67781e-5、FDR=0.000544714。TYMP普遍检出但配对差轻度降低（-0.110926，P=0.0308635，FDR=0.0586967），不能命名为TYMP富集候选。

去除SPP1/TYMP后的患者级通路差异见下表。胆固醇代谢和溶酶体五研究中位方向均升高；抗原呈递五研究均降低。部分吞噬相关上调基因有ORA富集，但患者整体程序不一致升高，两种端点分别解释。不同数据库的胆固醇程序并不等同。通路评分反映相对转录表达，不是实验功能或代谢通量。

{functional}

完整说明：06_functional_programs/run_v1/CONTENT_REVIEW_CN.md；原始程序差、ORA上下调、研究分层及覆盖表均保留。

## 2. 单核与巨噬轨迹输入已净化并冻结

轨迹输入33,098细胞、39患者、70有效文库，来自已审计75文库基础。对新增单核/C1边界输入复核复杂度、髓系骨架和非髓系/DC混合信号，排除533细胞；现有12巨噬候选及C4成员保持冻结。可信原始单核细胞提供起点资格，不用SPP1/TYMP表达挑起点。保留C1作为过渡候选，但不自动当作根。

主模型与去除C7/C11/C15应激/增殖群的重建模型分别对五研究拟合。HVG2500/PCA30为研究内重新计算，拟合特征去除SPP1/TYMP；full-library CP10K和log1p后缩放PCA。Monocle3保留分区、不闭环；每分区至少20个合格单核才设置根，没有根的分区保留无有效拟时序。起点为该分区单核PCA中位中心最近细胞；距离排序25%/50%位置的单核作为起点敏感性。

{partitions}

## 3. C4相对单核及其他巨噬的拟时序排序

效应为同患者、同连通分区内C4与比较群的拟时序中位数差，双方各至少20细胞。拟时序在根所在分区内缩放至0–1；不跨断开的分区直接比较。若患者有多个可比较分区，先取分区效应中位数再进入患者统计。全组织主模型：

{pt}

完整起点、去应激/增殖重建、DPT及肿瘤来源限定结果保存在07_trajectory/dynamics_v2/PSEUDOTIME_SUMMARY.tsv。DPT属于替代排序算法的敏感性，不是独立生物学验证。肿瘤限定统计沿用全组织图，不能称肿瘤独立重建。

## 4. 基因和功能程序沿相关分区的变化

以下为含C4的连通分区中，每患者至少50细胞/10个不同拟时序值的患者内Spearman相关中位数；同患者多个分区先合并。模型涵盖同一连通分区的其他分支，不能把全分区相关自动解释为C4专属分支或终末命运。SPP1/TYMP未进入本次轨迹拟合特征；程序评分同时排除这两基因。原始相关及控制原始总counts秩的偏相关分别保存。

{dynamics}

根变化、DPT、去应激/增殖及肿瘤组织限定完整结果见DYNAMIC_TREND_SUMMARY.tsv；逐研究方向见DYNAMIC_STUDY_DIRECTIONS.tsv。基因动态与功能程序在各自家族内校正，不与轨迹组间差或regulon检验混池。

## 5. 新pySCENIC网络及候选TF状态差异

网络从新净化细胞建立，按患者×当前候选分层每层至多100细胞，未使用旧SPP1/C1QC规则标签抽样。GRNBoost2全人TF列表，hg38 v10两套motif排名数据库剪枝，AUCell在全部33,098细胞计算。网络包含{regulon_n}个regulon、{edge_n}条目标关系。软件、来源数据库和SHA256见08_regulons/network_v1。

原方案三个TF的所有组间比较如下。TF RNA表达与regulon活性是不同测量层面；不能用RNA下降直接替代活性结论，也不要求三者全部上升。

{tf3}

完整网络中按C4相对其他巨噬的患者效应绝对值列出最多10个升高和10个降低候选，要求n>=10且名义P<0.05。此表为发现性内容候选，不表示已确认直接驱动或指定绘图顺序：

{shortlist}

## 6. 调控与轨迹的对应

JUN/JUNB/CEBPD在含C4分区中的患者内轨迹相关：

{tf3trajectory}

全regulon动态、测序深度控制、起点/方法/组织敏感性及留一患者状态效应见08_regulons/statistics_v1。REGULON_TARGET_MARKER_OVERLAP.tsv保留预测靶基因与C4上下调标记的交集及单独ORA家族；REGULON_TARGET_EDGES.tsv保留全部预测边。推断网络及靶基因重叠提供调控关联，不等于TF结合或实验因果。原定三TF的FDR单列，完整网络发现以全regulon家族q为对应校正口径；不存在的regulon记录不可评估。

## 7. 新旧状态对应与第三/四部分依赖

{crosswalk}

旧规则SPP1-state与无监督C4不能改名等同。旧广义细胞通信及独立bulk临床端点可保留；若将结论具体指向新C4，需要重新估计新标签下的患者级细分通信。旧cell2location的SPP1后验不能改名为C4丰度；新C4丰度需新参考覆盖审查和映射，旧表达程序可保留原端点名称。此次仅登记依赖，未运行第三/四部分新模型，未处理客户湿实验。

## 复核、失败记录与报告边界

新增技术复核{validation_n}项、通过{validation_pass}项，含冻结C4 ID、患者/文库、起点资格、分区有效性、拟时序范围、HVG端点排除、统计家族、原始AUCell独立再聚合及无图检查。功能模块另有72个家族检查及14项独立原始pseudobulk/ORA核验。检查数不是生物学验证样本数。

保留的技术失败：首次功能审阅JSON序列化numpy布尔失败，明确转Python bool修正；最初Monocle3输入将数字PCA列名误读为数据，细胞ID断言在拟合前拦截，显式header=TRUE后新建monocle3_v3；其下游空结果目录保留失败状态。未修改原矩阵或放宽科学门槛。

整体17簇分类仍探索性，C4核心稳定与全局分类稳定分开解释。同数据发现后的P/FDR为条件性证据。横断面轨迹不证明谱系转化或临床时间，AUCell不测量真实TF结合。名义显著按用户口径作为探索性支持；FDR不是继续分析或正文可用性的硬门槛。关键反向与不支持结果不可被下游阳性覆盖。

本稿不含图件，也未自动选择主图面板。内容经用户审核后再进入绘图。
'''
text=text.replace('{functional}',table(f[['library','term','n_patients','median_paired_program_difference','patients_positive','P','FDR']])).replace('{partitions}',table(pq)).replace('{pt}',table(primary_pt[['endpoint']+cols+['FDR','reporting_tier']])).replace('{dynamics}',table(primary_dt[['family','feature']+cols+['FDR','reporting_tier']])).replace('{tf3}',table(tf3[['scope','endpoint','regulon']+ar+['FDR_prespecified_3TF']])).replace('{shortlist}',table(selected[['regulon']+ar])).replace('{tf3trajectory}',table(tf3t[['scope','regulon']+ar+['FDR_prespecified_3TF']])).replace('{crosswalk}',table(ov)).replace('{regulon_n}',str(network['regulons'])).replace('{edge_n}',str(network['edges'])).replace('{validation_n}',str(validation['checks'])).replace('{validation_pass}',str(validation['passed']))
(O/'PART2_REMAINING_ANALYSIS_CONTENT_REVIEW_CN.md').write_text(text)
files=[F/'CONTENT_REVIEW_CN.md',F/'C4_CONTENT_REVIEW_KEY_RESULTS.tsv',D/'PSEUDOTIME_SUMMARY.tsv',D/'DYNAMIC_TREND_SUMMARY.tsv',S/'PRESPECIFIED_TF_CONTRASTS.tsv',S/'PRESPECIFIED_TF_TRAJECTORY.tsv',X/'DOWNSTREAM_DEPENDENCY_MATRIX.tsv',V/'STATUS.json']
manifest=[]
for p in files:manifest.append({'path':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
(O/'CONTENT_SOURCE_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2));(O/'STATUS.json').write_text(json.dumps({'status':'DRAFT_CONTENT_READY_FOR_MODEL_AND_USER_REVIEW','figures_created':False},indent=2));print('CONTENT_DRAFT_CREATED',str(O),flush=True)
