from pathlib import Path
import json
import numpy as np,pandas as pd
from scipy import sparse
from scipy.stats import wilcoxon,hypergeom
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/part2_rebuild_v1');O=R/'06_functional_programs/run_v1';M=R/'04_markers_annotation/run_v4_macrophage_only'
genes=pd.read_csv(M/'PSEUDOBULK_GENE_INDEX.tsv',sep='\t').gene.tolist();ix={g:i for i,g in enumerate(genes)};ri=pd.read_csv(M/'PSEUDOBULK_ROW_INDEX.tsv',sep='\t',dtype={'cluster':str});pb=sparse.load_npz(M/'PATIENT_CLUSTER_PSEUDOBULK_COUNTS.npz');d=np.load(M/'paired_patient_effects/cluster_4.npz');stats=pd.read_csv(O/'C4_ALL_PROGRAMS.tsv',sep='\t');checks=[]
lib={}
for line in (O/'resources/KEGG_2021_Human.gmt').read_text().splitlines():
 b=line.split('\t');lib[b[0]]=set(g.split(',')[0] for g in b[2:] if g)
obs=[];paireffects=[]
for patient in d['patient']:
 mask=ri.patient.eq(patient).to_numpy();inside=np.asarray(pb[mask&ri.cluster.eq('4').to_numpy()].sum(axis=0)).ravel().astype(float);outside=np.asarray(pb[mask&~ri.cluster.eq('4').to_numpy()].sum(axis=0)).ravel().astype(float)
 obs.append((inside>=10)|(outside>=10));paireffects.append(np.log2(inside/inside.sum()*1e6+1)-np.log2(outside/outside.sum()*1e6+1))
elig=set(np.asarray(genes)[np.asarray(obs).sum(axis=0)>=6]);raw=np.asarray(paireffects)
for term in ['Cholesterol metabolism','Lysosome','Antigen processing and presentation','Phagosome']:
 gs=(lib[term]&elig)-{'SPP1','TYMP'};v=raw[:,[ix[g] for g in sorted(gs)]].mean(axis=1);row=stats[(stats.library=='KEGG_2021_Human')&(stats.version=='without_SPP1_TYMP')&(stats.term==term)].iloc[0]
 checks.append({'check':term+'_raw_reaggregation_effect','pass':bool(np.isclose(np.median(v),row.median_paired_program_difference,atol=2e-6))});checks.append({'check':term+'_raw_reaggregation_direction','pass':int((v>0).sum())==row.patients_positive});checks.append({'check':term+'_raw_reaggregation_P','pass':bool(np.isclose(wilcoxon(v,method='approx',correction=True).pvalue,row.P,atol=1e-10))})
ora=pd.read_csv(O/'ALL_ORA_RESULTS.tsv.gz',sep='\t');rr=ora[(ora.cluster==4)&(ora.library=='KEGG_2021_Human')&(ora.version=='without_SPP1_TYMP')&(ora.term=='Lysosome')]
for _,row in rr.iterrows():checks.append({'check':'ORA_'+row.direction,'pass':bool(np.isclose(hypergeom.sf(row.overlap-1,row.universe,row.pathway_genes,row.selected_genes),row.P,atol=1e-10))})
assert all(x['pass'] for x in checks),checks
(O/'INDEPENDENT_RAW_VALIDATION.json').write_text(json.dumps({'status':'PASS','checks':checks},indent=2))
sel=stats[(stats.version=='without_SPP1_TYMP') & stats.term.isin(['Cholesterol metabolism','Lysosome','Phagosome','Antigen processing and presentation','Fatty Acid Metabolism','Hypoxia','Angiogenesis','Inflammatory Response','Interferon Gamma Response'])].copy();sel.to_csv(O/'C4_CONTENT_REVIEW_KEY_RESULTS.tsv',sep='\t',index=False)
text='''# 第二部分剩余分析第1步：C4功能定位（待内容审核，未绘图）

C4的突出特征是胆固醇处理及溶酶体相关表达程序升高，同时伴随炎症/干扰素和抗原呈递程序相对降低。C4不是所有脂质、吞噬、缺氧和血管相关程序同时增强的状态。

与其他巨噬候选在26名患者内比较，去除SPP1/TYMP后的KEGG胆固醇代谢程序差中位+0.624684、26/26正向；溶酶体+0.349194、26/26正向；二者P均8.80367e-6、各自库内FDR=2.00098e-5。抗原处理与呈递差-0.269322，25/26负向。吞噬体差+0.010696，14/26正向，P=0.542159，当前数据不支持其统一升高。

Hallmark脂肪酸代谢差+0.075275，22/26正向，P=0.002508、FDR=0.002850；Hallmark缺氧、炎症和干扰素响应程序降低；血管生成程序差-0.258267、22/26负向，P=0.00011316、FDR=0.00015292。胆固醇相关不同库的程序并不相同：Hallmark Cholesterol Homeostasis轻度下降，不能与KEGG Cholesterol metabolism混为同一端点。完整方向与基因覆盖需一并保留。

方法：下载并冻结Enrichr提供的Hallmark2020、KEGG2021Human、Reactome2022基因集，仅进行本地计算，不上传项目细胞数据。以患者内簇对其他成熟巨噬候选的log2(CPM+1)基因差取通路内均值，汇总患者差中位数；基因覆盖10–500。表达覆盖按原始pseudobulk任一侧>=10 counts且达到至少ceil(0.2*n)或3名患者计算；n=2描述结果以实际两名为上限，不作检验。不是按基因差异显著性筛选通路打分基因。完整程序及同时去除SPP1/TYMP两版本分别保存。

Wilcoxon双侧近似与连续性校正；每簇/数据库/版本的通路家族单独BH。名义P<0.05可作探索性支持，FDR完整保留。另对名义P<0.05且绝对患者中位基因差>=0.25的上下调基因做ORA，背景为表达可评估基因；上下调方向在每簇/数据库/版本共同BH。ORA是基因集合富集，不使用细胞或基因数充当患者样本量，也不与患者程序P值混淆。

已对全部12巨噬候选计算35978行患者程序比较与71956行ORA，并保存逐研究C4方向、逐患者程序差、基因覆盖和完整表。72个BH家族技术检查通过，关键四通路从原始pseudobulk独立再聚合的效应/方向/P及ORA检查通过。

解释边界：这些是同一发现数据中、相对于其他巨噬候选的转录程序差异，不是功能实验、代谢通量、绝对活性或血管效应。macrophage自身Angiogenesis评分下降不能直接代替第三部分对内皮响应的检验；同样不能用旧内皮阳性覆盖本项下降。后续轨迹和regulon需保留此实际方向。

来源：https://maayanlab.cloud/Enrichr/；RESOURCE_REGISTRY.json冻结具体下载端点/时间/SHA256。主要审阅表C4_CONTENT_REVIEW_KEY_RESULTS.tsv；完整结果C4_ALL_PROGRAMS.tsv及ALL_PATIENT_PROGRAM_COMPARISONS.tsv.gz；未生成图件。
'''
(O/'CONTENT_REVIEW_CN.md').write_text(text)
print('FUNCTION_REVIEW_READY',len(checks),'independent checks',flush=True)
