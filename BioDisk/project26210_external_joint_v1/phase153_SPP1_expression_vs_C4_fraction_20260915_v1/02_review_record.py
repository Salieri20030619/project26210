from pathlib import Path
import json,datetime,shutil
import pandas as pd,numpy as np
from scipy.stats import spearmanr,false_discovery_control
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');R=P/'part2_rebuild_v1';O=P/'phase153_SPP1_expression_vs_C4_fraction_20260915_v1'
d=pd.read_csv(O/'PATIENT_SPP1_C4_DECOMPOSITION.tsv',sep='\t')
old=pd.read_csv(P/'phase151_TYMP_SPP1_threeTF_role_20260915_v1/TUMOUR_TYMP_MACROPHAGE_MATCHED_PATIENTS.tsv',sep='\t').set_index('patient')
assert d.patient.is_unique and len(d)==17
assert np.allclose(d.SPP1_all,old.loc[d.patient,'SPP1_RNA_mean'])
assert np.allclose(d.C4_fraction,old.loc[d.patient,'C4_fraction'])
assert np.max(np.abs(d.SPP1_all-d.C4_weighted_component-d.nonC4_weighted_component))<1e-12
rows=[]
for col in ['SPP1_all','C4_fraction','SPP1_C4_mean','SPP1_nonC4_mean','C4_weighted_component','nonC4_weighted_component']:
 z=d.dropna(subset=[col])
 if col=='SPP1_C4_mean':z=z[z.n_C4>=20]
 if col=='SPP1_nonC4_mean':z=z[z.n_nonC4>=50]
 r,p=spearmanr(z.TYMP,z[col]);rows.append(dict(readout=col,n_patients=len(z),rho=r,P=p))
t=pd.DataFrame(rows);t['FDR_posthoc_six_readout_diagnostic']=false_discovery_control(t.P.to_numpy());t.to_csv(O/'DECOMPOSITION_DIAGNOSTIC_CORRELATIONS.tsv',sep='\t',index=False)
report="""# 为什么TYMP—SPP1名义显著，而TYMP—C4比例不显著

## 两个读出不同，且观测相关强度不同

前一项为17患者肿瘤TYMP RNA与全部成熟巨噬细胞SPP1平均表达的相关：rho=0.507353，P=0.037633，原五读出FDR=0.188164，仅名义显著。
后一项为相同17患者肿瘤TYMP RNA与新C4占成熟巨噬比例的相关：rho=0.308824，P=0.227786，原五读出FDR=0.569464。
患者数相同，第二项相关较弱，因而P值更大。不能把“一项P<0.05而另一项不满足”当作两项相关强度已被正式证明不同。本轮没有进行二者差异的正式检验。

## SPP1不是C4独占的标志，C4比例不能代替总体SPP1表达

在该17患者中，C4占成熟巨噬的比例中位数7.55%。SPP1阳性巨噬中属于C4的比例中位数14.29%。C4内SPP1表达均值的患者中位数0.730675，非C4为0.386303；C4相对富集SPP1，但非C4人数更多、也表达SPP1。不能把所有SPP1阳性细胞归为C4。

就本次采用的每细胞log1pCP10K表达指数而言，C4对总体SPP1表达指标的贡献比例中位数12.93%，说明多数该指标的贡献来自非C4。这是表达指数的组成分解，不是测得的RNA分子总量比例。

每名患者满足精确恒等式：
总体SPP1平均表达 = C4比例 × C4内SPP1平均表达 + 非C4比例 × 非C4内SPP1平均表达。

因此，总体SPP1可以通过C4比例、C4内表达、或其他巨噬内表达的变化而改变。仅检验C4比例只覆盖其中一项。

## 当前分解提示总体信号包含非C4贡献，没有证明C4与TYMP无关

TYMP与非C4巨噬SPP1平均表达的相关rho=0.441，n=17；与C4加权表达贡献rho=0.304，n=17；与满足至少20个C4细胞的群内SPP1平均表达rho=0.220，n=13。完整诊断性P值和事后六读出BH见DECOMPOSITION_DIAGNOSTIC_CORRELATIONS.tsv，不能将其重新解释为原Phase151的统计家族，也不据这些拆分选择新的显著结论。

这些结果解释了为什么“总体SPP1表达关联”不能直接等同“C4丰度关联”。它们不否定新C4这个候选状态，也不提供TYMP驱动C4的因果证据；目前的比例结果只能说当前数据未支持该17患者端点。

## 下一步仍按新C4优先计划执行

继续固定C4身份，在16名同肿瘤文库配对患者中直接检验恶性TYMP与C4比例；保留C4=0样本和全部成熟巨噬分母。SPP1总体表达与检出作为补充读出，区分“某一状态的细胞增多”与“细胞表达增强”。不因已有P值改变分群，也不以旧SPP1状态阴性替代新C4检验。

本次为已有17患者结果的组成解释，不是16患者计划的新结果。未重新训练模型、未绘图。患者数据与原Phase151逐一一致；表达分解最大误差小于1e-12。
"""
(O/'SPP1_EXPRESSION_VS_C4_FRACTION_CN.md').write_text(report,encoding='utf-8')
B=O/'state_before_update';B.mkdir(exist_ok=False)
for q in [P/'CODEX_HANDOFF.md',P/'PROJECT_STATUS.md',P/'DECISIONS.md',P/'RESULTS_MANIFEST.tsv']:shutil.copy2(q,B/q.name)
note="""## 2026-09-15 Phase153：SPP1表达与C4比例两个端点的差异

用户询问为何TYMP—SPP1名义显著而TYMP—C4比例不显著。已对Phase151相同17患者精确分解，表与原结果一致。C4比例中位7.55%；SPP1阳性巨噬中C4比例中位14.29%；C4对本次log表达均值指标贡献比例中位12.93%。总体SPP1为C4与非C4的加权均值，SPP1不是C4独占标志。
rho0.507/P0.0376和rho0.309/P0.228对应不同读出；原FDR分别0.188/0.569均未过。不能以一个P过线推断两个相关强度的差异显著。诊断性六读出另算事后BH，不改变旧家族。报告phase153_SPP1_expression_vs_C4_fraction_20260915_v1/SPP1_EXPRESSION_VS_C4_FRACTION_CN.md。
下一任务仍为Phase152v2的16患者同文库TYMP—新C4优先计划；当前只完成17患者解释，无新模型/新图。
"""
for fn in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
 q=P/fn;q.write_text(note+'\n'+q.read_text(),encoding='utf-8')
with (P/'RESULTS_MANIFEST.tsv').open('a') as f:f.write('\t'.join(['phase153_metric_decomposition','phase153','diagnostic',str(O/'SPP1_EXPRESSION_VS_C4_FRACTION_CN.md'),datetime.datetime.now(datetime.timezone.utc).isoformat(),'PASS','17 patient identity and weighted mean equality verified','No model refit; no figures'])+'\n')
(O/'STATUS.json').write_text(json.dumps(dict(status='COMPLETE_METRIC_DECOMPOSITION_REVIEW',patients=17,patient_identity_matches=True,weighted_mean_identity_pass=True,figures_created=False,model_rerun=False),indent=2))
print('COMPLETE',d.C4_fraction.median(),d.C4_share_SPP1_detected.median())

