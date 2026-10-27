from pathlib import Path
import json,math,datetime
import numpy as np,pandas as pd
P=Path(__file__).resolve().parents[1]
def bh(values):
 a=np.asarray(values,float);ix=np.argsort(a);v=np.minimum.accumulate((a[ix]*len(a)/np.arange(1,len(a)+1))[::-1])[::-1];out=np.empty(len(a));out[ix]=np.minimum(v,1);return out
def main():
 s=pd.read_csv(P/'tables/DIRECTION_STABILITY.tsv',sep='\t')
 cfg=json.loads((P/'CONFIG.json').read_text())
 rng=np.random.default_rng(262242);rows=[];eligibility=[]
 for mode in ['candidate_specific','common_coverage']:
  for molecule in ['Spermidine','Spermine']:
   for target in cfg['primary_RNA_targets']:
    g=s[(s.coverage_mode==mode)&(s.polarity=='pos')&(s.adduct=='M+H')&(s.molecule==molecule)&(s.RNA_target==target)].copy()
    assert g.patient.is_unique
    for patient in cfg['patient_ids']:
     sub=g[g.patient==patient]
     eligibility.append(dict(coverage_mode=mode,molecule=molecule,RNA_target=target,patient=patient,eligible=bool(len(sub) and sub.primary_rho.notna().all()),reason='mass_candidate_absent' if not len(sub) else 'preferred_mapping_or_peak_not_estimable' if sub.primary_rho.isna().any() else ''))
    v=g.primary_rho.dropna().to_numpy();n=len(v);positive=int((v>1e-12).sum());negative=int((v< -1e-12).sum());nz=positive+negative
    p=None;lo=None;hi=None
    if n>=10 and nz>=10:
     k=min(positive,negative);p=min(1.,2*sum(math.comb(nz,j) for j in range(k+1))/2**nz)
    if n>=10:
     boot=np.median(rng.choice(v,size=(10000,n),replace=True),axis=1);lo,hi=[float(x) for x in np.quantile(boot,[.025,.975])]
    rows.append(dict(family='primary_six' if mode=='candidate_specific' else 'common_coverage_sensitivity_six',coverage_mode=mode,molecule=molecule,RNA_target=target,eligible_patients=n,nonzero_patients=nz,positive_patients=positive,negative_patients=negative,median_patient_rho=float(np.median(v)) if n else None,min_patient_rho=float(v.min()) if n else None,max_patient_rho=float(v.max()) if n else None,median_bootstrap95_low=lo,median_bootstrap95_high=hi,P_exact_sign=p,all_candidates_positive_patients=int(((g.positive_candidates==g.candidates)&(g.estimable_candidates==g.candidates)).sum()),all_candidates_negative_patients=int(((g.negative_candidates==g.candidates)&(g.estimable_candidates==g.candidates)).sum()),orientation_changes_patients=int((g.orientation_stability=='direction_changes_across_candidates').sum()),single_candidate_patients=int((g.candidates==1).sum()),evidence='patient_level_exploratory_association' if p is not None else 'descriptive_insufficient_effective_n',anatomically_validated=False))
 for family in ['primary_six','common_coverage_sensitivity_six']:
  sub=[r for r in rows if r['family']==family];assert len(sub)==6
  q=bh([r['P_exact_sign'] if r['P_exact_sign'] is not None else 1. for r in sub])
  for r,v in zip(sub,q):r['FDR_BH']=float(v) if r['P_exact_sign'] is not None else None
 pd.DataFrame(rows).to_csv(P/'tables/PATIENT_LEVEL_PRIMARY_AND_SENSITIVITY.tsv',sep='\t',index=False)
 pd.DataFrame(eligibility).to_csv(P/'tables/PATIENT_ENDPOINT_ELIGIBILITY.tsv',sep='\t',index=False)
 # Reproduction check on original two cases using unchanged definitions and preferred morphology.
 prior=pd.read_csv(P.parent/'phase238_SpaMTP_regional_association_20260923_v1/tables/DIRECTION_STABILITY.tsv',sep='\t')
 oldkeys=['sample','polarity','molecule','observed_mz','RNA_target','coverage_mode']
 paired=prior.merge(s[s['sample'].isin(['R114_T','Y12_T'])],on=oldkeys,suffixes=('_old','_new'),validate='one_to_one')
 assert len(paired)==len(prior)
 paired['primary_rho_abs_difference']=(paired.primary_rho_old-paired.primary_rho_new).abs()
 paired.to_csv(P/'audit/PILOT_REPRODUCTION_COMPARISON.tsv',sep='\t',index=False)
 # Numerical rerun variation is reported, never used to select candidates or overwrite prior results.
 qa={'status':'PASS','patient_tests_planned':12,'separate_families':[6,6],'effective_n_gate':10,'max_patient_n':max(r['eligible_patients'] for r in rows),'pilot_comparison_rows':len(paired),'pilot_max_abs_rho_difference':float(paired.primary_rho_abs_difference.max()),'new_figures':0,'anatomically_validated':False}
 (P/'PATIENT_SUMMARY_QA.json').write_text(json.dumps(qa,indent=2))
 lines=['# 12例肿瘤空间代谢—转录对应分析','', '本轮统一纳入12个已登记肿瘤患者、24组正负离子采集；细胞/spot/谱图/区域及候选朝向均不扩增患者数。','', '## 患者层主要结果','']
 for r in rows:
  if r['family']!='primary_six':continue
  lines.append(f'- {r["molecule"]}—{r["RNA_target"]}：有效患者{r["eligible_patients"]}/12，患者rho中位数{r["median_patient_rho"]}，正/负方向{r["positive_patients"]}/{r["negative_patients"]}，P={r["P_exact_sign"]}，FDR={r["FDR_BH"]}；跨候选朝向变号患者{r["orientation_changes_patients"]}。')
 lines+=['','## 方法',
 '沿用Phase235谱图解析、固定候选式/加合物与5ppm质量匹配；按TIC归一化。沿用Phase237组织形态配准：8朝向种子、6种掩模设置、有限仿射对照；主相似变换Dice距离最佳值0.03以内的朝向全部保留。仅用组织形态评分，保留显式零基/一基、hires/fullres变换及逆变换验算。',
 '沿用Phase238固定ST四乘四区域；区域至少10映射spots、至少8区域、至少50映射spots及50%QC覆盖；分别报告各候选覆盖与共同覆盖。病例内Spearman仅报告效应，不对spot或区域做推断检验。预设主端点为正离子M+H亚精胺/精胺候选分别对应TYMP、SPP1绝对q05及巨噬内比例，共六项。',
 '患者层效应为各患者rho中位数。至少10例可估计且至少10个非零效应时作双侧精确符号检验；每个六项同类家族独立BH，主要与共同覆盖口径分开。预定但不可估计项保留，BH计算时按P=1保守保留家族大小，不将其报告为实测P。至少10患者时提供10000次患者重采样中位数区间；其余描述性呈现。','',
 '## 证据整合与边界',
 '全部患者/端点/朝向以及不可估计原因留在统计表，未按相关方向筛选样本或配准。原两例结果和本轮重复计算的差异列于audit/PILOT_REPRODUCTION_COMPARISON.tsv；既有Phase238结论不覆盖。',
 '峰为质量相容候选，尚非结构确证。配准是几何候选，仍无独立解剖标志验证。候选朝向的最小最大效应是配准变化范围，不是置信区间；只有一个候选时不能声称跨朝向验证。',
 '空间SPP1沿用旧SPP1-state定义，不改称新C4。polyamine_mta_context为转录背景，不等同实测MTA。不能将组织共定位写成TYMP—多胺—SPP1统一因果链；原血管/免疫主要阴性、既有变号结果及本轮反向结果一并保留。',
 '原RDS与历史结果只读；本轮仅形成文字统计和变换表，不生成新图。']
 (P/'REPORT_CN.md').write_text('\n'.join(lines))
 print(json.dumps(qa),flush=True)
if __name__=='__main__':main()
