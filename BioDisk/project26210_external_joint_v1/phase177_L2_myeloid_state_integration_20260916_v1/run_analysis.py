from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,shutil,datetime
r=Path("/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1")
d=r/"phase177_L2_myeloid_state_integration_20260916_v1";d.mkdir(exist_ok=False)
sources={
"dynamic_summary":"part2_rebuild_v1/07_trajectory/dynamics_v2/DYNAMIC_TREND_SUMMARY.tsv",
"dynamic_patient":"part2_rebuild_v1/07_trajectory/dynamics_v2/PATIENT_DYNAMIC_TRENDS.tsv",
"pseudotime_summary":"part2_rebuild_v1/07_trajectory/dynamics_v2/PSEUDOTIME_SUMMARY.tsv",
"partition_qc":"part2_rebuild_v1/07_trajectory/monocle3_v3/PARTITION_QC.tsv",
"program_genes":"part2_rebuild_v1/07_trajectory/dynamics_v2/PROGRAM_COVERAGE.tsv",
"same_library_bridge":"phase154_same_library_TYMP_C4_test_20260915_v1/TYMP_C4_AND_SPP1_RESULTS.tsv",
"old_bridge":"phase151_TYMP_SPP1_threeTF_role_20260915_v1/TUMOUR_TYMP_BRIDGE_SUMMARY.tsv",
"tf_summary":"phase151_TYMP_SPP1_threeTF_role_20260915_v1/TF_GENE_SUMMARY.tsv",
"composition_sensitivity":"phase159_metabolic_parts_C4_value_20260915_v1/MACROPHAGE_C4_COMPOSITION_SENSITIVITY.tsv",
"c4_contribution":"phase159_metabolic_parts_C4_value_20260915_v1/SAME_LIBRARY_C4_CONTRIBUTION.tsv"}
(d/"PLAN_BEFORE_SUMMARY_CN.md").write_text("""# 本轮范围与固定口径
先复核冻结结果，不重跑轨迹、不重聚类、不改变基因集、不画图。
肿瘤内深度调整为本轮叙事重点；全组织、未调整、替代根及DPT完整保留。
固定查看SPP1、TYMP、Cholesterol metabolism、Lysosome；它们均为前期已指定对象，不依据本轮P值筛选。
新增计算仅为共同患者的描述性对照：主Monocle3分别与去应激/周期、DPT配对，报告同向数、中位效应和逐研究留出方向；不新增P/FDR。
共同患者不等于共同细胞，不将覆盖率改变消除后的比较误称相同细胞敏感性。
从全部患者动态行复算历史摘要的样本数/中位数/方向，不改原P或FDR。历史统计完整留存。
接收端整体SPP1状态和C4身份分开，肿瘤发送端TYMP与巨噬侧TYMP分开。新旧17/16患者桥接不混用。
""",encoding="utf8")
snap=d/"source_tables";snap.mkdir()
manifest=[];tables={}
for k,rel in sources.items():
 p=r/rel;b=p.read_bytes();q=snap/(k+".tsv");q.write_bytes(b)
 manifest.append({"key":k,"source":str(p),"sha256":hashlib.sha256(b).hexdigest(),"bytes":len(b)})
 tables[k]=pd.read_csv(q,sep="\t")
(d/"INPUT_MANIFEST.json").write_text(json.dumps(manifest,indent=2),encoding="utf8")
p=tables["dynamic_patient"];s=tables["dynamic_summary"]
checks=[]
keys=["mode","scope","pt_version","family","feature"]
for i,row in s.iterrows():
 a=p.copy()
 for k in keys:a=a[a[k]==row[k]]
 vals=a[row.adjustment].dropna().to_numpy()
 checks.append({"row":int(i),"n_ok":len(vals)==row.n_patients,"median_ok":bool(np.isclose(np.median(vals),row.median_effect,rtol=1e-10,atol=1e-12)),"positive_ok":int((vals>0).sum())==row.positive,"negative_ok":int((vals<0).sum())==row.negative})
check=pd.DataFrame(checks);check.to_csv(d/"HISTORICAL_SUMMARY_VALIDATION.tsv",sep="\t",index=False)
features=["gene:SPP1","gene:TYMP","program:Cholesterol metabolism","program:Lysosome"]
selected=s[(s.scope=="tumor_only")&s.feature.isin(features)]
selected.to_csv(d/"TUMOR_DYNAMIC_EVIDENCE.tsv",sep="\t",index=False)
records=[];pairs=[];los=[]
for f in features:
 a=p[(p.scope=="tumor_only")&(p.feature==f)&(p["mode"]=="primary")&(p.pt_version=="pt_primary")][["patient","study","depth_adjusted_rho"]].rename(columns={"depth_adjusted_rho":"primary"})
 for label,mode,pt in [("without_stress_cycle","without_stress_cycle","pt_primary"),("DPT","primary","dpt_pseudotime")]:
  b=p[(p.scope=="tumor_only")&(p.feature==f)&(p["mode"]==mode)&(p.pt_version==pt)][["patient","study","depth_adjusted_rho"]].rename(columns={"depth_adjusted_rho":"sensitivity"})
  z=a.merge(b,on=["patient","study"],validate="one_to_one").dropna()
  z.insert(0,"comparison",label);z.insert(0,"feature",f);pairs.append(z)
  records.append({"feature":f,"comparison":label,"n_patients":len(z),"n_studies":z.study.nunique(),"primary_median":z.primary.median(),"sensitivity_median":z.sensitivity.median(),"primary_positive":int((z.primary>0).sum()),"sensitivity_positive":int((z.sensitivity>0).sum()),"both_positive":int(((z.primary>0)&(z.sensitivity>0)).sum()),"same_nonzero_sign":int(((np.sign(z.primary)==np.sign(z.sensitivity))&(z.primary!=0)).sum()),"median_paired_difference":(z.sensitivity-z.primary).median()})
  for study in sorted(z.study.unique()):
   sub=z[z.study!=study]
   los.append({"feature":f,"comparison":label,"omitted_study":study,"remaining_n":len(sub),"primary_median":sub.primary.median(),"sensitivity_median":sub.sensitivity.median()})
pd.DataFrame(records).to_csv(d/"COMMON_PATIENT_SENSITIVITY_SUMMARY.tsv",sep="\t",index=False)
pd.concat(pairs).to_csv(d/"COMMON_PATIENT_SENSITIVITY_VALUES.tsv",sep="\t",index=False)
pd.DataFrame(los).to_csv(d/"COMMON_PATIENT_LEAVE_STUDY.tsv",sep="\t",index=False)
qc=tables["partition_qc"].groupby("mode",as_index=False)[["cells","C4_cells","finite_primary","C4_finite"]].sum()
qc["C4_coverage"]=qc.C4_finite/qc.C4_cells
qc.to_csv(d/"TRAJECTORY_COVERAGE.tsv",sep="\t",index=False)
pg=tables["program_genes"];pg["contains_TYMP"]=pg.gene_symbols.str.split(";").map(lambda z:"TYMP" in z);pg["contains_SPP1"]=pg.gene_symbols.str.split(";").map(lambda z:"SPP1" in z)
pg[["program","genes","contains_TYMP","contains_SPP1"]].to_csv(d/"PROGRAM_TARGET_GENE_AUDIT.tsv",sep="\t",index=False)
tf=tables["tf_summary"];tf[(tf.scope=="tumor_only")&(tf.gene=="SPP1")].to_csv(d/"SPP1_TF_EVIDENCE.tsv",sep="\t",index=False)
c=tables["c4_contribution"]
print("checks",len(check),check.iloc[:,1:].all().to_dict())
print("COMMON\n"+pd.DataFrame(records).to_csv(sep="\t",index=False))
print("COVERAGE\n"+qc.to_csv(sep="\t",index=False))
print("C4 medians",c.C4_fraction.median(),c.C4_SPP1_RNA_share.median())
print("TF\n"+tf[(tf.scope=="tumor_only")&(tf.gene=="SPP1")&(tf.adjustment=="depth_adjusted")].to_csv(sep="\t",index=False))
print("LOS\n"+pd.DataFrame(los)[pd.DataFrame(los).feature=="gene:SPP1"].to_csv(sep="\t",index=False))
assert check.iloc[:,1:].all().all()
