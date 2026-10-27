from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
d=Path("/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase183_direct_myeloid_perturbation_search_20260916_v1")
(d/"NORMALIZATION_PLAN_BEFORE_TARGET_EXTRACTION_CN.md").write_text("""# 矩阵资格确认后的计算冻结
GSE261243两个矩阵为原始整数计数，LPS11059基因、TNF10811基因，各12文库=3供体×4条件。LPS供体D562/D563/D565，TNF供体D1/D2/D3；两批身份对应未知，不能合并称6独立供体。
样本为原代人CD14单核细胞，经MCSF20ng/ml培养过夜；IL10 100ng/ml过夜，LPS10ng/ml3h或TNF20ng/ml6h。不是成熟肿瘤巨噬，不含MTA。
仅提供filtered raw counts，未直接提供列名可识别的log2表达。每批全提供矩阵作median-of-ratios size factor标准化(所有样本>0的基因几何均值，逐样本比值中位数)，记录参与基因数。非零计数上四分位数标准化为预设敏感性。两者均为相对归一化，不假装获知全库深度。
主要SPP1；辅助LGALS3/SOCS3/STAT3。基因缺失记不可评估，不填0。主要对比IL10 vs Rest、IL10+刺激 vs 刺激；刺激vsRest、联合vsRest及差中差均为描述性诊断。供体内log2(归一化计数+1)差；每批n3，不作P/FDR。不建立或投射全基因响应签名。
""",encoding="utf8")
panel=["SPP1","LGALS3","SOCS3","STAT3"];meta=pd.read_csv(d/"GSE261243_SAMPLE_AUDIT.tsv",sep="\t")
values=[];effects=[];sizes=[];checks=[];manifest=[]
for bg,donors in [("LPS",["D562","D563","D565"]),("TNF",["D1","D2","D3"])]:
 file=d/f"sources/GSE261243_RNAseq_{bg}_filtered_raw_counts.txt.gz";f=pd.read_csv(file,sep="\t");cols=f.columns[2:].tolist();A=f.iloc[:,2:].to_numpy(float);assert f.gene_name.is_unique and f.gene_id.is_unique
 mappings=[]
 for c in cols:
  lookup=c if not(bg=="LPS" and c.startswith("Rest_")) else c.replace("Rest_","Rest_None_",1)
  hit=meta[meta.description.str.startswith(lookup+" |")];assert len(hit)==1
  mappings.append({"background":bg,"matrix_column":c,"sample":hit["sample"].iloc[0],"description":hit.description.iloc[0]})
 pd.DataFrame(mappings).to_csv(d/f"SAMPLE_MAPPING_{bg}.tsv",sep="\t",index=False)
 valid=(A>0).all(1);G=np.exp(np.log(A[valid]).mean(1));rat=A[valid]/G[:,None];sf=np.median(rat,axis=0);sf/=np.exp(np.log(sf).mean())
 # Independent logarithmic formulation of median-ratio uses median of ratios, not exp(median log) for even gene count.
 sf2=np.array([np.median(np.exp(np.log(A[valid,j])-np.log(G))) for j in range(A.shape[1])]);sf2/=np.exp(np.log(sf2).mean())
 checks.append({"check":bg+"_median_ratio_formula","passed":np.allclose(sf,sf2,atol=1e-12)})
 uq=np.array([np.quantile(a[a>0],.75) for a in A.T]);uq/=np.exp(np.log(uq).mean())
 for method,fac in [("median_ratio",sf),("upper_quartile",uq)]:
  norm=A/fac[None,:];expr=np.log2(norm+1)
  for c,size in zip(cols,fac):sizes.append({"background":bg,"method":method,"sample":c,"size_factor":size,"genes_for_median_ratio":int(valid.sum()),"provided_genes":len(f)})
  for gene in panel:
   hit=np.flatnonzero(f.gene_name.to_numpy()==gene)
   if not len(hit):
    effects.append(dict(background=bg,method=method,gene=gene,contrast="not_measured",donor="",effect=np.nan,measured=False));continue
   row=int(hit[0])
   for c in cols:
    j=cols.index(c);values.append(dict(background=bg,method=method,gene=gene,gene_id=f.gene_id.iloc[row],sample=c,raw_count=A[row,j],normalized_count=norm[row,j],log2_normalized_plus1=expr[row,j]))
   for donor in donors:
    names={"Rest":f"Rest_{donor}","IL10":f"IL10_{donor}","stim":f"{bg}_{donor}","both":f"IL10.plus.{bg}_{donor}"}
    ix={k:cols.index(n) for k,n in names.items()};e={k:expr[row,j] for k,j in ix.items()}
    for contrast,hi,low in [("IL10_vs_Rest","IL10","Rest"),("IL10_plus_stim_vs_stim","both","stim"),("stim_vs_Rest","stim","Rest"),("IL10_plus_stim_vs_Rest","both","Rest")]:
     diff=e[hi]-e[low];direct=np.log2((norm[row,ix[hi]]+1)/(norm[row,ix[low]]+1));checks.append(dict(check=f"{bg}_{method}_{gene}_{donor}_{contrast}",passed=abs(diff-direct)<1e-12))
     effects.append(dict(background=bg,method=method,gene=gene,contrast=contrast,donor=donor,effect=diff,measured=True,low_raw=A[row,ix[low]],high_raw=A[row,ix[hi]]))
    effects.append(dict(background=bg,method=method,gene=gene,contrast="IL10_by_stim_difference_in_differences",donor=donor,effect=(e["both"]-e["stim"])-(e["IL10"]-e["Rest"]),measured=True,low_raw=np.nan,high_raw=np.nan))
 manifest.append({"path":str(file),"sha256":hashlib.sha256(file.read_bytes()).hexdigest(),"bytes":file.stat().st_size})
E=pd.DataFrame(effects);E.to_csv(d/"TARGET_DONOR_EFFECTS.tsv",sep="\t",index=False);pd.DataFrame(values).to_csv(d/"TARGET_SAMPLE_VALUES.tsv",sep="\t",index=False);pd.DataFrame(sizes).to_csv(d/"NORMALIZATION_FACTORS.tsv",sep="\t",index=False)
summ=[]
for (bg,method,gene,contrast),p in E[E.measured].groupby(["background","method","gene","contrast"]):
 v=p.effect.to_numpy();summ.append(dict(background=bg,method=method,gene=gene,contrast=contrast,n_donors=len(v),positive=int(sum(v>0)),negative=int(sum(v<0)),median_delta=np.median(v),min_delta=min(v),max_delta=max(v)))
S=pd.DataFrame(summ);S.to_csv(d/"TARGET_EFFECT_SUMMARY.tsv",sep="\t",index=False)
pd.DataFrame(checks).to_csv(d/"VALIDATION.tsv",sep="\t",index=False);assert all(x["passed"] for x in checks)
for file in (d/"sources").iterdir():
 if file.suffix in [".soft",".gz"]:manifest.append(dict(path=str(file),sha256=hashlib.sha256(file.read_bytes()).hexdigest(),bytes=file.stat().st_size))
(d/"INPUT_MANIFEST.json").write_text(json.dumps(manifest,indent=2))
(d/"COMPUTE_STATUS.json").write_text(json.dumps(dict(status="COMPLETE",checks=len(checks),P_values_generated=0,FDR_values_generated=0),indent=2))
print(S[(S.method=="median_ratio")&(S.contrast.isin(["IL10_vs_Rest","IL10_plus_stim_vs_stim"]))].to_csv(sep="\t",index=False))
print("SPP1_SAMPLES\n"+pd.DataFrame(values).query("gene=='SPP1' and method=='median_ratio'").to_csv(sep="\t",index=False))
