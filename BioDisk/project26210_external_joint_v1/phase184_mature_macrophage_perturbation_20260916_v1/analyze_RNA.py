from pathlib import Path
import pandas as pd,numpy as np,tarfile,gzip,io,json,hashlib
d=Path("/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase184_mature_macrophage_perturbation_20260916_v1")
mapping={"SPP1":"6696","LGALS3":"3958","SOCS3":"9021","STAT3":"6774"}
(d/"GENE_ID_MAP.json").write_text(json.dumps({g:{"entrez":i,"url":"https://www.ncbi.nlm.nih.gov/gene/"+i} for g,i in mapping.items()},indent=2))
(d/"RNA_ANALYSIS_PLAN_BEFORE_TARGETS_CN.md").write_text("""# 两个RNAseq干预的最终表达口径
公开附件为Entrez编号FPKM，不是原始counts；读取原FPKM，以log2(FPKM+1)作描述。主要SPP1，辅助LGALS3/SOCS3/STAT3。
GSE181250比较anti-IL10+LPS与同型IgG+LPS各3独立制备；GSE180897分别在未刺激和LPS背景比较siSTAT3与siCNT各3独立制备。公开编号未明确跨组供体身份，不能凭编号进行供体配对。报告组中位差、全部样本范围及处理/对照是否完全分离，不报告3/3配对变化，不作P/FDR。
GSE46903具备BC58/76/78严格72h配对，需另取得完整归一化表达和平台映射；不读取损坏或未完整下载矩阵作为结果。
""",encoding="utf8")
rows=[];man=[]
for acc in ["GSE181250","GSE180897"]:
 md=pd.read_csv(d/(acc+"_METADATA.tsv"),sep="\t").set_index("GSM")
 p=d/"sources"/(acc+"_RAW.tar")
 with tarfile.open(p) as t:
  for member in t.getmembers():
   if not member.isfile():continue
   gsm=member.name.split("_")[0];assert gsm in md.index
   txt=gzip.decompress(t.extractfile(member).read()).decode()
   vals={}
   for ln in txt.splitlines():
    c=ln.strip().split("\t")
    if len(c)>=2:assert c[0] not in vals;vals[c[0]]=float(c[1])
   descr=str(md.loc[gsm,"characteristics_ch1"])
   if acc=="GSE181250":group="anti_IL10_LPS" if "Blocking" in descr else "IgG_LPS"
   else:group=("siSTAT3" if "STAT3 knockdown" in descr else "siCNT")+("_LPS" if "Exposed to LPS" in descr else "_rest")
   for gene,eid in mapping.items():
    fpkm=vals.get(eid,np.nan)
    rows.append(dict(study=acc,GSM=gsm,group=group,gene=gene,entrez=eid,FPKM=fpkm,log2_FPKM_plus1=np.log2(fpkm+1) if np.isfinite(fpkm) else np.nan,provided=eid in vals))
 man.append(dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
R=pd.DataFrame(rows);R.to_csv(d/"RNA_TARGET_SAMPLE_VALUES.tsv",sep="\t",index=False)
out=[]
for acc,tr,co in [("GSE181250","anti_IL10_LPS","IgG_LPS"),("GSE180897","siSTAT3_rest","siCNT_rest"),("GSE180897","siSTAT3_LPS","siCNT_LPS")]:
 for gene in mapping:
  a=R[(R.study==acc)&(R.group==tr)&(R.gene==gene)].log2_FPKM_plus1.dropna().to_numpy()
  b=R[(R.study==acc)&(R.group==co)&(R.gene==gene)].log2_FPKM_plus1.dropna().to_numpy()
  assert len(a)==len(b)==3
  out.append(dict(study=acc,contrast=tr+"_vs_"+co,gene=gene,n_treated=3,n_control=3,paired=False,median_delta=np.median(a)-np.median(b),treated_min=min(a),treated_max=max(a),control_min=min(b),control_max=max(b),all_treated_below_all_controls=max(a)<min(b),all_treated_above_all_controls=min(a)>max(b)))
S=pd.DataFrame(out);S.to_csv(d/"RNA_TARGET_SUMMARY.tsv",sep="\t",index=False)
(d/"RNA_INPUT_MANIFEST.json").write_text(json.dumps(man,indent=2))
print(S.to_csv(sep="\t",index=False));print("SPP1\n"+R[R.gene=="SPP1"].to_csv(sep="\t",index=False))
