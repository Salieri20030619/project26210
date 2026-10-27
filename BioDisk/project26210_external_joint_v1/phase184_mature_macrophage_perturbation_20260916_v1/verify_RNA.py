from pathlib import Path
import pandas as pd,numpy as np,tarfile,gzip,csv,json,statistics
d=Path("/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase184_mature_macrophage_perturbation_20260916_v1")
R=pd.read_csv(d/"RNA_TARGET_SAMPLE_VALUES.tsv",sep="\t");S=pd.read_csv(d/"RNA_TARGET_SUMMARY.tsv",sep="\t");checks=[]
for acc in ["GSE181250","GSE180897"]:
 with tarfile.open(d/"sources"/(acc+"_RAW.tar")) as t:
  for mm in t.getmembers():
   if not mm.isfile():continue
   gsm=mm.name.split("_")[0];pairs=list(csv.reader(gzip.decompress(t.extractfile(mm).read()).decode().splitlines(),delimiter="\t"));vals={a[0]:float(a[1]) for a in pairs if len(a)>1}
   for _,row in R[R.GSM==gsm].iterrows():
    x=vals[str(row.entrez)];checks.append(dict(check=gsm+"_"+row.gene,passed=x==row.FPKM and abs(np.log2(x+1)-row.log2_FPKM_plus1)<1e-12))
for i,row in S.iterrows():
 tr,co=row.contrast.split("_vs_");a=R[(R.study==row.study)&(R.group==tr)&(R.gene==row.gene)].log2_FPKM_plus1.to_numpy();b=R[(R.study==row.study)&(R.group==co)&(R.gene==row.gene)].log2_FPKM_plus1.to_numpy()
 checks.append(dict(check="summary_"+str(i),passed=len(a)==len(b)==3 and abs(statistics.median(a)-statistics.median(b)-row.median_delta)<1e-12))
pd.DataFrame(checks).to_csv(d/"RNA_INDEPENDENT_VALIDATION.tsv",sep="\t",index=False);assert all(c["passed"] for c in checks)
(d/"RNA_VALIDATION_STATUS.json").write_text(json.dumps({"checks":len(checks),"all_pass":True},indent=2));print(len(checks))
