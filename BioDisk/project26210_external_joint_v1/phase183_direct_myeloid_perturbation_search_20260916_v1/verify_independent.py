from pathlib import Path
import pandas as pd,numpy as np,gzip,json,csv,statistics,hashlib
d=Path("/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase183_direct_myeloid_perturbation_search_20260916_v1")
checks=[]
def ck(n,v,e=0):checks.append({"check":n,"passed":bool(v),"error":float(e)})
V=pd.read_csv(d/"TARGET_SAMPLE_VALUES.tsv",sep="\t");E=pd.read_csv(d/"TARGET_DONOR_EFFECTS.tsv",sep="\t");S=pd.read_csv(d/"TARGET_EFFECT_SUMMARY.tsv",sep="\t")
for bg in ["LPS","TNF"]:
 p=d/f"sources/GSE261243_RNAseq_{bg}_filtered_raw_counts.txt.gz"
 with gzip.open(p,"rt") as f:rows=list(csv.DictReader(f,delimiter="\t"))
 lookup={x["gene_name"]:x for x in rows}
 for _,row in V[V.background==bg].iterrows():
  ck(f"raw_{bg}_{row.method}_{row.gene}_{row['sample']}",float(lookup[row.gene][row["sample"]])==row.raw_count)
 for method in ["median_ratio","upper_quartile"]:
  vv=V[(V.background==bg)&(V.method==method)]
  for _,row in E[(E.background==bg)&(E.method==method)&E.measured].iterrows():
   if row.contrast=="IL10_by_stim_difference_in_differences":continue
   names={"IL10_vs_Rest":("IL10","Rest"),"IL10_plus_stim_vs_stim":("IL10.plus."+bg,bg),"stim_vs_Rest":(bg,"Rest"),"IL10_plus_stim_vs_Rest":("IL10.plus."+bg,"Rest")}
   hi,low=names[row.contrast];a=vv[(vv.gene==row.gene)&(vv["sample"]==hi+"_"+row.donor)].normalized_count.iloc[0];b=vv[(vv.gene==row.gene)&(vv["sample"]==low+"_"+row.donor)].normalized_count.iloc[0]
   err=abs(np.log2((a+1)/(b+1))-row.effect);ck(f"effect_{bg}_{method}_{row.gene}_{row.donor}_{row.contrast}",err<1e-12,err)
for i,row in S.iterrows():
 p=E[(E.background==row.background)&(E.method==row.method)&(E.gene==row.gene)&(E.contrast==row.contrast)]
 vals=p.effect.tolist();ck("summary_"+str(i),len(vals)==3 and statistics.median(vals)==row.median_delta and sum(x>0 for x in vals)==row.positive)
for item in json.loads((d/"INPUT_MANIFEST.json").read_text()):ck("hash_"+item["path"],hashlib.sha256(Path(item["path"]).read_bytes()).hexdigest()==item["sha256"])
pd.DataFrame(checks).to_csv(d/"INDEPENDENT_VALIDATION.tsv",sep="\t",index=False)
# CSV float round-trip medians allow numerical equality tolerance, preserve failure rows if any.
if not all(c["passed"] for c in checks):
 print([c for c in checks if not c["passed"]][:10]);raise RuntimeError("Independent validation failed")
n1=len(pd.read_csv(d/"VALIDATION.tsv",sep="\t"));status={"independent_checks":len(checks),"compute_checks":n1,"total_checks":n1+len(checks),"all_pass":True,"max_error":max(c["error"] for c in checks)}
(d/"VALIDATION_STATUS.json").write_text(json.dumps(status,indent=2));print(status)
