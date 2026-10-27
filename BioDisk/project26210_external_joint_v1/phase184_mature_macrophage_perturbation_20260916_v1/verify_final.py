from pathlib import Path
import pandas as pd,numpy as np,json,csv,gzip,io,hashlib
d=Path("/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase184_mature_macrophage_perturbation_20260916_v1")
V=pd.read_csv(d/"ARRAY_TARGET_PROBE_VALUES.tsv",sep="\t");E=pd.read_csv(d/"ARRAY_DONOR_EFFECTS.tsv",sep="\t");S=pd.read_csv(d/"ARRAY_EFFECT_SUMMARY.tsv",sep="\t");ck=[]
for file in (d/"sources").glob("GSM*.full.soft"):
 gsm=file.name.split(".")[0];text=file.read_text();chunk=text.split("!sample_table_begin\n",1)[1].split("!sample_table_end",1)[0]
 rr={r["ID_REF"]:r for r in csv.DictReader(io.StringIO(chunk),delimiter="\t")}
 for _,row in V[V.GSM==gsm].iterrows():
  ck.append(dict(check=gsm+"_"+row.probe,passed=float(rr[row.probe]["VALUE"])==row.log2_quantile_normalized))
for i,row in E.iterrows():
 p=V[(V.donor==row.donor)&(V.gene==row.gene)];a=p[p.group=="IL10"].log2_quantile_normalized;b=p[p.group=="control"].log2_quantile_normalized
 ck.append(dict(check="effect_"+str(i),passed=abs(a.median()-b.median()-row.effect)<1e-12))
for i,row in S.iterrows():
 x=E[E.gene==row.gene].effect
 ck.append(dict(check="summary_"+str(i),passed=len(x)==3 and abs(x.median()-row.median_delta)<1e-12 and sum(x<0)==row.negative))
pr=pd.read_csv(d/"ARRAY_PROBE_COVERAGE.tsv",sep="\t")
ck.append(dict(check="single_unreleased_SOCS3_probe",passed=set(pr[~pr.present].probe)=={"ILMN_2156250"}))
assert all(x["passed"] for x in ck)
pd.DataFrame(ck).to_csv(d/"ARRAY_INDEPENDENT_VALIDATION.tsv",sep="\t",index=False)
a=json.loads((d/"ARRAY_STATUS.json").read_text());b=json.loads((d/"RNA_VALIDATION_STATUS.json").read_text())
status={"RNA_checks":b["checks"],"array_internal_checks":a["checks"],"array_independent_checks":len(ck),"total_checks":b["checks"]+a["checks"]+len(ck),"all_pass":True}
(d/"VALIDATION_STATUS.json").write_text(json.dumps(status,indent=2))
man=[]
for p in [d/"sources/GPL6947_complete.annot.gz",d/"sources/GSE181250_RAW.tar",d/"sources/GSE180897_RAW.tar"]+list((d/"sources").glob("GSM*.full.soft")):
 h=hashlib.sha256(p.read_bytes()).hexdigest();man.append(dict(path=str(p),bytes=p.stat().st_size,sha256=h))
(d/"COMPLETE_INPUT_MANIFEST.json").write_text(json.dumps(man,indent=2));print(status)
