from pathlib import Path
import gzip,io,pandas as pd,numpy as np,json,hashlib,re
d=Path("/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase184_mature_macrophage_perturbation_20260916_v1")
p=d/"sources/GPL6947_complete.annot.gz";txt=gzip.decompress(p.read_bytes()).decode();lines=txt.splitlines();i=next(i for i,l in enumerate(lines) if l.startswith("ID\t"))
table="\n".join(l for l in lines[i:] if not l.startswith(("!","#")) and l.strip())
ann=pd.read_csv(io.StringIO(table),sep="\t",dtype=str).fillna("");assert ann.ID.is_unique
genes={"SPP1":"6696","LGALS3":"3958","SOCS3":"9021","STAT3":"6774"}
selected=[]
for gene,eid in genes.items():
 hits=ann[(ann["Gene symbol"]==gene)&(ann["Gene ID"]==eid)]
 for _,r in hits.iterrows():selected.append({"gene":gene,"entrez":eid,"probe":r.ID})
A=pd.DataFrame(selected);A.to_csv(d/"ARRAY_PROBE_MAPPING_ALL_ANNOTATED.tsv",sep="\t",index=False)
presence=[]
for file in sorted((d/"sources").glob("GSM*.full.soft")):
 ss=file.read_text();assert "!sample_table_end" in ss
 ff=pd.read_csv(io.StringIO(ss.split("!sample_table_begin\n",1)[1].split("!sample_table_end",1)[0]),sep="\t")
 for probe in A.probe:presence.append({"GSM":file.name.split(".")[0],"probe":probe,"present":probe in set(ff.ID_REF)})
pr=pd.DataFrame(presence);pr.to_csv(d/"ARRAY_PROBE_COVERAGE.tsv",sep="\t",index=False)
keep=pr.groupby("probe").present.all();A=A[A.probe.isin(keep[keep].index)]
A.to_csv(d/"ARRAY_PROBE_MAPPING.tsv",sep="\t",index=False);assert set(A.gene)==set(genes)
pairs={"BC58":("GSM1140525","GSM1140514"),"BC76":("GSM1140530","GSM1140523"),"BC78":("GSM1140486","GSM1140532")}
metadata=pd.read_csv(d/"GSE46903_SAMPLE_METADATA.tsv",sep="\t").set_index("GSM");values=[];frames={};checks=[]
for donor,(tr,co) in pairs.items():
 for group,gsm in [("IL10",tr),("control",co)]:
  info=metadata.loc[gsm,"characteristics_ch1"]
  assert f"donor: {donor}" in info and "time: 72h" in info and "initial differentiation: GM-CSF" in info
  assert ("activation stimuli: IL10" if group=="IL10" else "activation stimuli: con") in info
  text=(d/f"sources/{gsm}.full.soft").read_text();assert "!sample_table_end" in text
  tt=text.split("!sample_table_begin\n",1)[1].split("!sample_table_end",1)[0];f=pd.read_csv(io.StringIO(tt),sep="\t").set_index("ID_REF");assert f.index.is_unique
  frames[gsm]=f
  for _,row in A.iterrows():
   assert row.probe in f.index
   q=f.loc[row.probe];values.append(dict(donor=donor,group=group,GSM=gsm,gene=row.gene,probe=row.probe,log2_quantile_normalized=float(q.iloc[0]),technical_detection_p=float(q.iloc[1])))
V=pd.DataFrame(values);V.to_csv(d/"ARRAY_TARGET_PROBE_VALUES.tsv",sep="\t",index=False)
eff=[];probeeff=[]
for donor,(tr,co) in pairs.items():
 for gene,aa in A.groupby("gene"):
  t=V[(V.GSM==tr)&(V.gene==gene)].set_index("probe");c=V[(V.GSM==co)&(V.gene==gene)].set_index("probe").loc[t.index]
  delta=t.log2_quantile_normalized-c.log2_quantile_normalized
  eff.append(dict(donor=donor,gene=gene,probes=len(t),effect=t.log2_quantile_normalized.median()-c.log2_quantile_normalized.median(),positive_probes=int(sum(delta>0)),negative_probes=int(sum(delta<0)),control_probes_detected=int(sum(c.technical_detection_p<.05)),treated_probes_detected=int(sum(t.technical_detection_p<.05))))
  for probe in t.index:
   x=t.loc[probe,"log2_quantile_normalized"]-c.loc[probe,"log2_quantile_normalized"];probeeff.append(dict(donor=donor,gene=gene,probe=probe,effect=x))
   direct=frames[tr].loc[probe,"VALUE"]-frames[co].loc[probe,"VALUE"];checks.append(dict(check=donor+"_"+probe,passed=abs(x-direct)<1e-12))
E=pd.DataFrame(eff);E.to_csv(d/"ARRAY_DONOR_EFFECTS.tsv",sep="\t",index=False);pd.DataFrame(probeeff).to_csv(d/"ARRAY_PROBE_DONOR_EFFECTS.tsv",sep="\t",index=False)
s=[]
for gene,e in E.groupby("gene"):
 s.append(dict(gene=gene,n_donors=len(e),median_delta=e.effect.median(),positive=int(sum(e.effect>0)),negative=int(sum(e.effect<0)),min_delta=e.effect.min(),max_delta=e.effect.max(),probes=int(e.probes.iloc[0])))
S=pd.DataFrame(s);S.to_csv(d/"ARRAY_EFFECT_SUMMARY.tsv",sep="\t",index=False)
pd.DataFrame(checks).to_csv(d/"ARRAY_VALIDATION.tsv",sep="\t",index=False);assert all(x["passed"] for x in checks)
(d/"ARRAY_STATUS.json").write_text(json.dumps({"status":"COMPLETE","donors":3,"libraries":6,"platform_annot_rows":len(ann),"checks":len(checks),"new_P_values":0,"annotation_sha256":hashlib.sha256(p.read_bytes()).hexdigest()},indent=2))
print("PROBES\n"+A.to_csv(sep="\t",index=False));print("SUMMARY\n"+S.to_csv(sep="\t",index=False));print("DONORS\n"+E.to_csv(sep="\t",index=False))
