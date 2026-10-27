from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
r=Path("/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1");d=r/"phase179_STAT3_external_activity_validation_20260916_v1"
genes=set(pd.read_csv(r/"part2_rebuild_v1/08_regulons/input_v1/GENES.tsv",sep="\t").gene)
pub=Path("/mnt/BioDisk/wangmingyuan/project26210_phase98_public_paired_metabolomics_v1/downloads_v1")
for c in ["RC18","RC20","CPTAC","CPTAC_val"]:genes &=set(pd.read_csv(pub/("matched_tpm_"+c+".csv"),nrows=0).columns)
rows=[]
for res in ["collectri","dorothea"]:
 f=pd.read_csv(d/"sources"/(res+".tsv"),sep="\t")
 for _,v in f.iterrows():
  gene=v.target_genesymbol;direction=bool(v.is_stimulation)!=bool(v.is_inhibition)
  why="keep"
  if v.source_genesymbol!="STAT3":why="not_STAT3"
  elif not direction:why="ambiguous_or_unsigned"
  elif gene in ["TYMP","SPP1","STAT3"]:why="outcome_exposure_self_removed"
  elif gene not in genes:why="not_in_common_measured_universe"
  elif pd.isna(v.references):why="no_reference"
  elif res=="dorothea" and not set(str(v.dorothea_level).split(";"))&set("ABC"):why="confidence"
  rows.append({"resource":res,"gene":gene,"sign":1 if v.is_stimulation else -1,"decision":why,"references":v.references})
a=pd.DataFrame(rows);a.to_csv(d/"TARGET_AUDIT.tsv",sep="\t",index=False)
k=a[a.decision=="keep"].drop_duplicates(["resource","gene"]);assert not k.duplicated(["resource","gene"]).any()
k.to_csv(d/"FROZEN_STAT3_TARGETS.tsv",sep="\t",index=False)
print(k.groupby(["resource","sign"]).size().to_string())
print(a.groupby(["resource","decision"]).size().to_string())
assert len(k[k.resource=="collectri"])>=20 and len(k[k.resource=="dorothea"])>=10
(d/"TARGET_FREEZE.json").write_text(json.dumps({"sha256":hashlib.sha256((d/"FROZEN_STAT3_TARGETS.tsv").read_bytes()).hexdigest(),"common_universe":len(genes),"primary":"collectri","sensitivity":"dorotheaABC","excluded":["SPP1","TYMP","STAT3"],"score":"signed mean of gene-wise z(log expression); gene-wise rank-z sensitivity","biological_interpretation":"transcriptional activity proxy, not phosphoprotein","sign_policy":"stimulation XOR inhibition; ambiguous excluded","source_overlap":"two resources share prior evidence; not independent replication"},indent=2))
