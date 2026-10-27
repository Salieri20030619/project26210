from pathlib import Path
import pandas as pd,json,subprocess,hashlib
P=Path(__file__).resolve().parents[1];R=P.parent;H=R.parent
q=subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True)
(P/'logs/inventory_memory.log').write_text(q.stdout+q.stderr);assert q.returncode==0
rows=[]
for ds in ['GSE264586','GSE294109']:
 f=H/'project26210_external_joint_v1_phase67_direct_program_replication_v1/inputs_v1/metadata_v2'/f'{ds}_phase67_metadata_v1.tsv.gz'
 d=pd.read_csv(f,sep='\t')
 print(ds, d.shape,list(d.columns),flush=True)
 for c in d:
  if any(k in c.lower() for k in ['patient','sample','broad','fine','tissue','treat','response']):
   print(c,d[c].value_counts().head(35).to_dict(),flush=True)
 d.head(3).to_csv(P/'resources'/f'{ds}_schema_example.tsv',sep='\t',index=False)
 for (b,fi),g in d.groupby(['phase7_broad_annotation','phase7_fine_annotation']):
  rows.append(dict(dataset=ds,broad=b,fine=fi,cells=len(g)))
pd.DataFrame(rows).to_csv(P/'tables/EXTERNAL_LABEL_INVENTORY.tsv',sep='\t',index=False)
