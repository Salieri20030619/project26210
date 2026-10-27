from pathlib import Path
import json,datetime,hashlib,subprocess,sys,os
import numpy as np,pandas as pd
P=Path(__file__).resolve().parents[1];R=P.parent
H=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/03_malignancy_cnv/cnv_annotated_h5ad_v1/part1_scvi_integrated_cell_annotated_cnv_v1.h5ad')
def read(f,**kw):return pd.read_csv(f,sep='\t',float_precision='round_trip',**kw)
def tab(d,n):pd.DataFrame(d).to_csv(P/'tables'/n,sep='\t',index=False)
def dump(f,o):Path(f).write_text(json.dumps(o,ensure_ascii=False,indent=2,default=lambda x:x.item() if hasattr(x,'item') else str(x)))
def sha(f):return hashlib.sha256(Path(f).read_bytes()).hexdigest()
def gate(stage):
 q=subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True)
 with (P/'logs/memory_gate.log').open('a') as f:f.write(stage+'\n'+q.stdout+q.stderr)
 if q.returncode:raise RuntimeError('MEMORY_GATE')
