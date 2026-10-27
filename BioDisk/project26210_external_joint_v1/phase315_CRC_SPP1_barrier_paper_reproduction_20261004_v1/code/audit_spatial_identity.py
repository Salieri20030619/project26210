from pathlib import Path
import pandas as pd,json,subprocess,numpy as np,hashlib
P=Path(__file__).resolve().parents[1]
subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],check=True)
x=pd.ExcelFile(P/'downloads/source.xlsx')
selected=[s for s in x.sheet_names if s.startswith('Figure 6') or s.startswith('Supplementary Figure 9') or s.startswith('Supplementary Figure 8')]
d={s:pd.read_excel(x,sheet_name=s) for s in selected}
print('GOs',d['Figure 6l'].to_string(index=False))
for s in ['Figure 6b','Figure 6e','Supplementary Figure 9b','Supplementary Figure 9e','Supplementary Figure 9g','Supplementary Figure 9h']:
 a=d[s];col='DefineTypes' if 'DefineTypes' in a else 'ident';print('LABELS',s,a[col].value_counts().to_dict())
for s in ['Figure 6j','Figure 6k','Figure 6h']:
 a=d[s]
 print('CHECK',s)
 for t in ['Supplementary Figure 9g','Supplementary Figure 9h','Supplementary Figure 9k']:
  b=d[t];ac=[c for c in a if 'Macrophage' in c][0];bc=[c for c in b if 'macrophage' in c.lower()][0]
  if len(a)==len(b):print(t,'equal_order_max_abs_diff',float(np.max(abs(a[ac].to_numpy()-b[bc].to_numpy()))))
print('METADATA')
w=pd.ExcelFile(P/'downloads/41467_2022_29366_MOESM7_ESM.xlsx')
for s in w.sheet_names:
 a=pd.read_excel(w,sheet_name=s,nrows=2);print(s,list(a.columns))
for s,a in d.items():a.to_csv(P/'tables'/('SOURCE_'+s.replace(' ','_')+'.tsv'),sep='\t',index=False)
