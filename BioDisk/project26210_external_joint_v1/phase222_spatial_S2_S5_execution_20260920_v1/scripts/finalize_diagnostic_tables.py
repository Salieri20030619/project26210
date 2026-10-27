from pathlib import Path
import json,shutil
import numpy as np,pandas as pd
P=Path(__file__).resolve().parents[1];archive=P/'audit/diagnostic_table_corrections';archive.mkdir(exist_ok=True)
f=P/'S5/LEAVE_ONE_PATIENT_COHORT.tsv';shutil.copy2(f,archive/f.name);pat=pd.read_csv(P/'ALL_PATIENT_EFFECTS.tsv',sep='\t');rows=[]
for ep,g in pat[pat.model=='primary'].groupby('endpoint'):
 pp=g[g.effect.notna()]
 for kind,field in [('patient','patient'),('cohort','cohort')]:
  for val in pp[field].unique():
   a=pp.loc[pp[field]!=val,'effect'].to_numpy();rows.append(dict(endpoint=ep,drop_type=kind,dropped=val,n=len(a),median=float(np.median(a)) if len(a) else np.nan,positive=int((a>0).sum()),negative=int((a<0).sum())))
new=pd.DataFrame(rows);new.to_csv(f,sep='\t',index=False)
for ep,g in new[new.drop_type=='patient'].groupby('endpoint'):assert len(g)==pat[(pat.endpoint==ep)&(pat.model=='primary')].effect.notna().sum()
f=P/'S4/TLS_ELIGIBILITY.tsv';shutil.copy2(f,archive/f.name);t=pd.read_csv(f,sep='\t');t['distance_computed']=t.S1_candidate
for col in ['TLS_source','peri_source','far_source','unreachable_or_unknown_source']:t.loc[~t.distance_computed,col]=np.nan
t.to_csv(f,sep='\t',index=False)
(P/'audit/DIAGNOSTIC_TABLE_CORRECTION.json').write_text(json.dumps({'S5':'Leave-one only includes patients/cohorts with an estimable primary effect; NA-only records are not deletion units. Main effects, P and FDR unchanged.','S4':'Flag distance_computed; uncomputed coverage fields are NA rather than artificial zeros/infinities. Eligibility and effects unchanged.'},indent=2))
print(new.groupby(['endpoint','drop_type']).size().to_string())
print('LOO_RANGES',new[new.drop_type=='patient'].groupby('endpoint')['median'].agg(['min','max']).to_string())
print('LEAVE_COHORT',new[new.drop_type=='cohort'].to_string(index=False))
