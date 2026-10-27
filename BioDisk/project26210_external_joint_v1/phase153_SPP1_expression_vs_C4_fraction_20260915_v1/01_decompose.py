from pathlib import Path
import pandas as pd,numpy as np,json
from scipy.stats import spearmanr
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');R=P/'part2_rebuild_v1'
O=P/'phase153_SPP1_expression_vs_C4_fraction_20260915_v1';O.mkdir(exist_ok=False)
b=pd.read_csv(P/'phase151_TYMP_SPP1_threeTF_role_20260915_v1/TUMOUR_TYMP_MACROPHAGE_MATCHED_PATIENTS.tsv',sep='\t')
m=pd.read_csv(R/'08_regulons/input_v1/CELL_METADATA.tsv.gz',sep='\t',index_col=0,dtype={'candidate':str},low_memory=False)
g=pd.read_csv(R/'07_trajectory/dynamics_v2/CELL_EXPRESSION_PROGRAM_SCORES.tsv.gz',sep='\t',index_col=0,usecols=['cell_id','gene:SPP1'])
m=m[(~m.candidate.isin(['1','Monocyte']))&(m.source_tissue.str.lower()=='tumor')]
rows=[]
for _,r in b.iterrows():
 z=m[m.curated_patient_id==r.patient];a=g.loc[z.index,'gene:SPP1'];inside=z.candidate=='4';n4=int(inside.sum());n0=int((~inside).sum())
 mu4=a[inside].mean() if n4 else np.nan;mu0=a[~inside].mean() if n0 else np.nan
 p=n4/len(z);v4=float(a[inside].sum()/len(z));v0=float(a[~inside].sum()/len(z))
 rows.append(dict(patient=r.patient,study=r.study,TYMP=r.z_TYMP_RNA,n_macro=len(z),n_C4=n4,n_nonC4=n0,C4_fraction=p,SPP1_all=a.mean(),SPP1_C4_mean=mu4,SPP1_nonC4_mean=mu0,C4_weighted_component=v4,nonC4_weighted_component=v0,C4_share_expression_index=v4/a.mean(),C4_share_SPP1_detected=int((a[inside]>0).sum())/int((a>0).sum())))
d=pd.DataFrame(rows)
d.to_csv(O/'PATIENT_SPP1_C4_DECOMPOSITION.tsv',sep='\t',index=False)
(O/'STATUS.json').write_text(json.dumps(dict(status='DECOMPOSITION_COMPLETE',patients=len(d),identity_max_error=float(np.max(np.abs(d.SPP1_all-d.C4_weighted_component-d.nonC4_weighted_component))),figures_created=False,model_rerun=False),indent=2))
print('DECOMPOSITION_MEDIANS');print(d.drop(columns=['patient','study']).median().to_string())
for col in ['SPP1_all','C4_fraction','SPP1_C4_mean','SPP1_nonC4_mean','C4_weighted_component','nonC4_weighted_component']:
 z=d.dropna(subset=[col])
 if col=='SPP1_C4_mean':z=z[z.n_C4>=20]
 if col=='SPP1_nonC4_mean':z=z[z.n_nonC4>=50]
 r,p=spearmanr(z.TYMP,z[col]);print(col,'n',len(z),'rho',r,'P_descriptive_unadjusted',p)
print('IDENTITY_MAX_ERROR',np.max(np.abs(d.SPP1_all-d.C4_weighted_component-d.nonC4_weighted_component)))
print('C4_NONC4_SIZES');print(d.groupby('study').agg(patients=('patient','size'),C4_fraction_median=('C4_fraction','median'),C4_expression_share_median=('C4_share_expression_index','median')).to_string())

