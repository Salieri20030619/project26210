from pathlib import Path
import pandas as pd,numpy as np,json,subprocess,datetime
from scipy.stats import rankdata
P=Path(__file__).resolve().parents[1]
subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],check=True)
reg=pd.read_csv(P/'resources/EXTERNAL_SECTION_REGISTER.tsv',sep='\t');sec=pd.read_csv(P/'tables/SECTION_EFFECTS_ALL.tsv',sep='\t',float_precision='round_trip');pa=pd.read_csv(P/'tables/PATIENT_EFFECTS_ALL.tsv',sep='\t',float_precision='round_trip');rows=[]
for a in reg.itertuples():
 m=pd.read_csv(P/'inputs'/a.sample/'METADATA.tsv',sep='\t',float_precision='round_trip')
 for model in ['unadjusted','depth_adjusted']:
  x=rankdata(np.round(m.g_FAP.to_numpy(),12));y=rankdata(np.round(m.g_SPP1.to_numpy(),12))
  if model=='depth_adjusted':
   z=rankdata(m.log1p_UMI);D=np.column_stack([np.ones(len(z)),z]);x=x-D@np.linalg.lstsq(D,x,rcond=None)[0];y=y-D@np.linalg.lstsq(D,y,rcond=None)[0]
  ef=float(np.corrcoef(x,y)[0,1]);old=float(sec[(sec['sample']==a.sample)&(sec.endpoint=='g_FAP__g_SPP1')&(sec.model==model)].effect.iloc[0]);rows.append({'sample':a.sample,'dataset':a.dataset,'patient':a.patient,'preparation':a.preparation,'model':model,'original_rho':old,'round12_rho':ef,'abs_delta':abs(ef-old),'direction_changed':bool(np.sign(ef)!=np.sign(old))})
df=pd.DataFrame(rows);df.to_csv(P/'qc/DIRECT_GENE_ROUND12_SECTION_AUDIT.tsv',sep='\t',index=False);pat=[]
for a in pa[pa.endpoint.eq('g_FAP__g_SPP1')].itertuples():
 d=df[(df.dataset==a.dataset)&(df.patient==a.patient)&(df.model==a.model)]
 if a.stratum!='overall':d=d[d.preparation.eq(a.stratum.removeprefix('preparation_'))]
 v=np.tanh(np.arctanh(d.round12_rho).mean());pat.append({'dataset':a.dataset,'stratum':a.stratum,'patient':a.patient,'model':a.model,'original_rho':a.effect,'round12_rho':v,'abs_delta':abs(v-a.effect),'direction_changed':bool(np.sign(v)!=np.sign(a.effect))})
pd.DataFrame(pat).to_csv(P/'qc/DIRECT_GENE_ROUND12_PATIENT_AUDIT.tsv',sep='\t',index=False)
o={'status':'PASS' if not any(a['direction_changed'] for a in pat) else 'DIRECTION_CHANGED','time':datetime.datetime.now().astimezone().isoformat(),'section_comparisons':len(rows),'patient_stratum_comparisons':len(pat),'section_direction_changes':int(df.direction_changed.sum()),'patient_direction_changes':sum(a['direction_changed'] for a in pat),'max_section_rho_difference':float(df.abs_delta.max()),'max_patient_rho_difference':max(a['abs_delta'] for a in pat),'P_FDR_changes':0 if not any(a['direction_changed'] for a in pat) else 'requires_reassessment','reason':'Patient P is exact sign test; all patient directions and all n unchanged imply all raw P and BH unchanged. Original results retained, no rounding adopted as new model.'}
(P/'qc/FLOAT_TIE_SENSITIVITY.json').write_text(json.dumps(o,indent=2));print(json.dumps(o,indent=2))
