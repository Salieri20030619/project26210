from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='2'
import json,subprocess
import numpy as np,pandas as pd
from scipy.stats import rankdata,t,false_discovery_control
P=Path(__file__).resolve().parents[1];pr=json.loads((P/'resources/ANALYSIS_PROTOCOL.json').read_text());src=pd.read_csv(P/'tables/SOURCE_PATIENT_MEASURES.tsv',sep='\t');rec=pd.read_csv(P/'tables/CD8_PATIENT_MEASURES.tsv',sep='\t')
q=subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True);(P/'logs/memory_gate_analyze.log').write_text(q.stdout+q.stderr);assert q.returncode==0
results=[];pairs=[];sens=[];loo=[]
def evaluate(g,model,depth=False):
 n=len(g)
 if n==0:return np.nan,0,-1
 D=[np.ones(n)]
 if model=='batch_stage_adjusted' and g.dataset.iloc[0]=='OBRADOVIC':
  for c in ['technical_batch','pT_category']:D.extend(pd.get_dummies(g[c].astype(str),drop_first=True,dtype=float).to_numpy().T)
 if depth:
  for c in ['source_median_UMI','receiver_median_UMI']:
   if g[c].isna().any():return np.nan,0,-1
   D.append(rankdata(np.round(g[c].to_numpy(float),12)))
 D=np.column_stack(D);dr=int(np.linalg.matrix_rank(D));df=n-dr-1
 if n<3 or df<1:return np.nan,dr,df
 v=np.column_stack([rankdata(np.round(g[c].to_numpy(float),12)) for c in ['source_value','receiver_value']]);z=v-D@np.linalg.lstsq(D,v,rcond=None)[0]
 if (np.std(z,axis=0)<1e-9).any():return np.nan,dr,df
 return float(np.corrcoef(z.T)[0,1]),dr,df
for ds,scope in [('OBRADOVIC','Untreated'),('HCRN','Baseline'),('HCRN','Progression')]:
 for source in pr['sources']:
  a=src[(src.dataset==ds)&(src.stratum==scope)&(src.source==source)].copy();a=a.rename(columns={'value':'source_value','median_UMI':'source_median_UMI','eligible_primary':'source_primary','eligible_strict':'source_strict','cells':'source_cells'})
  for endpoint in list(pr['programs'])+['CD8_share_of_CD8_plus_macrophages']:
   b=rec[(rec.dataset==ds)&(rec.stratum==scope)&(rec.endpoint==endpoint)].copy();b=b.rename(columns={'value':'receiver_value','eligible_primary':'receiver_primary','eligible_strict':'receiver_strict'})
   columns=['analysis_unit','patient','endpoint','family','receiver_value','receiver_primary','receiver_strict','CD8_cells','macrophage_cells','receiver_median_UMI','reason']
   g=a.merge(b[columns],on=['analysis_unit','patient'],validate='one_to_one');g=g[g.source_primary&g.receiver_primary].copy();assert not g.patient.duplicated().any()
   for model in pr['models']:
    rho,dr,df=evaluate(g,model);n=len(g);pv=float(2*t.sf(abs(rho)*np.sqrt(df/max(1-rho*rho,1e-15)),df)) if n>=10 and df>=5 and np.isfinite(rho) else np.nan
    family='CD8_Hallmark' if endpoint in pr['programs'] else 'capture_composition';eid=f'{ds}|{scope}|{source}|{endpoint}|{model}'
    row=dict(edge_id=eid,dataset=ds,stratum=scope,source=source,endpoint=endpoint,family=family,model=model,patients=n,rho=rho,design_rank=dr,df=df,P=pv,status='inferential_exploratory' if np.isfinite(pv) else 'descriptive' if np.isfinite(rho) else 'not_evaluable',adjustment_applied=bool(model=='batch_stage_adjusted' and ds=='OBRADOVIC'),family_planned_n=300 if family=='CD8_Hallmark' else 6)
    results.append(row)
    for z in g.to_dict('records'):pairs.append({k:row[k] for k in ['edge_id','dataset','stratum','source','endpoint','model']}|{k:z[k] for k in ['analysis_unit','patient','technical_batch','pT_category','source_value','receiver_value','source_cells','CD8_cells','macrophage_cells','source_median_UMI','receiver_median_UMI','source_strict','receiver_strict']})
    for name,sub,dep in [('strict_qualification',g[g.source_strict&g.receiver_strict],False),('source_receiver_depth_adjusted',g,True)]:
     rr,rank,dd=evaluate(sub,model,dep);sens.append(dict(edge_id=eid,sensitivity=name,patients=len(sub),rho=rr,design_rank=rank,df=dd,P=np.nan))
    if n>=4:
     for patient in g.patient:
      sub=g[g.patient!=patient];rr,rank,dd=evaluate(sub,model);loo.append(dict(edge_id=eid,left_out=patient,patients=len(sub),rho=rr,design_rank=rank,df=dd))
 print('ANALYZED',ds,scope,flush=True)
e=pd.DataFrame(results);e['FDR']=np.nan;e['family_tested_n']=0
for key,ids in e.groupby(['dataset','stratum','family']).groups.items():
 valid=[i for i in ids if np.isfinite(e.at[i,'P'])];e.loc[list(ids),'family_tested_n']=len(valid)
 if valid:e.loc[valid,'FDR']=false_discovery_control(e.loc[valid,'P'].to_numpy())
l=pd.DataFrame(loo);ss=pd.DataFrame(sens);details=[]
for row in e.to_dict('records'):
 a=l[l.edge_id==row['edge_id']].rho.dropna();d=ss[ss.edge_id==row['edge_id']]
 details.append(dict(edge_id=row['edge_id'],loo_n=len(a),loo_same_direction=int(((a*row['rho'])>0).sum()) if np.isfinite(row['rho']) else 0,loo_min=a.min() if len(a) else np.nan,loo_max=a.max() if len(a) else np.nan,strict_patients=int(d[d.sensitivity=='strict_qualification'].patients.iloc[0]),strict_rho=d[d.sensitivity=='strict_qualification'].rho.iloc[0],depth_rho=d[d.sensitivity=='source_receiver_depth_adjusted'].rho.iloc[0]))
e=e.merge(pd.DataFrame(details),on='edge_id',validate='one_to_one');e['nominal_candidate']=e.P<.05;e['FDR_significant']=e.FDR<.05
e.to_csv(P/'tables/ALL_ASSOCIATIONS.tsv',sep='\t',index=False);pd.DataFrame(pairs).to_csv(P/'tables/ALL_PATIENT_PAIRS.tsv.gz',sep='\t',index=False);ss.to_csv(P/'tables/ALL_SENSITIVITIES.tsv',sep='\t',index=False);l.to_csv(P/'tables/LEAVE_ONE_PATIENT.tsv',sep='\t',index=False)
c=e[(e.dataset=='OBRADOVIC')&(e.model=='batch_stage_adjusted')&e.nominal_candidate].sort_values(['family','P']);c.to_csv(P/'tables/ADJUSTED_NOMINAL_CANDIDATES.tsv',sep='\t',index=False)
(P/'qc/ANALYSIS_COMPLETE.json').write_text(json.dumps({'status':'COMPLETE','associations':len(e),'patient_pairs':len(pairs),'sensitivities':len(ss),'leave_one_rows':len(l),'adjusted_nominal_candidates':len(c),'adjusted_FDR_significant':int(c.FDR_significant.sum())},indent=2));print('SUMMARY',len(e),len(c),'FDR',int(c.FDR_significant.sum()),flush=True);print(c[['source','endpoint','patients','rho','P','FDR','loo_same_direction','loo_n','strict_patients','strict_rho','depth_rho']].to_string(index=False))
