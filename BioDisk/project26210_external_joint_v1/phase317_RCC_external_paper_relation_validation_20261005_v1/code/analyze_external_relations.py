from pathlib import Path
import json,subprocess,datetime
import pandas as pd,numpy as np
from scipy.stats import rankdata,binomtest
P=Path(__file__).resolve().parents[1]
reg=pd.read_csv(P/'resources/EXTERNAL_SECTION_REGISTER.tsv',sep='\t');rows=[]
def boolean(v):return v.astype(str).str.lower().isin(['true','1'])
def corr(x,y,z=None):
 x=np.asarray(x,float);y=np.asarray(y,float);ok=np.isfinite(x)&np.isfinite(y)
 if z is not None:ok&=np.isfinite(z)
 n=int(ok.sum())
 if n<100:return np.nan,n,'fewer_than_100_spots'
 x=rankdata(x[ok]);y=rankdata(y[ok])
 if z is not None:
  zz=rankdata(np.asarray(z)[ok]);D=np.column_stack([np.ones(n),zz]);x=x-D@np.linalg.lstsq(D,x,rcond=None)[0];y=y-D@np.linalg.lstsq(D,y,rcond=None)[0]
 else:x=x-x.mean();y=y-y.mean()
 den=np.linalg.norm(x)*np.linalg.norm(y)
 if den<1e-10:return np.nan,n,'constant_variable'
 return float(x@y/den),n,'OK'
pairs=[('FAP_CAF25','SPP1_MAC25'),('FAP_CAF25_disjoint','SPP1_MAC25_disjoint'),('g_FAP','g_SPP1'),('FAP_CAF25','mapped_SPP1_share')]
for a in reg.itertuples():
 subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],check=True)
 d=P/'inputs'/a.sample;m=pd.read_csv(d/'METADATA.tsv',sep='\t',float_precision='round_trip');scores=pd.read_csv(d/'MODULE_SCORES.tsv',sep='\t',float_precision='round_trip');assert m.spot_id.tolist()==scores.spot_id.tolist();z=m.merge(scores,on='spot_id',validate='one_to_one')
 base={'sample':a.sample,'dataset':a.dataset,'patient':a.patient,'preparation':a.preparation}
 for family,pp in [('coupling',pairs),('process',[(x,y) for x in ['FAP_CAF25_disjoint','SPP1_MAC25_disjoint'] for y in ['ECM_org_heldout','Collagen_org_heldout','TGFb_response_heldout']])]:
  for x,y in pp:
   for model in ['unadjusted','depth_adjusted']:
    ef,n,why=corr(z[x],z[y],None if model=='unadjusted' else z.log1p_UMI)
    rows.append(dict(**base,family=family,endpoint=x+'__'+y,model=model,effect_type='spearman_rho',effect=ef,n_spots=n,reason=why,expected_direction=1))
 for definition in ['morphology_supported','expression_sensitivity']:
  valid=boolean(z.boundary_input_valid);inside=valid&boolean(z['inside_'+definition]);outside=valid&~boolean(z['inside_'+definition]);ni=int(inside.sum());no=int(outside.sum())
  for gene in ['CD3D','CD8A']:
   for typ in ['delta_mean_log_expression','delta_detection_fraction']:
    ok=ni>=50 and no>=50;column='g_'+gene if typ=='delta_mean_log_expression' else 'det_'+gene;v=z[column].to_numpy(float)
    reason='OK' if ok else 'fewer_than_50_spots_in_compartment'
    if ok and not np.isfinite(v[valid]).all():ok=False;reason='gene_not_measured'
    if ok and np.max(z.loc[valid,'g_'+gene])==0:ok=False;reason='all_zero_gene'
    effect=float(v[inside].mean()-v[outside].mean()) if ok else np.nan
    rows.append(dict(**base,family='immune',endpoint=gene+'__'+typ,model=definition,effect_type=typ,effect=effect,n_spots=int(valid.sum()),n_inside=ni,n_outside=no,reason=reason,expected_direction=-1,mean_inside=float(v[inside].mean()) if ni and np.isfinite(v[inside]).all() else np.nan,mean_outside=float(v[outside].mean()) if no and np.isfinite(v[outside]).all() else np.nan))
 print('ANALYZED',a.sample,flush=True)
sections=pd.DataFrame(rows);sections.to_csv(P/'tables/SECTION_EFFECTS_ALL.tsv',sep='\t',index=False)
subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],check=True)
# Distinct summary strata are not added together as extra patients.
allsec=sections.copy();allsec['stratum']='overall';extra=sections[sections.dataset.eq('GSE175540')].copy();extra['stratum']='preparation_'+extra.preparation;allsec=pd.concat([allsec,extra],ignore_index=True)
keys=['dataset','stratum','family','endpoint','model','effect_type','expected_direction'];patients=[]
for kk,g in allsec.groupby(keys+['patient'],dropna=False,sort=True):
 v=g.effect.dropna().to_numpy();eff=float(np.tanh(np.arctanh(np.clip(v,-1+1e-12,1-1e-12)).mean())) if len(v) and kk[5]=='spearman_rho' else float(v.mean()) if len(v) else np.nan
 patients.append(dict(zip(keys+['patient'],kk),effect=eff,n_sections=int(g.effect.notna().sum()),sections=';'.join(g.loc[g.effect.notna(),'sample']),ineligible_sections=';'.join(g.loc[g.effect.isna(),'sample'])))
pa=pd.DataFrame(patients);pa.to_csv(P/'tables/PATIENT_EFFECTS_ALL.tsv',sep='\t',index=False);summary=[]
for kk,g in pa.groupby(keys,dropna=False,sort=True):
 v=g.effect.dropna().to_numpy();n=len(v);positive=int((v>0).sum());negative=int((v<0).sum());zero=n-positive-negative;nt=positive+negative
 pv=float(binomtest(positive,nt,.5,alternative='two-sided').pvalue) if n>=10 and nt>0 else (1. if n>=10 else np.nan)
 summary.append(dict(zip(keys,kk),n_patients=n,n_positive=positive,n_negative=negative,n_zero=zero,sign_test_nonzero_n=nt,median_effect=float(np.median(v)) if n else np.nan,p_raw=pv,FDR=np.nan,family_planned_size={'coupling':8,'process':12,'immune':8}[kk[2]],evidence_type='exploratory_patient_direction_test' if n>=10 else 'descriptive_n_lt10'))
su=pd.DataFrame(summary)
for kk,ix in su.groupby(['dataset','stratum','family']).groups.items():
 ix=list(ix);size=int(su.loc[ix[0],'family_planned_size']);assert len(ix)==size
 valid=su.loc[ix,'p_raw'].notna();vals=su.loc[ix,'p_raw'].fillna(1).to_numpy();o=np.argsort(vals);f=np.empty(size);f[o]=np.minimum(1,np.minimum.accumulate((vals[o]*size/np.arange(1,size+1))[::-1])[::-1]);su.loc[np.array(ix)[valid], 'FDR']=f[valid]
su['nominal']=su.p_raw.lt(.05);su['FDR_significant']=su.FDR.lt(.05);su['expected_median_direction']=su.median_effect*su.expected_direction>0
su.to_csv(P/'tables/COHORT_SUMMARY_ALL.tsv',sep='\t',index=False)
(P/'qc/ANALYSIS_COMPLETE.json').write_text(json.dumps({'status':'COMPLETE','time':datetime.datetime.now().astimezone().isoformat(),'section_rows':len(sections),'patient_rows':len(pa),'summary_rows':len(su),'n_sections':reg.sample.nunique() if False else reg['sample'].nunique(),'n_patients':reg[['dataset','patient']].drop_duplicates().shape[0],'reference_model_calibration_pass':False,'evidence':'exploratory published-marker transfer'},indent=2))
print(su[su.stratum.eq('overall')][['dataset','family','endpoint','model','n_patients','n_positive','n_negative','median_effect','p_raw','FDR']].to_json(orient='records'),flush=True)
