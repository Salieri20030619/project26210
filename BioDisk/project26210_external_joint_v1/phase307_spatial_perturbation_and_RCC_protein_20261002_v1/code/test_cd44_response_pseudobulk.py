from pathlib import Path
import json,datetime,hashlib
import numpy as np,pandas as pd
from scipy.stats import rankdata,t as td,false_discovery_control
from scipy.linalg import orth
P=Path(__file__).resolve().parents[1];Q=P.parent/'phase305_larger_immune_sequential_20261002_v1'
proto={'created':datetime.datetime.now().astimezone().isoformat(),'version':'patient_pseudobulk_v2','reason':'V1 inherited small-marker-program requirement of detection in5%of individual cells excluded every patient for large402/694gene bulk-derived programs. No V1 association was estimable. The measurement unit here is patient pseudobulk, so use inherited Phase305 single-gene measurement support>=5CD8cells per gene, with at least50%of each directional panel supported. This is a recorded measurement-model revision made after coverage audit and before any evaluable association, not effect/P optimization. V1 full outputs retained.',
'qualification':'CD8>=100; each directional panel >=80%measured in every library AND >=50%of its full genes detected in>=5CD8cells; at least10genes/direction. No genes removed or weights changed. Patient5 Static expected to remain excluded from coverage audit. All other genes, source gates, scores, directions and two models exactly V1.',
'inference':'Per cohort/timepoint/source family4 (Coculture/Static xbatch-stage/batch-stage-CD8depth); n>=10,df>=5 two-sidedt approximation; full rawP/BH; smaller strata noP. Expected negative source-SPP1 to KO-like response; TYMP secondary. Exploratory same-cohort RNA alignment with external CD44 loss, not SPP1-specific perturbation or independent cohort replication.'}
f=P/'resources/CD44_PSEUDOBULK_PROTOCOL_V2.json';assert not f.exists();f.write_text(json.dumps(proto,indent=2))
d=pd.read_csv(P/'tables/CD44_RESPONSE_PATIENT_SCORES.tsv',sep='\t');co=pd.read_csv(P/'tables/CD44_RESPONSE_GENE_COVERAGE.tsv',sep='\t');audit=[];checks=[]
def ck(name,ok):
 checks.append({'check':name,'passed':bool(ok)})
 if not ok:raise AssertionError(name)
for i,r in d.iterrows():
 reasons=[]
 if r.cells<100:reasons.append('CD8_CELL_GATE')
 for direction in ['up','down']:
  g=co[co.analysis_unit.eq(r.analysis_unit)&co.program.eq(r.program)&co.direction.eq(direction)];measured=g.measured_all_libraries.to_numpy(bool);supported=measured&g.detected_cells.ge(5).to_numpy();n=len(g)
  if measured.sum()<.8*n:reasons.append(direction+'_MEASUREMENT')
  if supported.sum()<.5*n:reasons.append(direction+'_PATIENT_DETECTION')
  audit.append({'analysis_unit':r.analysis_unit,'program':r.program,'direction':direction,'genes':n,'measured':int(measured.sum()),'supported_in5cells':int(supported.sum()),'cells':r.cells})
 d.at[i,'eligible']=not reasons;d.at[i,'reason']=';'.join(reasons)
d.to_csv(P/'tables/CD44_PSEUDOBULK_PATIENT_SCORES_V2.tsv',sep='\t',index=False);pd.DataFrame(audit).to_csv(P/'tables/CD44_PSEUDOBULK_QUALIFICATION_V2.tsv',sep='\t',index=False)
old=pd.read_csv(Q/'tables/PATIENT_MEASURES.tsv',sep='\t');source=old[old.lineage.eq('Macrophage')&old.eligible];out=[];pairs=[]
for ds,st in [('OBRADOVIC','Untreated'),('HCRN','Baseline'),('HCRN','Progression')]:
 for src in ['SPP1','TYMP']:
  ss=source[source.dataset.eq(ds)&source.stratum.eq(st)&source.measure.eq(src)][['analysis_unit','value','pT_category']].rename(columns={'value':'source_value'})
  for program in ['Coculture','Static']:
   g=d[d.dataset.eq(ds)&d.stratum.eq(st)&d.program.eq(program)&d.eligible].merge(ss,on='analysis_unit',validate='one_to_one');n=len(g);ck(ds+st+src+program+'_independent_patient_n',n==g.patient.nunique())
   for rr in g.to_dict('records'):pairs.append(rr|{'source':src})
   for model in ['batch_stage_adjusted','batch_stage_CD8depth_adjusted']:
    cols=[np.ones(n)]
    if ds=='OBRADOVIC':cols.extend(pd.get_dummies(g[['technical_batch','pT_category']],drop_first=True,dtype=float).to_numpy().T)
    if model.endswith('CD8depth_adjusted') and n and np.ptp(g.CD8_median_UMI)>0:cols.append(rankdata(g.CD8_median_UMI))
    D=np.column_stack(cols);rank=int(np.linalg.matrix_rank(D)) if n else 0;df=n-rank-1;rho=np.nan;pv=np.nan
    if n>=3:
     xy=np.column_stack([rankdata(g.source_value),rankdata(g.value)]);res=xy-D@np.linalg.lstsq(D,xy,rcond=None)[0]
     if min(np.std(res,axis=0))>1e-9:
      rho=float(np.corrcoef(res.T)[0,1]);U=orth(D);other=xy-U@U.T@xy;ck(ds+st+src+program+model+'_independent_projection',abs(rho-other[:,0]@other[:,1]/np.sqrt((other[:,0]@other[:,0])*(other[:,1]@other[:,1])))<1e-10)
     if n>=10 and df>=5 and np.isfinite(rho):pv=float(2*td.sf(abs(rho)*np.sqrt(df/max(1-rho*rho,1e-15)),df))
    out.append({'dataset':ds,'stratum':st,'source':src,'program':program,'model':model,'patients':n,'df':df,'rho':rho,'P':pv,'family_planned_n':4,'expected_direction':'negative','reason':'' if np.isfinite(pv) else 'SMALL_N_OR_NOT_ESTIMABLE'})
e=pd.DataFrame(out);e['FDR']=np.nan;e['family_tested_n']=0
for key,ids in e.groupby(['dataset','stratum','source']).groups.items():
 valid=[i for i in ids if pd.notna(e.at[i,'P'])];e.loc[list(ids),'family_tested_n']=len(valid)
 if valid:
  pv=e.loc[valid,'P'].to_numpy();q=false_discovery_control(pv);order=np.argsort(pv);manual=np.empty(len(pv));manual[order]=np.minimum(np.minimum.accumulate((pv[order]*len(pv)/np.arange(1,len(pv)+1))[::-1])[::-1],1);ck(str(key)+'_BH',np.allclose(q,manual,atol=1e-12));e.loc[valid,'FDR']=q
e.to_csv(P/'tables/CD44_PSEUDOBULK_ASSOCIATIONS_V2.tsv',sep='\t',index=False);pd.DataFrame(pairs).to_csv(P/'tables/CD44_PSEUDOBULK_PAIRS_V2.tsv',sep='\t',index=False);ck('all24_models',len(e)==24);ck('smalln_noP',e.loc[e.patients<10,'P'].isna().all());(P/'qc/CD44_PSEUDOBULK_VALIDATION_V2.json').write_text(json.dumps({'status':'PASS','checks':checks,'check_n':len(checks),'time':datetime.datetime.now().astimezone().isoformat(),'protocol_sha256':hashlib.sha256(f.read_bytes()).hexdigest()},indent=2));print(e[e.dataset.eq('OBRADOVIC')].to_string(index=False),flush=True)
