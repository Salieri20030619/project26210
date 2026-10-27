from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='4'
import json,datetime,subprocess,hashlib,zipfile
import numpy as np,pandas as pd,openpyxl,anndata as ad
from scipy.stats import rankdata,t as tdist,ttest_rel,false_discovery_control
P=Path(__file__).resolve().parents[1];R=P.parent;Q=R/'phase305_larger_immune_sequential_20261002_v1'
protocol={
 'created':datetime.datetime.now().astimezone().isoformat(),
 'scope':'External source-data reanalysis plus literature-directed exploratory recipient endpoints in existing OBRADOVIC/HCRN single cells. Source workbook values inspected before protocol; external results are reproduction, not preregistered discovery.',
 'singlecell_endpoints':{'PDCD1_ENTPD1_double_detected_fraction':['PDCD1','ENTPD1'],'TOX_detected_fraction':['TOX'],'IFNG_GZMB_effector_mean_logCP10K':['IFNG','GZMB']},
 'endpoint_basis':'CIR-24-1146 Fig5D/E; protein phenotypes translated to exploratory RNA readouts. RNA detection is not flow-protein positivity, definitive exhaustion, antigen specificity, or direct function.',
 'populations':'Reuse frozen Phase305 CD8/mature-macrophage identities, patients and treatment strata. No new clustering, no expression-based patient selection, no cohort pooling.',
 'qualification':'CD8 >=100 and all endpoint genes measured in every contributing library. Fractions use count>0 with zeros retained, no minimum positive-cell requirement. Two-gene effector endpoint additionally uses inherited program gate: >=1 of2 genes detected in >=5% CD8. Source TYMP/SPP1 eligibility inherited unchanged.',
 'models':['batch_stage_adjusted','batch_stage_CD8depth_adjusted'],
 'adjustment':'OBRADOVIC author batch1/2 and pT1b/3a; sensitivity adds rank median CD8 UMI. HCRN separate Baseline/Progression, sensitivity adds rank median CD8 UMI.',
 'inference':'Partial Spearman with two-sided t approximation df=n-rank(D)-1, only n>=10 and df>=5; smaller strata effects only. No cell-level P. BH within cohort/stratum/source: 3 endpoints x2 models, planned6. Report all directions, all rawP/FDR and all exclusions.',
 'external_spatial':'Trehan Fig6F-O: reproduce paired t tests and effects using source rows paired as deposited, no region-as-n. Nine distance hypotheses G-O excluding density F form a9 endpoint family, density separate1. Source L/M have14 rows versus caption15: exploratory descriptive only for these discrepant endpoints; do not infer missing patient identities.',
 'external_perturbation':'Trehan Fig4I-M has2 deposited scores per arm, donor mapping not supplied; group mean differences and fractional attenuation only, no P or paired replicate inference. Do not treat normalized scores as fold changes. Fig5A-D KO6/WT7 show all8 contrasts descriptively, retain authorP separately; Fig5E source9/group vs caption5/group flagged, no tumor-weight inference.',
 'boundaries':'External liver/mouse effects are not ccRCC validation or TYMP perturbation. Existing Phase305 broad CD8 negative results retained. None of these establishes SPP1-driven collagen geometry or physical exclusion.'}
fp=P/'resources/DOWNSTREAM_PROTOCOL.json'
if fp.exists():
 old=json.loads(fp.read_text());assert {k:v for k,v in old.items() if k!='created'}=={k:v for k,v in protocol.items() if k!='created'}
else:fp.write_text(json.dumps(protocol,ensure_ascii=False,indent=2))
assert zipfile.ZipFile(P/'downloads/Trehan_SourceData.xlsx').testzip() is None
w=openpyxl.load_workbook(P/'downloads/Trehan_SourceData.xlsx',read_only=True,data_only=True)
pert=[];raw=[]
for panel,row,ep in [('I',74,'T_cell_exhaustion'),('J',78,'JAK_STAT'),('K',82,'MAPK'),('L',86,'NFkB'),('M',90,'TLR')]:
 vals={g:np.array([w['Figure 4'].cell(i,c).value for i in [row,row+1]],float) for g,c in [('IgG',3),('Spp1',4),('Spp1_aCD44',5)]}
 for g,v in vals.items():
  for k,x in enumerate(v):raw.append({'panel':'4'+panel,'endpoint':ep,'group':g,'source_row':row+k,'score':x,'biological_id':'not_provided'})
 means={g:float(v.mean()) for g,v in vals.items()};inc=means['Spp1']-means['IgG'];res=means['Spp1_aCD44']-means['Spp1']
 pert.append({'panel':'4'+panel,'endpoint':ep,'deposited_values_per_arm':2,**means,'stimulation_difference':inc,'CD44_blockade_difference':res,'attenuation_fraction':-res/inc,'remaining_vs_IgG':means['Spp1_aCD44']-means['IgG'],'P':np.nan,'FDR':np.nan,'reason':'biological_replicates_unresolved;descriptive'})
pd.DataFrame(raw).to_csv(P/'tables/EXTERNAL_STIMULATION_SOURCE_VALUES.tsv',sep='\t',index=False);pd.DataFrame(pert).to_csv(P/'tables/EXTERNAL_STIMULATION_EFFECTS.tsv',sep='\t',index=False)
sp=[];sr=[]
endpoints=[('F','Macrophage_density',3,3,17,15),('G','Macro_to_CD8_distance',7,3,17,15),('H','Macro_to_Ki67_CD8_distance',11,3,17,15),('I','Macro_to_CD45RA_CD8_distance',15,3,17,15),('J','Macro_to_CD45RO_CD8_distance',19,3,17,15),('K','Macro_to_CD39_CD8_distance',23,3,17,15),('L','Macro_to_NK_distance',3,21,34,15),('M','CD8_to_NK_distance',7,21,34,15),('N','Macro_to_Treg_distance',11,21,34,14),('O','CD8_to_Treg_distance',15,21,34,14)]
for panel,ep,c,lo,hi,ncaption in endpoints:
 xy=np.array([[w['Figure 6'].cell(i,c).value,w['Figure 6'].cell(i,c+1).value] for i in range(lo,hi+1)],float);assert np.isfinite(xy).all();x,y=xy.T;diff=x-y;n=len(x);ok=n==ncaption
 for i,(a,b) in enumerate(xy):sr.append({'panel':'6'+panel,'endpoint':ep,'source_row':lo+i,'patient_row_index':i+1,'alphaSMA_positive':a,'alphaSMA_negative':b})
 pv=float(ttest_rel(x,y).pvalue) if ok else np.nan
 sp.append({'panel':'6'+panel,'endpoint':ep,'source_pairs':n,'caption_patients':ncaption,'mean_positive':x.mean(),'mean_negative':y.mean(),'mean_paired_difference':diff.mean(),'median_paired_difference':np.median(diff),'median_ratio_positive_over_negative':np.median(x/y),'positive_direction':int((diff>0).sum()),'negative_direction':int((diff<0).sum()),'P':pv,'family':'density' if panel=='F' else 'distances_G_to_O','family_planned_n':1 if panel=='F' else 9,'reason':'' if ok else 'SOURCE_CAPTION_N_DISCREPANCY_NO_INFERENCE'})
sp=pd.DataFrame(sp);sp['FDR']=np.nan
for fam,ids in sp.groupby('family').groups.items():
 good=[i for i in ids if pd.notna(sp.at[i,'P'])];m=1 if fam=='density' else 9
 pv=sp.loc[good,'P'].to_numpy();padded=np.r_[pv,np.ones(m-len(pv))];sp.loc[good,'FDR']=false_discovery_control(padded)[:len(pv)]
sp.to_csv(P/'tables/EXTERNAL_SPATIAL_EFFECTS.tsv',sep='\t',index=False);pd.DataFrame(sr).to_csv(P/'tables/EXTERNAL_SPATIAL_SOURCE_VALUES.tsv',sep='\t',index=False)
ko=[];kr=[]
for panel,ep,ks,ws,pp,pn in [('A','CD39_positive',3,10,.0165,.4001),('B','TIM3_positive',20,27,.0464,.9834),('C','CD44pos_CD62Lneg',37,44,.0338,.3264),('D','CD44neg_CD62Lpos',53,60,.0229,.0986)]:
 for col,pop,ap in [(4,'TRP2_positive',pp),(5,'TRP2_negative',pn)]:
  groups={g:np.array([w['Figure 5'].cell(i,col).value for i in range(start,start+n)],float) for g,start,n in [('KO',ks,6),('WT',ws,7)]}
  for g,v in groups.items():
   for k,x in enumerate(v):kr.append({'panel':'5'+panel,'endpoint':ep,'population':pop,'group':g,'source_row':(ks if g=='KO' else ws)+k,'percent':x})
  x=groups['KO'];y=groups['WT'];ko.append({'panel':'5'+panel,'endpoint':ep,'population':pop,'KO_n':len(x),'WT_n':len(y),'KO_mean':x.mean(),'WT_mean':y.mean(),'difference_percentage_points':x.mean()-y.mean(),'relative_mean_change':x.mean()/y.mean()-1,'author_P':ap,'new_P':np.nan,'new_FDR':np.nan,'reason':'small_groups;external_mouse;authorP_not_new_test'})
pd.DataFrame(ko).to_csv(P/'tables/EXTERNAL_KO_EFFECTS.tsv',sep='\t',index=False);pd.DataFrame(kr).to_csv(P/'tables/EXTERNAL_KO_SOURCE_VALUES.tsv',sep='\t',index=False)
print('EXTERNAL_COMPLETE',flush=True)
subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],check=True,stdout=subprocess.DEVNULL)
rows=[];coverage=[]
for ds in ['OBRADOVIC','HCRN']:
 a=ad.read_h5ad(Q/f'inputs/{ds}_counts_qc.h5ad',backed='r');meta=pd.read_csv(Q/f'tables/{ds}_QUALIFIED_CELL_METADATA.tsv.gz',sep='\t',index_col=0);assert a.obs_names.equals(meta.index)
 cells=np.flatnonzero(meta.lineage_CD8.to_numpy(bool));x=a[cells,:].X.tocsr();m=meta.iloc[cells].copy();names=a.var_names.copy();a.file.close();present={}
 for lib in m.library.unique():
  b=ad.read_h5ad(Q/f'inputs/{ds}_{lib}_qc.h5ad',backed='r');present[lib]=set(b.var_names);b.file.close()
 for unit,g in m.groupby('analysis_unit',observed=True):
  sel=m.analysis_unit.eq(unit).to_numpy();xx=x[sel];n=xx.shape[0];avail=set.intersection(*(present[l] for l in g.library.unique()));summ=np.asarray(xx.sum(0)).ravel();pb=np.log1p(1e4*summ/summ.sum());genes=['PDCD1','ENTPD1','TOX','IFNG','GZMB','CD44'];z=xx[:,names.get_indexer(genes)].toarray();det=(z>0).sum(0)
  base={'dataset':ds,'analysis_unit':unit,'patient':g.patient.iloc[0],'stratum':g.stratum.iloc[0],'technical_batch':g.technical_batch.iloc[0],'cells':n,'CD8_median_UMI':float(g.total_counts.median())}
  for j,gene in enumerate(genes):coverage.append(base|{'gene':gene,'measured_all_libraries':gene in avail,'detected_cells':int(det[j]),'detected_fraction':det[j]/n,'pseudobulk_logCP10K':pb[names.get_loc(gene)]})
  values=[np.mean((z[:,0]>0)&(z[:,1]>0)),np.mean(z[:,2]>0),pb[names.get_indexer(['IFNG','GZMB'])].mean()]
  for (ep,gs),value in zip(protocol['singlecell_endpoints'].items(),values):
   reasons=[]
   if n<100:reasons.append('CD8_CELL_GATE')
   if not set(gs)<=avail:reasons.append('GENE_MEASUREMENT_GATE')
   if ep.startswith('IFNG') and (det[3:5]>=.05*n).sum()<1:reasons.append('INHERITED_PROGRAM_DETECTION_GATE')
   rows.append(base|{'endpoint':ep,'value':value,'eligible':not reasons,'reason':';'.join(reasons)})
 print('RECIPIENT_INPUTS',ds,flush=True)
d=pd.DataFrame(rows);d.to_csv(P/'tables/RECIPIENT_PATIENT_ENDPOINTS.tsv',sep='\t',index=False);pd.DataFrame(coverage).to_csv(P/'tables/RECIPIENT_GENE_COVERAGE.tsv',sep='\t',index=False)
old=pd.read_csv(Q/'tables/PATIENT_MEASURES.tsv',sep='\t');source=old[old.lineage.eq('Macrophage')&old.eligible].copy();results=[];pairs=[]
for ds,st in [('OBRADOVIC','Untreated'),('HCRN','Baseline'),('HCRN','Progression')]:
 for src in ['SPP1','TYMP']:
  ss=source[source.dataset.eq(ds)&source.stratum.eq(st)&source.measure.eq(src)][['analysis_unit','value','pT_category']].rename(columns={'value':'source_value'})
  for ep in protocol['singlecell_endpoints']:
   gg=d[d.dataset.eq(ds)&d.stratum.eq(st)&d.endpoint.eq(ep)&d.eligible].merge(ss,on='analysis_unit',validate='one_to_one');n=len(gg);assert gg.patient.nunique()==n
   for model in protocol['models']:
    cols=[np.ones(n)]
    if ds=='OBRADOVIC':
     for cat in ['technical_batch','pT_category']:cols.extend(pd.get_dummies(gg[cat],drop_first=True,dtype=float).to_numpy().T)
    if model.endswith('CD8depth_adjusted') and n and np.ptp(gg.CD8_median_UMI)>0:cols.append(rankdata(gg.CD8_median_UMI))
    D=np.column_stack(cols);rank=int(np.linalg.matrix_rank(D)) if n else 0;df=n-rank-1;rho=np.nan;pval=np.nan
    if n>=3:
     xy=np.column_stack([rankdata(gg.source_value),rankdata(gg.value)]);res=xy-D@np.linalg.lstsq(D,xy,rcond=None)[0]
     if min(np.std(res,axis=0))>1e-9:rho=float(np.corrcoef(res.T)[0,1])
     if n>=10 and df>=5 and np.isfinite(rho):pval=float(2*tdist.sf(abs(rho)*np.sqrt(df/max(1-rho*rho,1e-15)),df))
    results.append({'dataset':ds,'stratum':st,'source':src,'endpoint':ep,'model':model,'patients':n,'design_rank':rank,'df':df,'rho':rho,'P':pval,'family_planned_n':6,'reason':'' if np.isfinite(pval) else 'SMALL_N_OR_NOT_ESTIMABLE'})
   for rr in gg.to_dict('records'):pairs.append(rr|{'source':src})
res=pd.DataFrame(results);res['FDR']=np.nan;res['family_tested_n']=0
for key,ids in res.groupby(['dataset','stratum','source']).groups.items():
 good=[i for i in ids if pd.notna(res.at[i,'P'])];res.loc[list(ids),'family_tested_n']=len(good)
 if good:res.loc[good,'FDR']=false_discovery_control(res.loc[good,'P'].to_numpy())
res.to_csv(P/'tables/RECIPIENT_ASSOCIATIONS_ALL.tsv',sep='\t',index=False);pd.DataFrame(pairs).to_csv(P/'tables/RECIPIENT_ASSOCIATION_PAIRS.tsv',sep='\t',index=False)
(P/'qc/ANALYSIS_COMPLETE.json').write_text(json.dumps({'status':'COMPLETE','time':datetime.datetime.now().astimezone().isoformat(),'association_rows':len(res),'external_spatial_rows':len(sp),'protocol_sha256':hashlib.sha256(fp.read_bytes()).hexdigest()},indent=2))
print('OBRADOVIC_RESULTS',res[res.dataset.eq('OBRADOVIC')].to_string(index=False),flush=True)
print('ANALYSIS_COMPLETE',flush=True)
