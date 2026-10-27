from common import *
import sys,json
from scipy import sparse
from scipy.stats import spearmanr,rankdata
sys.path.insert(0,str(R/'phase157_ordered_bridge_reassessment_20260915_v1'))
from bridge_stats import bridge,bootstrap_ab,rankz,corrp
B=pd.read_csv(R/'phase154_same_library_TYMP_C4_test_20260915_v1/SAME_LIBRARY_PATIENT_METRICS.tsv',sep='\t').set_index('patient')
C=pd.read_csv(R/'phase154_same_library_TYMP_C4_test_20260915_v1/SELECTED_CELL_TWO_GENE_COUNTS.tsv.gz',sep='\t',index_col=0,dtype={'candidate':str})
assert len(B)==16 and B.index.is_unique
m=pd.read_csv(D/'inputs/MACRO_METADATA.tsv.gz',sep='\t',index_col=0,dtype={'candidate':str})
defs={'U':['4','10','15'],'C4':['4'],'C10':['10'],'C15':['15']}
metrics=[]
for pat,b in B.iterrows():
 z=C[(C.sample_id==b.sample_id)&(C.lineage=='Macrophage')]
 assert (z.candidate.astype(str)==m.loc[z.index,'candidate']).all()
 assert len(z)==b.macrophage_cells
 for label,cs in defs.items():
  a=z[z.candidate.isin(cs)];other=z[~z.candidate.isin(cs)]
  raw=float(a.SPP1_counts.sum());tot=float(a.recomputed_library_counts.sum())
  metrics.append(dict(patient=pat,study=b.study,sample_id=b.sample_id,group=label,cells=len(a),other_cells=len(other),fraction=len(a)/len(z),SPP1_RNA_share=raw/z.SPP1_counts.sum() if z.SPP1_counts.sum()>0 else np.nan,SPP1_log2CPM=np.log2(raw/tot*1e6+1) if tot else np.nan,SPP1_detection=float((a.SPP1_counts>0).mean()) if len(a) else np.nan))
M=pd.DataFrame(metrics)
M.to_csv(D/'bridge/SAME_LIBRARY_GROUP_METRICS.tsv',sep='\t',index=False)
wide=M.pivot(index='patient',columns='group',values='fraction').loc[B.index]
assert np.allclose(wide.C4,B.C4_fraction)
rng=np.random.default_rng(262200);nperm=9999
groups=[np.flatnonzero(B.study.to_numpy()==s) for s in sorted(B.study.unique())]
perms=np.tile(np.arange(16),(nperm,1))
for ix in groups:perms[:,ix]=np.array([rng.permutation(ix) for _ in range(nperm)])
direct=[];loo=[]
for label in defs:
 for ex in ['TYMP_pseudobulk','TYMP_mean_log1pCP10K']:
  x=rankdata(B[ex]);xr=x.copy()
  endpoints={'fraction':wide[label].to_numpy(),'whole_macro_SPP1':B.SPP1_pseudobulk.to_numpy(),'whole_macro_SPP1_detection':B.SPP1_detection_fraction.to_numpy()}
  for target,y in endpoints.items():
   y=rankdata(y);yr=y.copy()
   for ix in groups:xr[ix]-=xr[ix].mean();yr[ix]-=yr[ix].mean()
   den=np.sqrt(np.dot(xr,xr)*np.dot(yr,yr))
   rho=float(np.dot(xr,yr)/den) if den else np.nan
   pp=float((np.sum(np.abs(xr[perms]@yr/den)>=abs(rho)-1e-14)+1)/(nperm+1)) if den else np.nan
   rr,p=spearmanr(x,y)
   direct.extend([dict(group=label,exposure=ex,target=target,adjustment='raw',n=16,rho=rr,P=p),dict(group=label,exposure=ex,target=target,adjustment='study_permutation',n=16,rho=rho,P=pp)])
  if ex=='TYMP_pseudobulk':
   for typ,vals in [('patient',B.index),('study',B.study.unique())]:
    for omit in vals:
     keep=(B.index!=omit) if typ=='patient' else (B.study!=omit).to_numpy()
     for adj in ['raw','study']:
      cov=np.ones((keep.sum(),1)) if adj=='raw' else pd.get_dummies(B.loc[keep,'study'],dtype=float).to_numpy()
      rr,pp=corrp(rankz(B.loc[keep,ex]),rankz(wide.loc[keep,label]),cov)
      loo.append(dict(group=label,omission_type=typ,omitted=omit,adjustment=adj,n=int(keep.sum()),rho=rr))
F=add_fdr(pd.DataFrame(direct),['group','exposure','adjustment'])
F.to_csv(D/'bridge/TYMP_DIRECT_RESULTS.tsv',sep='\t',index=False)
pd.DataFrame(loo).to_csv(D/'bridge/TYMP_LEAVE_ONE.tsv',sep='\t',index=False)
# Matched-correlation difference with reranking in every stratified patient bootstrap.
X=pd.read_csv(R/'phase157_ordered_bridge_reassessment_20260915_v1/04_other_X_screen/PATIENT_CANDIDATE_VALUES.tsv',sep='\t',index_col=0).loc[B.index]
preds={'TYMP':B.TYMP_pseudobulk,**{c:X[c] for c in ['ligand:IGFBP3','program:Cholesterol Homeostasis','program:Glycolysis']}}
rng=np.random.default_rng(262201);draw=np.tile(np.arange(16),(2000,1))
for ix in groups:draw[:,ix]=rng.choice(ix,size=(2000,len(ix)),replace=True)
diffs=[]
for name,x in preds.items():
 values=[]
 for ix in draw:
  rx=rankdata(x.to_numpy()[ix]);u=rankdata(wide.U.to_numpy()[ix]);c=rankdata(wide.C4.to_numpy()[ix])
  if np.std(rx)>0 and np.std(u)>0 and np.std(c)>0:values.append(np.corrcoef(rx,u)[0,1]-np.corrcoef(rx,c)[0,1])
 fullu=spearmanr(x,wide.U).statistic;fullc=spearmanr(x,wide.C4).statistic
 diffs.append(dict(predictor=name,n=16,U_rho=fullu,C4_rho=fullc,delta_rho=fullu-fullc,bootstrap_low=np.quantile(values,.025),bootstrap_high=np.quantile(values,.975),valid=len(values)))
pd.DataFrame(diffs).to_csv(D/'bridge/U_MINUS_C4_RHO_BOOTSTRAP.tsv',sep='\t',index=False)
# Retain original screen and decomposition registries; no re-selection by union results.
registries=[
 ('original116',X,pd.read_csv(R/'phase157_ordered_bridge_reassessment_20260915_v1/04_other_X_screen/CANDIDATE_REGISTRY.tsv',sep='\t')),
 ('metabolic_components',pd.read_csv(R/'phase159_metabolic_parts_C4_value_20260915_v1/PATIENT_SENDER_VALUES.tsv',sep='\t',index_col=0).loc[B.index],pd.read_csv(R/'phase159_metabolic_parts_C4_value_20260915_v1/CANDIDATE_REGISTRY.tsv',sep='\t'))]
mp=pd.read_csv(R/'phase157_ordered_bridge_reassessment_20260915_v1/03_MTA_to_new_C4/MATCHED_MTA_PROXY_PATIENTS.tsv',sep='\t').set_index('patient').loc[B.index]
mcols=['direct_MTA_response_RNA_proxy','polyamine_synthesis_core_RNA','polyamine_SAM_supply_extended_RNA','MTAP_clearance_RNA']
registries.append(('MTA_proxies',mp[mcols],pd.DataFrame({'candidate':mcols,'family':'MTA_proxy'})))
results=[];sensitivity=[];receiver=[]
for screen,values,registry in registries:
 assert set(registry.candidate)==set(values.columns)
 for label in ['U','C4']:
  yy=wide[label].to_numpy()
  for row in registry.itertuples():
   name=row.candidate;family=str(row.family)+(('_'+str(row.level)) if hasattr(row,'level') else '')
   for adj in ['raw','study_adjusted']:
    out=bridge(B.TYMP_pseudobulk,values[name],yy,B.study,adj)
    out.update(bootstrap_ab(B.TYMP_pseudobulk,values[name],yy,B.study,adj,seed=262202,B=2000))
    results.append(dict(screen=screen,group=label,candidate=name,family=family,adjustment=adj,**out))
    for covname in ['program:Hypoxia','program:G2-M Checkpoint']:
     cov=np.ones((16,1)) if adj=='raw' else pd.get_dummies(B.study,dtype=float).to_numpy()
     cov=np.column_stack([cov,rankz(X[covname])]);A=rankz(B.TYMP_pseudobulk);V=rankz(values[name]);Y=rankz(yy)
     ar,ap=corrp(A,V,cov);yr,yp=corrp(V,Y,cov);br,bp=corrp(V,Y,np.column_stack([cov,A]))
     sensitivity.append(dict(screen=screen,group=label,candidate=name,family=family,adjustment=adj,covariate=covname,n=16,a_rho=ar,a_P=ap,XY_rho=yr,XY_P=yp,b_rho=br,b_P=bp))
 out=pd.DataFrame([v for v in results if v['screen']==screen])
 print('SCREEN_DONE',screen,len(out),flush=True)
res=pd.DataFrame(results);sen=pd.DataFrame(sensitivity)
for df,keys in [(res,['screen','group','family','adjustment']),(sen,['screen','group','family','adjustment','covariate'])]:
 for edge in ['a','XY','b']:add_fdr(df,keys,p=edge+'_P',name=edge+'_FDR')
res.to_csv(D/'bridge/ALL_FIXED_CANDIDATE_BRIDGES.tsv',sep='\t',index=False)
sen.to_csv(D/'bridge/ALL_FIXED_CANDIDATE_HYPOXIA_CYCLE.tsv',sep='\t',index=False)
assert np.nanmax(abs(res.identity_error))<1e-9
# MTA response score already frozen over each whole compartment; only group membership changes.
score=pd.read_csv(R/'phase157_ordered_bridge_reassessment_20260915_v1/02_inputs/CELL_MTA_RESPONSE_SCORES.tsv.gz',sep='\t',index_col=0)
rows=[]
for scope in ['all','tumor']:
 md=m if scope=='all' else m[m.source_tissue.str.lower()=='tumor']
 for label,cs in defs.items():
  for pat,z in md.groupby('curated_patient_id'):
   yes=z.candidate.isin(cs);a=z.index[yes];b=z.index[~yes]
   if min(len(a),len(b))<20:continue
   rows.append(dict(scope=scope,group=label,patient=pat,study=z.study.iloc[0],effect=score.loc[a,'MTA_response_transcript_score'].mean()-score.loc[b,'MTA_response_transcript_score'].mean()))
pt=pd.DataFrame(rows);ss=[]
for key,z in pt.groupby(['scope','group']):ss.append(dict(zip(['scope','group'],key),**summary(z.effect)))
ss=add_fdr(pd.DataFrame(ss),['scope'])
pt.to_csv(D/'receiver/MTA_RESPONSE_PATIENTS.tsv',sep='\t',index=False);ss.to_csv(D/'receiver/MTA_RESPONSE_SUMMARY.tsv',sep='\t',index=False)
# Paired receiver expression/detection families from the new exact raw pseudobulk.
panels={'MTA_receptor':['ADORA2A','ADORA2B','STAT3','SOCS3','IL10','VEGFA','ARG1'],
'cholesterol':['LDLR','SCARB1','NR1H3','NR1H2','RXRA','ABCA1','ABCG1','CYP27A1'],
'lactate':['SLC16A1','SLC16A3','BSG','HCAR1','GPR132','LDHA','LDHB'],
'IL6_IL10_LGALS3':['IL6','IL6R','IL6ST','IL10','IL10RA','IL10RB','LGALS3','SPP1']}
for label in defs:
 z=np.load(D/'core/contrasts'/f'tumor_{label}.npz',allow_pickle=True);gi={g:i for i,g in enumerate(z['genes'])};eff=z['effects'];a=z['target_detection'];b=z['reference_detection']
 for family,genes in panels.items():
  for gene in genes:
   j=gi[gene];v=eff[:,j]
   receiver.append(dict(group=label,family=family,gene=gene,target_median_detection=float(np.median(a[:,j])),reference_median_detection=float(np.median(b[:,j])),target_patients_ge10pct=int((a[:,j]>=.1).sum()),reference_patients_ge10pct=int((b[:,j]>=.1).sum()),**summary(v)))
rc=add_fdr(pd.DataFrame(receiver),['group','family']);rc.to_csv(D/'receiver/RECEPTOR_GENE_PANELS.tsv',sep='\t',index=False)
(D/'qc/BRIDGE_STATUS.json').write_text(json.dumps({'status':'PASS','same_library_patients':16,'original_C4_fractions_reproduced':True,'candidate_models':len(res),'conditional_identity_max':float(abs(res.identity_error).max()),'bootstrap_rho_differences':len(diffs),'no_client_data':True,'no_figures':True},indent=2))
print('BRIDGE_COMPLETE',flush=True)
