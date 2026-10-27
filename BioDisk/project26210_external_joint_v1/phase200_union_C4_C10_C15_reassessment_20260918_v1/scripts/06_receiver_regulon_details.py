from common import *
import json
from scipy import sparse
m=pd.read_csv(D/'inputs/MACRO_METADATA.tsv.gz',sep='\t',index_col=0,dtype={'candidate':str})
genes=pd.read_csv(D/'inputs/GENES.tsv',sep='\t').gene.tolist();gi={g:i for i,g in enumerate(genes)}
raw=sparse.load_npz(D/'inputs/MACRO_COUNTS.npz').tocsr()
tot=m.library_counts.to_numpy();norm=raw.multiply(1e4/tot[:,None]).tocsr();norm.data=np.log1p(norm.data)
defs={'U':['4','10','15'],'C4':['4'],'C10':['10'],'C15':['15']}
tum=m.source_tissue.str.lower().eq('tumor')
features=['SPP1','TYMP','LGALS3','IL6','IL6R','IL6ST','IL10','IL10RA','IL10RB','SLC16A3','LDHA','JUN','JUNB','CEBPD','STAT3','SRM','SMS','MTAP']
E=pd.DataFrame(norm[:,[gi[g] for g in features]].toarray(),index=m.index,columns=features)
# Per-patient effects, no cell-level P values.
effects=[];elig=[];coverage=[]
edge_defs=[('LGALS3','SPP1','LGALS3'),('IL6','SPP1','IL6'),('JUN','SPP1','threeTF_RNA'),('JUNB','SPP1','threeTF_RNA'),('CEBPD','SPP1','threeTF_RNA'),('TYMP','SPP1','TYMP_RNA')]
for label,cs in dict(defs,WholeMac=list(m.candidate.unique())).items():
 for pat,z in m[tum & m.candidate.isin(cs)].groupby('curated_patient_id'):
  e=E.loc[z.index];depth=np.log1p(z.library_counts.to_numpy());ng=z.n_genes_raw.to_numpy()
  for gene in ['IL6R','IL6ST','IL10RA','IL10RB','SPP1']:
   coverage.append(dict(group=label,patient=pat,feature=gene,cells=len(z),detected=int((e[gene]>0).sum()),fraction=float((e[gene]>0).mean())))
  for name,a,b in [('IL6R_AND_IL6ST','IL6R','IL6ST'),('IL10RA_AND_IL10RB','IL10RA','IL10RB')]:
   hit=(e[a]>0)&(e[b]>0);coverage.append(dict(group=label,patient=pat,feature=name,cells=len(z),detected=int(hit.sum()),fraction=float(hit.mean())))
  for x,y,family in edge_defs:
   ok=len(z)>=50 and (e[x]>0).sum()>=10 and (e[y]>0).sum()>=10
   elig.append(dict(group=label,patient=pat,edge=x+'_'+y,cells=len(z),detected_x=int((e[x]>0).sum()),detected_y=int((e[y]>0).sum()),eligible=bool(ok)))
   if not ok:continue
   adjustments={'raw':[],'depth':[depth],'depth_genes':[depth,ng]}
   if label=='WholeMac':
    adjustments['depth_C4']=[depth,z.candidate.eq('4').to_numpy(int)]
    adjustments['depth_U']=[depth,z.candidate.isin(defs['U']).to_numpy(int)]
   for adj,cv in adjustments.items():
    rr=partial(np.round(e[x],12),np.round(e[y],12),cv)
    effects.append(dict(module='RNA',group=label,patient=pat,study=z.study.iloc[0],edge=x+'_'+y,family=family,adjustment=adj,rho=rr))
# Exact Phase161 same-library source expression/hypoxia cache; re-aggregate only masks.
c=pd.read_csv(R/'phase161_SLC16A3_cell_origin_20260915_v1/CELL_TARGET_EXPRESSION.tsv.gz',sep='\t')
h=np.load(R/'phase161_SLC16A3_cell_origin_20260915_v1/HYPoxia_CELL_INPUT.npz',allow_pickle=True)
assert np.array_equal(c.cell_id.to_numpy(),h['cell_id'])
candidate=m.candidate.reindex(c.cell_id).to_numpy()
for label,cs in dict(defs,WholeMac=list(m.candidate.unique())).items():
 for pat,z in c[np.isin(candidate,cs)].groupby('curated_patient_id'):
  for x,y in [('SLC16A3','SPP1'),('LDHA','SPP1')]:
   if len(z)<50 or (z['counts:'+x]>0).sum()<10 or (z['counts:'+y]>0).sum()<10:continue
   hx=h['expression'][z.index];sd=hx.std(0,ddof=1);ok=sd>0
   score=((hx[:,ok]-hx[:,ok].mean(0))/sd[ok]).mean(1)
   for adj,cv in [('raw',[]),('depth',[np.log1p(z.total_counts)]),('depth_hypoxia',[np.log1p(z.total_counts),score])]:
    effects.append(dict(module='same_library_lactate',group=label,patient=pat,study=z.study.iloc[0],edge=x+'_'+y,family='lactate',adjustment=adj,rho=partial(z['expr:'+x],z['expr:'+y],cv)))
# Frozen Phase179 activity score; reuse it rather than train/redefine target sets.
cache=pd.read_csv(R/'phase179_STAT3_external_activity_validation_20260916_v1/run_v2/SC_PRIMARY_SCORE_CACHE.tsv.gz',sep='\t')
activity=[]
for resource,ss in cache.groupby('resource'):
 for pat,z in ss.groupby('patient'):
  z=z.set_index('cell_id');zz=m.loc[z.index];dep=z.depth.to_numpy()
  assert np.allclose(z.SPP1,E.loc[z.index,'SPP1'])
  for label,cs in defs.items():
   hit=zz.candidate.isin(cs).to_numpy()
   if min(hit.sum(),(~hit).sum())>=20:
    activity.append(dict(resource=resource,group=label,patient=pat,study=zz.study.iloc[0],effect=float(z.score.to_numpy()[hit].mean()-z.score.to_numpy()[~hit].mean())))
  for label,cs in dict(defs,WholeMac=list(m.candidate.unique())).items():
   hit=zz.candidate.isin(cs).to_numpy()
   if hit.sum()<50 or (z.SPP1.to_numpy()[hit]>0).sum()<10:continue
   for adj,cv in [('raw',[]),('depth',[dep[hit]]),('depth_genes',[dep[hit],zz.n_genes_raw.to_numpy()[hit]])]:
    effects.append(dict(module='STAT3_frozen_'+resource,group=label,patient=pat,study=zz.study.iloc[0],edge='STAT3_SPP1',family='STAT3_SPP1',adjustment=adj,rho=partial(z.score.to_numpy()[hit],z.SPP1.to_numpy()[hit],cv)))
# LGALS3-independent STAT3 scores for the previously studied three-edge mechanism family.
K=pd.read_csv(R/'phase182_IL10_LGALS3_receiver_audit_20260916_v1/STAT3_TARGETS_WITHOUT_LGALS3.tsv',sep='\t')
target=sorted(set(K.gene));assert not set(target)&{'SPP1','TYMP','STAT3','LGALS3'}
targetE=norm[:,[gi[g] for g in target]].toarray()
for pat,z in m[tum].groupby('curated_patient_id'):
 if len(z)<50:continue
 ix=m.index.get_indexer(z.index);mat=targetE[ix];baseE=E.loc[z.index];dep=np.log1p(z.library_counts.to_numpy())
 for resource,kk in K.groupby('resource'):
  j=[target.index(g) for g in kk.gene];a=mat[:,j]
  for scoring in ['zmean','rankmean']:
   b=rankdata(a,axis=0) if scoring=='rankmean' else a.copy();sd=b.std(0,ddof=1);ok=sd>1e-12
   if ok.sum()<(20 if resource=='collectri' else 10):continue
   score=np.round((((b[:,ok]-b[:,ok].mean(0))/sd[ok])*kk.sign.to_numpy()[ok]).mean(1),12)
   for label,cs in dict(defs,WholeMac=list(m.candidate.unique())).items():
    hit=z.candidate.isin(cs).to_numpy()
    if hit.sum()<50:continue
    vals={'STAT3':score[hit],'LGALS3':baseE.LGALS3.to_numpy()[hit],'SPP1':baseE.SPP1.to_numpy()[hit]}
    for xx,yy in [('STAT3','LGALS3'),('STAT3','SPP1'),('LGALS3','SPP1')]:
     if (vals[yy]>0).sum()<10 or (xx=='LGALS3' and (vals[xx]>0).sum()<10):continue
     for adj,cv in [('depth',[dep[hit]]),('depth_genes',[dep[hit],z.n_genes_raw.to_numpy()[hit]])]:
      effects.append(dict(module='STAT3_noLGALS3_'+resource+'_'+scoring,group=label,patient=pat,study=z.study.iloc[0],edge=xx+'_'+yy,family='three_mechanism_edges',adjustment=adj,rho=partial(vals[xx],vals[yy],cv)))
ef=pd.DataFrame(effects);ef.to_csv(D/'receiver/PATIENT_WITHIN_CELL_EFFECTS.tsv.gz',sep='\t',index=False)
pd.DataFrame(elig).to_csv(D/'receiver/EXPRESSION_EDGE_ELIGIBILITY.tsv',sep='\t',index=False)
pd.DataFrame(coverage).to_csv(D/'receiver/JOINT_RECEPTOR_DETECTION.tsv',sep='\t',index=False)
summaries=[];leave=[]
keys=['module','group','family','edge','adjustment']
for key,z in ef.groupby(keys):
 v=z.rho.dropna().to_numpy();s=summary(v)
 if len(v):
  rng=np.random.default_rng(262203);boots=np.median(v[rng.integers(0,len(v),(2000,len(v)))],axis=1)
  s.update(CI_low=float(np.quantile(boots,.025)),CI_high=float(np.quantile(boots,.975)))
 summaries.append(dict(zip(keys,key),**s))
 for st in z.study.unique():leave.append(dict(zip(keys,key),omitted_study=st,**summary(z.loc[z.study!=st,'rho'])))
ss=add_fdr(pd.DataFrame(summaries),['module','group','family','adjustment'])
ss.to_csv(D/'receiver/WITHIN_CELL_SUMMARY.tsv',sep='\t',index=False)
pd.DataFrame(leave).to_csv(D/'receiver/WITHIN_CELL_LEAVE_STUDY.tsv',sep='\t',index=False)
# Shared-patient rho effects: U versus C4, preserving both estimates.
diff=[]
for key,z in ef[ef.group.isin(['U','C4'])].groupby(['module','family','edge','adjustment']):
 w=z.pivot(index='patient',columns='group',values='rho').dropna()
 if not {'U','C4'}<=set(w):continue
 diff.append(dict(zip(['module','family','edge','adjustment'],key),U_median=float(w.U.median()),C4_median=float(w.C4.median()),**summary(w.U-w.C4)))
if diff:add_fdr(pd.DataFrame(diff),['module','family','adjustment']).to_csv(D/'receiver/SHARED_PATIENT_U_MINUS_C4.tsv',sep='\t',index=False)
ac=pd.DataFrame(activity);ac.to_csv(D/'receiver/STAT3_ACTIVITY_PATIENTS.tsv',sep='\t',index=False)
asum=[dict(resource=res,group=label,**summary(z.effect)) for (res,label),z in ac.groupby(['resource','group'])]
add_fdr(pd.DataFrame(asum),['resource']).to_csv(D/'receiver/STAT3_ACTIVITY_SUMMARY.tsv',sep='\t',index=False)
# All regulon contrasts for constituents, union and direct shared-patient difference.
auc=pd.read_csv(R/'part2_rebuild_v1/08_regulons/network_v2_20260915/AUC_ALL_CELLS.csv',index_col=0).loc[m.index]
pairs=[]
for scope in ['all','tumor']:
 md=m if scope=='all' else m[tum]
 for pat,z in md.groupby('curated_patient_id'):
  a=auc.loc[z.index]
  for label,cs in defs.items():
   hit=z.candidate.isin(cs).to_numpy()
   if min(hit.sum(),(~hit).sum())<20:continue
   v=a.iloc[hit].mean()-a.iloc[~hit].mean()
   for reg,delta in v.items():pairs.append(dict(scope=scope,patient=pat,study=z.study.iloc[0],group=label,regulon=reg,TF=reg.split('(')[0],effect=float(delta)))
pa=pd.DataFrame(pairs);pa.to_csv(D/'regulons/ALL_GROUP_PATIENT_CONTRASTS.tsv.gz',sep='\t',index=False)
rows=[]
for key,z in pa.groupby(['scope','group','regulon','TF']):rows.append(dict(zip(['scope','group','regulon','TF'],key),**summary(z.effect)))
rs=add_fdr(pd.DataFrame(rows),['scope','group']);rs.to_csv(D/'regulons/ALL_GROUP_AUC_SUMMARY.tsv',sep='\t',index=False)
rd=[]
for key,z in pa[pa.group.isin(['U','C4'])].groupby(['scope','regulon','TF']):
 w=z.pivot(index='patient',columns='group',values='effect').dropna()
 if {'U','C4'}<=set(w):
  rd.append(dict(zip(['scope','regulon','TF'],key),U_median=float(w.U.median()),C4_median=float(w.C4.median()),**summary(w.U-w.C4)))
add_fdr(pd.DataFrame(rd),['scope']).to_csv(D/'regulons/U_MINUS_C4_SHARED_PATIENT.tsv',sep='\t',index=False)
(D/'qc/RECEIVER_REGULON_STATUS.json').write_text(json.dumps({'status':'PASS','patient_effect_rows':len(ef),'receiver_summary_rows':len(ss),'regulon_contrast_rows':len(rs),'frozen_network_retrained':False,'cell_P_values_exported':False},indent=2))
print('RECEIVER_REGULON_DETAILS_COMPLETE',flush=True)
