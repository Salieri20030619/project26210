from common import *
import json
from scipy import sparse
m=pd.read_csv(D/'inputs/MACRO_METADATA.tsv.gz',sep='\t',index_col=0,dtype={'candidate':str})
genes=pd.read_csv(D/'inputs/GENES.tsv',sep='\t').gene.tolist();gi={g:i for i,g in enumerate(genes)}
raw=sparse.load_npz(D/'inputs/MACRO_COUNTS.npz').tocsr();tot=m.library_counts.to_numpy()
defs={'U':['4','10','15'],'C4':['4'],'C10':['10'],'C15':['15']}
effects=[]
# Existing whole-macrophage STAT3 association: only replace the composition adjustment.
c=pd.read_csv(R/'phase179_STAT3_external_activity_validation_20260916_v1/run_v2/SC_PRIMARY_SCORE_CACHE.tsv.gz',sep='\t')
for (resource,pat),z in c.groupby(['resource','patient']):
 z=z.set_index('cell_id');md=m.loc[z.index]
 for label,cs in [('C4',['4']),('U',['4','10','15'])]:
  if len(z)<50 or (z.SPP1>0).sum()<10:continue
  effects.append(dict(module='STAT3_frozen_'+resource,group='WholeMac',patient=pat,study=md.study.iloc[0],family='STAT3_SPP1',edge='STAT3_SPP1',adjustment='depth_'+label,rho=partial(z.score,z.SPP1,[z.depth,md.candidate.isin(cs).to_numpy(int)])))
# Reuse exact Phase182 score caches and cell order for LGALS3-independent activity.
oldmeta=pd.read_csv(R/'part2_rebuild_v1/08_regulons/input_v1/CELL_METADATA.tsv.gz',sep='\t',index_col=0,dtype={'candidate':str},low_memory=False)
oldmeta=oldmeta[(~oldmeta.candidate.isin(['1','Monocyte']))&oldmeta.source_tissue.str.lower().eq('tumor')]
cache=np.load(R/'phase182_IL10_LGALS3_receiver_audit_20260916_v1/SCORE_EDGE_CACHE.npz')
K=pd.read_csv(R/'phase182_IL10_LGALS3_receiver_audit_20260916_v1/STAT3_TARGETS_WITHOUT_LGALS3.tsv',sep='\t')
assert 'IL6' not in set(K.gene)
for key in cache.files:
 resource,scoring,pat=key.split('__',2);z=oldmeta[oldmeta.curated_patient_id==pat]
 v=cache[key];assert len(v)==len(z)
 ix=m.index.get_indexer(z.index);assert (ix>=0).all()
 sp=np.log1p(raw[ix,gi['SPP1']].toarray().ravel()/tot[ix]*1e4)
 assert np.allclose(sp,v[:,2])
 il6=np.log1p(raw[ix,gi['IL6']].toarray().ravel()/tot[ix]*1e4)
 values={'STAT3':v[:,0],'LGALS3':v[:,1],'SPP1':v[:,2],'IL6':il6}
 for x,y,family in [('STAT3','LGALS3','three_edges'),('STAT3','SPP1','three_edges'),('LGALS3','SPP1','three_edges'),('IL6','STAT3','IL6_STAT3')]:
  if (values[y]>0).sum()<10 and y!='STAT3':continue
  if x in ['LGALS3','IL6'] and (values[x]>0).sum()<10:continue
  for label,cs in [('C4',['4']),('U',['4','10','15'])]:
   effects.append(dict(module='STAT3_noLGALS3_'+resource+'_'+scoring,group='WholeMac',patient=pat,study=z.study.iloc[0],family=family,edge=x+'_'+y,adjustment='depth_'+label,rho=partial(values[x],values[y],[v[:,3],z.candidate.isin(cs).to_numpy(int)])))
ef=pd.DataFrame(effects);ef.to_csv(D/'receiver/COMPOSITION_ADJUSTED_PATIENT_EFFECTS.tsv',sep='\t',index=False)
rows=[]
keys=['module','group','family','edge','adjustment']
for key,z in ef.groupby(keys):rows.append(dict(zip(keys,key),**summary(z.rho)))
add_fdr(pd.DataFrame(rows),['module','group','family','adjustment']).to_csv(D/'receiver/COMPOSITION_ADJUSTED_SUMMARY.tsv',sep='\t',index=False)
# Focused regulon activity versus independent TYMP/SPP1 RNA within each candidate pool.
auc=pd.read_csv(R/'part2_rebuild_v1/08_regulons/network_v2_20260915/AUC_ALL_CELLS.csv',index_col=0).loc[m.index]
tfcols=[c for c in auc if c.split('(')[0] in ['JUN','JUNB','CEBPD']]
es={g:np.log1p(raw[:,gi[g]].toarray().ravel()/tot*1e4) for g in ['TYMP','SPP1','IL10']}
rs=[];il10=[]
for label,cs in dict(defs,WholeMac=list(m.candidate.unique())).items():
 for pat,z in m[m.candidate.isin(cs)&m.source_tissue.str.lower().eq('tumor')].groupby('curated_patient_id'):
  if len(z)<50:continue
  ix=m.index.get_indexer(z.index);dep=np.log1p(tot[ix])
  if (es['IL10'][ix]>0).sum()>=10 and (es['SPP1'][ix]>0).sum()>=10:
   il10.append(dict(group=label,patient=pat,study=z.study.iloc[0],rho=partial(es['IL10'][ix],es['SPP1'][ix],[dep])))
  for target in ['TYMP','SPP1']:
   if (es[target][ix]>0).sum()<10:continue
   for reg in tfcols:
    for adj,cv in [('raw',[]),('depth',[dep])]:
     rs.append(dict(group=label,patient=pat,study=z.study.iloc[0],regulon=reg,target=target,adjustment=adj,rho=partial(auc.loc[z.index,reg],es[target][ix],cv)))
tfp=pd.DataFrame(rs);tfp.to_csv(D/'regulons/THREE_TF_RNA_PATIENT_EFFECTS.tsv',sep='\t',index=False)
tfs=[dict(zip(['group','regulon','target','adjustment'],key),**summary(z.rho)) for key,z in tfp.groupby(['group','regulon','target','adjustment'])]
add_fdr(pd.DataFrame(tfs),['group','target','adjustment']).to_csv(D/'regulons/THREE_TF_RNA_SUMMARY.tsv',sep='\t',index=False)
if il10:
 il=pd.DataFrame(il10);il.to_csv(D/'receiver/IL10_SPP1_PATIENTS.tsv',sep='\t',index=False)
 pd.DataFrame([dict(group=k,**summary(z.rho)) for k,z in il.groupby('group')]).to_csv(D/'receiver/IL10_SPP1_SUMMARY.tsv',sep='\t',index=False)
else:(D/'receiver/IL10_SPP1_NOT_ESTIMABLE.txt').write_text('No patient/group meets 50 cells and 10 IL10-positive / 10 SPP1-positive cells; no test generated.')
# Original fixed ligand catalogue: update only candidate-pool receptor RNA availability.
lr=pd.read_csv(R/'phase35_part1_c01_c04_cell_communication_v1/resource_lock_v1/cellchat_interactions_v1.tsv',sep='\t',low_memory=False)
registry=pd.read_csv(R/'phase157_ordered_bridge_reassessment_20260915_v1/04_other_X_screen/CANDIDATE_REGISTRY.tsv',sep='\t')
ligands=registry.loc[registry.family=='secreted_ligand','candidate'].str.replace('ligand:','',regex=False).tolist()
out=[]
for label in defs:
 z=np.load(D/'core/contrasts'/f'tumor_{label}.npz',allow_pickle=True);gn=list(z['genes']);dt=z['target_detection'];expr=z['target_log2CPM']
 for ligand in ligands:
  choices=lr[(lr.annotation=='Secreted Signaling')&(lr.ligand==ligand)]
  for receptor,row in choices.groupby('receptor'):
   symbols=set()
   for val in row['receptor.symbol'].dropna().astype(str):symbols.update(x.strip() for x in val.split(',') if x.strip())
   mapped=symbols<=set(gn) and bool(symbols)
   idx=[gn.index(g) for g in sorted(symbols)] if mapped else []
   fractions=np.min(dt[:,idx],axis=1) if mapped else np.full(len(dt),np.nan)
   out.append(dict(group=label,ligand=ligand,receptor=receptor,subunits=';'.join(sorted(symbols)),n=len(dt),all_subunits_mapped=mapped,patients_all_subunits_ge10pct=int((fractions>=.1).sum()) if mapped else 0,median_min_subunit_detection=float(np.nanmedian(fractions)) if mapped else np.nan,median_mean_subunit_log2CPM=float(np.median(expr[:,idx].mean(1))) if mapped else np.nan))
pd.DataFrame(out).to_csv(D/'receiver/FIXED_LIGAND_RECEPTOR_AVAILABILITY.tsv',sep='\t',index=False)
(D/'qc/DEPENDENCY_COMPLETION_STATUS.json').write_text(json.dumps({'status':'PASS','composition_adjusted_patient_effects':len(ef),'threeTF_RNA_patient_effects':len(tfp),'ligand_receptor_rows':len(out),'source_cache_SPP1_alignment_verified':True},indent=2))
print('DEPENDENCY_COMPLETION_DONE',flush=True)
