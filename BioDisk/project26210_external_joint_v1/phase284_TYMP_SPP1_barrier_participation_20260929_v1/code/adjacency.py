from common import *
OLD=R/'phase234_part4_replan_execution_20260921_v1'
H=R/'phase277_metabolic_macrophage_barrier_cooperation_20260928_v1'
CORE=['source_depth','receiver_depth','source_macro','receiver_macro','source_malignant','receiver_malignant','source_SPP1','log_neighbors']
def main():
 gate('adjacency_start');reg=read(H/'tables/INPUT_SECTION_REGISTER.tsv');quant=read(PRE/'tables/SPATIAL_FROZEN_QUANTILES.tsv').set_index('sample');effects=[];qual=[];manifest=[]
 for r in reg.itertuples():
  paths=[OLD/'inputs'/f'{r.sample}_INPUTS.tsv.gz',OLD/'graphs'/f'{r.sample}_distance800.npz',H/'inputs'/f'{r.sample}_SOURCE_MEMBERSHIP.tsv.gz']
  manifest += [dict(path=str(f),sha256=sha(f)) for f in paths]
  d=read(paths[0]);D=sparse.load_npz(paths[1]).tocsr();mem=read(paths[2]);outside=d.malignant_fraction.lt(.2).to_numpy()
  band=np.zeros(len(d),bool);band[np.unique(mem.spot_index.to_numpy(int))]=True
  q=float(quant.at[r.sample,'SPP1_q75'])
  sources=np.flatnonzero(band&outside&d.source_eligible.to_numpy(bool)&d.total_macro_fraction.ge(.05).to_numpy()&d.SPP1_share.gt(q).to_numpy())
  for radius in [200,400]:
   records=[]
   for i in sources:
    lo,hi=D.indptr[i:i+2];js=D.indices[lo:hi];dist=D.data[lo:hi];keep=(dist>0)&(dist<=radius+1e-6)&outside[js]&(js!=i);nei=js[keep]
    if len(nei)<3:continue
    row=dict(spot_index=int(i),spot_id=d.spot_id.iloc[i],source_TYMP=float(d.TYMP_RNA.iloc[i]),neighbors=len(nei),neighbor_indices=';'.join(map(str,nei)),source_depth=float(d.log1p_UMI.iloc[i]),source_macro=float(d.total_macro_fraction.iloc[i]),source_malignant=float(d.malignant_fraction.iloc[i]),source_SPP1=float(d.SPP1_share.iloc[i]),source_endo=float(d.endothelial_fraction.iloc[i]),log_neighbors=float(np.log(len(nei))))
    for dest,src in [('receiver_depth','log1p_UMI'),('receiver_macro','total_macro_fraction'),('receiver_malignant','malignant_fraction'),('receiver_endo','endothelial_fraction')]:row[dest]=float(d[src].iloc[nei].mean())
    for target in ['POSTNCAF','ECM4']:row['source_'+target]=float(d[target].iloc[i]);row['receiver_'+target]=float(d[target].iloc[nei].mean())
    records.append(row)
   z=pd.DataFrame(records)
   if len(z)==0:z=pd.DataFrame(columns=['spot_index','source_TYMP']+CORE+['source_endo','receiver_endo','source_POSTNCAF','receiver_POSTNCAF','source_ECM4','receiver_ECM4'])
   z.to_csv(P/'inputs'/f'{r.sample}_ADJACENCY_{radius}.tsv.gz',sep='\t',index=False)
   qual.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,radius=radius,candidate_source_spots=len(sources),qualified_source_spots=len(z),minimum_receiver_spots=3,minimum_model_spots=50,SPP1_threshold=q,zero_self_pairs=True))
   for target in ['POSTNCAF','ECM4']:
    for model,covs in [('composition_adjusted',CORE),('local_target_vascular_sensitivity',CORE+['source_'+target,'source_endo','receiver_endo'])]:
     e,n,reason=fit(z.source_TYMP.to_numpy(float),z['receiver_'+target].to_numpy(float),z[covs].to_numpy(float),min_n=50)
     effects.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,radius=radius,target=target,model=model,effect=e,source_spots=n,reason=reason,family='distinct_spot_organization_8'))
  print('ADJACENCY_SECTION',r.sample,len(sources),flush=True)
 e=pd.DataFrame(effects);tab(e,'ADJACENCY_SECTION_EFFECTS.tsv');tab(qual,'ADJACENCY_QUALIFICATION.tsv');keys=['dataset','family','radius','target','model'];rows=[]
 for key,g in e.groupby(keys+['patient']):
  v=g.effect.dropna().to_numpy();value=float(np.tanh(np.arctanh(np.clip(v,-.99999999,.99999999)).mean())) if len(v) else np.nan
  rows.append(dict(zip(keys+['patient'],key),effect=value,n_sections=len(v),sections=';'.join(g.loc[g.effect.notna(),'sample'])))
 pa=pd.DataFrame(rows);tab(pa,'ADJACENCY_PATIENT_EFFECTS.tsv');summary(pa,keys,['dataset','family'],'ADJACENCY_COHORT_SUMMARY.tsv')
 for rec in manifest:assert sha(rec['path'])==rec['sha256']
 dump(P/'qc/ADJACENCY_INPUTS.json',manifest);dump(P/'qc/ADJACENCY_COMPLETE.json',dict(status='COMPUTED',section_effects=len(e),patient_effects=len(pa),section_radius_qualifications=len(qual)))
 print('ADJACENCY_COMPLETE',flush=True)
if __name__=='__main__':main()
