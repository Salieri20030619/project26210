from source import *
H=R/'phase277_metabolic_macrophage_barrier_cooperation_20260928_v1'
HP=H/'harmonized_platform_v1'
CORE=['inside_depth','outside_depth','inside_macro','outside_macro','inside_malignant','outside_malignant','log_side_count_ratio_new']
OUT=['outside_C4','outside_POSTN','outside_ECM','SPP1_POSTN_coverage','C1QC_POSTN_coverage','CD8_RNA_delta','T_abundance_log2_delta']
def interaction(z,y,cov):
 cols=['inside_TYMP','SPP1_POSTN_coverage',y]+cov
 a=z.loc[z.split_eligible,cols].dropna().to_numpy();n=len(a)
 o=dict(effect=np.nan,n=n,reason='',condition=np.nan,beta_TYMP=np.nan,beta_coverage=np.nan,coverage_slope_TYMP_low=np.nan,coverage_slope_TYMP_high=np.nan,quadrant_min=0)
 if n<50:o['reason']='FEWER_THAN_50_SEGMENTS';return o
 sd=a.std(0)
 if np.any(sd==0):o['reason']='CONSTANT_VARIABLE';return o
 a=(a-a.mean(0))/sd;x,f,yy=a[:,0],a[:,1],a[:,2];q=[int((((x>=np.median(x))==i)&((f>=np.median(f))==j)).sum()) for i in [False,True] for j in [False,True]];o['quadrant_min']=min(q)
 if min(q)<5:o['reason']='QUADRANT_LT5';return o
 D=np.column_stack([np.ones(n),x,f,x*f,a[:,3:]])
 o['condition']=float(np.linalg.cond(D))
 if np.linalg.matrix_rank(D)<D.shape[1] or o['condition']>=30:o['reason']='RANK_OR_CONDITION';return o
 beta=np.linalg.lstsq(D,yy,rcond=None)[0];o.update(effect=float(beta[3]),beta_TYMP=float(beta[1]),beta_coverage=float(beta[2]),coverage_slope_TYMP_low=float(beta[2]-beta[3]),coverage_slope_TYMP_high=float(beta[2]+beta[3]))
 return o
def main():
 gate('spatial_start');reg=read(H/'tables/INPUT_SECTION_REGISTER.tsv');ef=[];qual=[];hashes=[];details=[];thresholds=[]
 for i,r in enumerate(reg.itertuples()):
  paths=[OLD/'inputs'/f'{r.sample}_INPUTS.tsv.gz',HP/'inputs'/f'{r.sample}_PROGRAMS.tsv.gz',HP/'inputs'/f'{r.sample}_MODEL_INPUT.tsv.gz',H/'inputs'/f'{r.sample}_SOURCE_MEMBERSHIP.tsv.gz']
  hashes.extend(dict(path=str(f),sha256=sha(f)) for f in paths)
  d=read(paths[0]);sc=read(paths[1]);z=read(paths[2]);mem=read(paths[3])
  assert np.array_equal(d.spot_id,sc.spot_id) and np.array_equal(z.segment,np.arange(len(z)))
  inside=d.malignant_fraction.ge(.2).to_numpy()
  qsp=float(d.loc[~inside,'SPP1_share'].quantile(.75));qc1=float(d.loc[~inside,'C1QC_fraction'].quantile(.75));qf=float(d.loc[~inside,'POSTNCAF'].quantile(.75))
  d['jointSPP1']=(d.SPP1_share>qsp)&(d.POSTNCAF>qf);d['jointC1QC']=(d.C1QC_fraction>qc1)&(d.POSTNCAF>qf)
  thresholds.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,outside_spots=int((~inside).sum()),SPP1_q75=qsp,C1QC_q75=qc1,POSTN_q75=qf,coverage_strict_greater=True))
  cols=['inside_TYMP','outside_TYMP','inside_n','outside_n']+CORE+['inside_endothelial','outside_endothelial']+OUT[:5]
  for col in cols:z[col]=np.nan
  for sid,g in mem.groupby('segment'):
   ids=g.spot_index.to_numpy(int);a=ids[inside[ids]];b=ids[~inside[ids]];z.at[sid,'inside_n']=len(a);z.at[sid,'outside_n']=len(b)
   if not len(a) or not len(b):continue
   for side,jj in [('inside',a),('outside',b)]:
    for dest,src in [('TYMP','TYMP_RNA'),('depth','log1p_UMI'),('macro','total_macro_fraction'),('malignant','malignant_fraction'),('endothelial','endothelial_fraction')]:z.at[sid,side+'_'+dest]=float(d[src].iloc[jj].mean())
   z.at[sid,'log_side_count_ratio_new']=float(np.log((len(a)+.5)/(len(b)+.5)))
   z.at[sid,'outside_C4']=float(sc.C4_platform17.iloc[b].mean())
   for dest,src in [('outside_POSTN','POSTNCAF'),('outside_ECM','ECM4'),('SPP1_POSTN_coverage','jointSPP1'),('C1QC_POSTN_coverage','jointC1QC')]:z.at[sid,dest]=float(d[src].iloc[b].mean())
  z['split_eligible']=z.source_eligible&(z.inside_n>=5)&(z.outside_n>=5)
  z.to_csv(P/'inputs'/f'{r.sample}_SPLIT_INTERFACE.tsv.gz',sep='\t',index=False)
  configs=[('unadjusted',[]),('core_adjusted',CORE),('outside_TYMP_adjusted',CORE+['outside_TYMP']),('vascular_sensitivity',CORE+['outside_TYMP','inside_endothelial','outside_endothelial'])]
  for y in OUT:
   for model,cov in configs:
    e,n,reason=fit(z.inside_TYMP,z[y],z[cov].to_numpy() if cov else None,z.split_eligible.to_numpy(bool),min_n=50)
    ef.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,family='spatial_cross_compartment',endpoint=y,model=model,effect_type='rho',effect=e,n_segments=n,reason=reason))
  for y in OUT[-2:]:
   for model,cov in configs[1:3]:
    o=interaction(z,y,cov);ef.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,family='spatial_TYMP_coverage_interaction',endpoint=y,model=model,effect_type='interaction_beta',effect=o['effect'],n_segments=o['n'],reason=o['reason']))
    details.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,endpoint=y,model=model,**o))
  qual.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,segments=len(z),old_source_eligible=int(z.source_eligible.sum()),split_eligible=int(z.split_eligible.sum()),min_each_side=5,nonoverlap_asserted=True))
  print('SPATIAL_SECTION',i+1,r.sample,int(z.split_eligible.sum()),flush=True)
 ef=pd.DataFrame(ef);tab(ef,'SPATIAL_SECTION_EFFECTS.tsv');tab(qual,'SPATIAL_QUALIFICATION.tsv');tab(details,'SPATIAL_INTERACTION_DETAILS.tsv');tab(thresholds,'SPATIAL_FROZEN_QUANTILES.tsv')
 keys=['dataset','family','endpoint','model','effect_type'];pa=[]
 for key,g in ef.groupby(keys+['patient']):
  v=g.effect.dropna().to_numpy();e=(float(np.tanh(np.arctanh(np.clip(v,-.99999999,.99999999)).mean())) if key[4]=='rho' else float(v.mean())) if len(v) else np.nan
  pa.append(dict(zip(keys+['patient'],key),effect=e,sections=';'.join(g.loc[g.effect.notna(),'sample']),n_sections=len(v)))
 pa=pd.DataFrame(pa);tab(pa,'SPATIAL_PATIENT_EFFECTS.tsv');su=[]
 for key,g in pa.groupby(keys):su.append(dict(zip(keys,key),**infer(g.effect,2624280+len(su))))
 su=bh(pd.DataFrame(su),['dataset','family']);tab(su,'SPATIAL_COHORT_SUMMARY.tsv');loo=[]
 for key,g in pa.groupby(keys):
  g=g[g.effect.notna()]
  if len(g)<10:continue
  for r in g.itertuples():
   v=g.loc[g.patient!=r.patient,'effect'].to_numpy();loo.append(dict(zip(keys,key),omitted_patient=r.patient,n=len(v),median=float(np.median(v)),positive=int((v>0).sum()),negative=int((v<0).sum())))
 tab(loo,'SPATIAL_LEAVE_PATIENT.tsv')
 for rec in hashes:assert sha(rec['path'])==rec['sha256']
 dump(P/'qc/SPATIAL_INPUTS.json',hashes);dump(P/'qc/SPATIAL_COMPLETED.json',dict(at=now(),sections=len(reg),section_effects=len(ef),cohort_endpoints=len(su),figures=0))
 print('SPATIAL_COMPLETE',flush=True)
if __name__=='__main__':main()
