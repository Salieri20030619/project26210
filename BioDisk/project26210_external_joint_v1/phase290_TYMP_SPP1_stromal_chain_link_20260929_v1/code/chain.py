from common import *
from scipy import stats
PRE=R/'phase289_state_source_boundary_process_20260929_v1'
OLD=R/'phase234_part4_replan_execution_20260921_v1'
HARM=R/'phase277_metabolic_macrophage_barrier_cooperation_20260928_v1/harmonized_platform_v1'
CORE=['in_log1p_UMI','out_log1p_UMI','in_total_macro_fraction','out_total_macro_fraction','in_malignant_fraction','out_malignant_fraction','log_side_count_ratio']
VASC=['in_endothelial_fraction','out_endothelial_fraction','in_pericyte_fraction','out_pericyte_fraction']
def effect(x,y,C,mask,modifier=None):
 good=mask&np.isfinite(x)&np.isfinite(y)&np.isfinite(C).all(1)
 if modifier is not None:good &= np.isfinite(modifier)
 n=int(good.sum())
 if n<20:return np.nan,n,'fewer_than_20_edges'
 xx=stats.rankdata(x[good]);yy=stats.rankdata(y[good]);cov=stats.rankdata(C[good],axis=0)
 cov=cov[:,cov.std(0)>1e-9];cov=(cov-cov.mean(0))/cov.std(0) if cov.shape[1] else cov
 if modifier is not None:
  mm=stats.rankdata(modifier[good])
  if xx.std()<1e-9 or mm.std()<1e-9:return np.nan,n,'constant_interaction_main_effect'
  xx=(xx-xx.mean())/xx.std();mm=(mm-mm.mean())/mm.std()
  cov=np.column_stack([cov,xx,mm]);xx=xx*mm
 D=np.column_stack([np.ones(n),cov])
 if n-np.linalg.matrix_rank(D)<5:return np.nan,n,'residual_dimension_lt5'
 xres=xx-D@np.linalg.lstsq(D,xx,rcond=None)[0];yres=yy-D@np.linalg.lstsq(D,yy,rcond=None)[0]
 if xres.std()<1e-8 or yres.std()<1e-8:return np.nan,n,'constant_residual'
 return float(np.corrcoef(xres,yres)[0,1]),n,'OK'
def main():
 gate('chain_link_start')
 cfg=json.loads((HARM/'ADDENDUM_BEFORE_HARMONIZED_MODELS.json').read_text())
 c4=cfg['C4_platform15_noAPOEshared'];caf=json.loads((PRE/'PANELS_FROZEN.json').read_text())['CAF_activation']
 assert len(c4)==15 and not set(c4)&({'TYMP','SPP1'}|set(caf))
 dump(P/'resources/READOUT_DEFINITIONS.json',dict(primary='Phase234 V3_SPP1_macro_share: old SPP1 macrophage mapped share within macrophages, not current C4',sensitivity='Phase277 C4_platform15_noAPOEshared: current C4-related RNA program, not absolute abundance',C4_genes=c4,CAF_activation=caf,shared_genes=[],source_region='Phase289 frozen macro_out members',unit='local region RNA/mapping proxy'))
 q=read(PRE/'tables/BOUNDARY_QUALIFICATION.tsv');q=q[q.definition.eq('morphology_supported')];q.to_csv(P/'tables/INHERITED_BOUNDARY_QUALIFICATION.tsv',sep='\t',index=False)
 manifest=[];effects=[];qual=[]
 for r in q.itertuples():
  old=OLD/'inputs'/f'{r.sample}_INPUTS.tsv.gz';hp=HARM/'inputs'/f'{r.sample}_PROGRAMS.tsv.gz'
  d=read(old);h=read(hp)
  assert d.spot_id.is_unique and h.spot_id.is_unique
  hh=h.set_index('spot_id').reindex(d.spot_id)
  assert hh.index.tolist()==d.spot_id.tolist() and hh.C4_platform15_noAPOEshared.notna().all()
  spot=np.column_stack([d.V3_SPP1_macro_share.to_numpy(),hh.C4_platform15_noAPOEshared.to_numpy()])
  mf=PRE/'inputs'/f'{r.sample}_morphology_supported_MEMBERSHIPS.json'
  if not mf.exists():
   qual.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,radius=0,readout='none',boundary_eligible=bool(r.eligible),reason='no_boundary_memberships'));continue
  mem=json.loads(mf.read_text());lookup={(rr['radius'],rr['edge']):rr['macro_out'] for rr in mem}
  manifest.extend(dict(path=str(f),sha256=sha(f)) for f in [old,hp,mf])
  for radius in [200,400]:
   f=PRE/'inputs'/f'{r.sample}_morphology_supported_{radius}_EDGES.tsv.gz'
   if not f.exists():continue
   e=read(f);manifest.append(dict(path=str(f),sha256=sha(f)))
   aggregate=np.full((len(e),2),np.nan);counts=np.zeros((len(e),2),int)
   for i,k in enumerate(e.edge):
    ids=lookup[(radius,int(k))]
    for j in range(2):
     vals=spot[ids,j];vals=vals[np.isfinite(vals)];counts[i,j]=len(vals)
     if len(vals)>=3:aggregate[i,j]=float(vals.mean())
   for j,readout in enumerate(['SPP1_mapped_share','C4_disjoint_program']):
    m=aggregate[:,j];a=e.out_p_CAF_activation.to_numpy();mask=e.geometry_eligible.to_numpy()&bool(r.eligible)
    common=mask&np.isfinite(m)&np.isfinite(a);valid_joint=int(common.sum())>=20 and len(np.unique(m[common]))>=3 and len(np.unique(a[common]))>=3 and np.any(m[common]>0) and np.any(a[common]>0)
    joint=np.full(len(e),np.nan)
    if valid_joint:joint[common]=(stats.rankdata(m[common])/(common.sum()+1))*(stats.rankdata(a[common])/(common.sum()+1))
    out=e.copy();out['M_state']=m;out['A_activation']=a;out['J_state_activation']=joint;out['state_valid_spots']=counts[:,j];out['joint_eligible']=common&valid_joint
    out.to_csv(P/'inputs'/f'{r.sample}_{radius}_{readout}_CHAIN.tsv.gz',sep='\t',index=False)
    qual.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,radius=radius,readout=readout,boundary_eligible=bool(r.eligible),edges=len(e),joint_edges=int(common.sum()),joint_estimable=bool(valid_joint),reason='OK' if valid_joint else 'geometry_or_state_activation_variation_insufficient'))
    for model,covs in [('core',CORE),('vascular_sensitivity',CORE+VASC)]:
     plans=[('organization','tumor_TYMP_to_joint','in_g_TYMP','J_state_activation',covs,None),('organization','macro_TYMP_to_joint','macro_out_g_TYMP','J_state_activation',covs,None),('immune_distribution','joint_to_inner_CD8','J_state_activation','in_CD8_expression',covs+['out_CD8_expression','out_Tcell_fraction'],None),('immune_distribution','joint_to_inner_T','J_state_activation','in_Tcell_fraction',covs+['out_Tcell_fraction'],None),('interaction','tumor_TYMP_x_state_to_activation','in_g_TYMP','A_activation',covs,'M_state'),('interaction','macro_TYMP_x_state_to_activation','macro_out_g_TYMP','A_activation',covs,'M_state')]
     for family,endpoint,x,y,cc,z in plans:
      eff,n,reason=effect(out[x].to_numpy(),out[y].to_numpy(),out[cc].to_numpy(),mask,modifier=out[z].to_numpy() if z else None)
      effects.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,radius=radius,readout=readout,model=model,family=family,endpoint=endpoint,effect=eff,edges=n,reason=reason,estimand='rank_product_interaction_partial_rho' if z else 'partial_Spearman'))
  print('CHAIN_SECTION_DONE',r.sample,flush=True)
 tab(qual,'CHAIN_QUALIFICATION.tsv');ef=pd.DataFrame(effects);tab(ef,'CHAIN_SECTION_EFFECTS.tsv')
 keys=['dataset','radius','readout','model','family','endpoint'];rows=[]
 for key,g in ef.groupby(keys+['patient']):
  vals=g.effect.dropna().to_numpy();value=float(np.tanh(np.mean(np.arctanh(np.clip(vals,-.999999,.999999))))) if len(vals) else np.nan
  rows.append(dict(zip(keys+['patient'],key),effect=value,n_sections=len(vals),sections=';'.join(g.loc[g.effect.notna(),'sample'])))
 pt=pd.DataFrame(rows);tab(pt,'CHAIN_PATIENT_EFFECTS.tsv');summary=[]
 for key,g in pt.groupby(keys):
  vals=g.effect.dropna().to_numpy();n=len(vals);p=float(stats.wilcoxon(vals).pvalue) if n>=10 and np.any(vals!=0) else (1. if n>=10 else np.nan)
  summary.append(dict(zip(keys,key),patients=n,positive=int((vals>0).sum()),negative=int((vals<0).sum()),median_effect=float(np.median(vals)) if n else np.nan,P=p))
 su=pd.DataFrame(summary);su['FDR']=np.nan
 for _,g in su.groupby(['dataset','radius','readout','model','family']):
  vv=g[g.P.notna()]
  if len(vv):su.loc[vv.index,'FDR']=stats.false_discovery_control(vv.P.to_numpy(),method='bh')
 tab(su,'CHAIN_COHORT_SUMMARY.tsv')
 chains=[];chainpatients=[]
 chainkeys=['dataset','radius','readout','model']
 for key,g in pt.groupby(chainkeys):
  wide=g.pivot(index='patient',columns='endpoint',values='effect')
  for source,up in [('tumor','tumor_TYMP_to_joint'),('macro','macro_TYMP_to_joint')]:
   for immune,down in [('CD8','joint_to_inner_CD8'),('T','joint_to_inner_T')]:
    zz=wide[[up,down]].dropna();a=zz[up];b=zz[down]
    for patient,row in zz.iterrows():
     pattern=('positive' if row[up]>0 else 'negative' if row[up]<0 else 'zero')+'_'+('positive' if row[down]>0 else 'negative' if row[down]<0 else 'zero')
     chainpatients.append(dict(zip(chainkeys,key),source=source,immune=immune,patient=patient,TYMP_to_joint=row[up],joint_to_immune=row[down],sign_pattern=pattern,compatible_direction=bool(row[up]>0 and row[down]<0)))
    chains.append(dict(zip(chainkeys,key),source=source,immune=immune,paired_patients=len(zz),compatible=int(((a>0)&(b<0)).sum()),positive_positive=int(((a>0)&(b>0)).sum()),negative_negative=int(((a<0)&(b<0)).sum()),negative_positive=int(((a<0)&(b>0)).sum()),zero_any=int(((a==0)|(b==0)).sum()),median_upstream=float(a.median()) if len(zz) else np.nan,median_downstream=float(b.median()) if len(zz) else np.nan,P=np.nan,interpretation='same_patient_directional_concordance_not_mediation'))
 tab(chains,'SAME_PATIENT_CHAIN_CONCORDANCE.tsv');tab(chainpatients,'SAME_PATIENT_CHAIN_DETAILS.tsv')
 consistency=[]
 for key,g in pt.groupby(['dataset','radius','model','family','endpoint']):
  w=g.pivot(index='patient',columns='readout',values='effect').dropna()
  for readout in w.columns:
   vals=w[readout].to_numpy();consistency.append(dict(zip(['dataset','radius','model','family','endpoint'],key),readout=readout,common_patients=len(vals),positive=int((vals>0).sum()),negative=int((vals<0).sum()),median_effect=float(np.median(vals)) if len(vals) else np.nan))
 tab(consistency,'COMMON_PATIENT_READOUT_COMPARISON.tsv')
 dump(P/'qc/INPUT_MANIFEST.json',manifest)
 dump(P/'qc/CHAIN_ANALYSIS_COMPLETE.json',dict(status='COMPUTED',section_effects=len(ef),patient_effects=len(pt),cohort_models=len(su),chain_combinations=len(chains),no_new_boundaries=True,no_new_genes=True,source_objects_unchanged=True,causal_mediation_estimated=False,new_figures=0))
 print('CHAIN_ANALYSIS_COMPLETE',flush=True)
if __name__=='__main__':main()
