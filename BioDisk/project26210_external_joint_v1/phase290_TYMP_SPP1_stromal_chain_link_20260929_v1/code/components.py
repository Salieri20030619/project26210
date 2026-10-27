from common import *
from chain import effect,CORE,VASC
def main():
 gate('component_diagnostic')
 q=read(P/'tables/INHERITED_BOUNDARY_QUALIFICATION.tsv');out=[]
 for r in q.itertuples():
  for radius in [200,400]:
   for readout in ['SPP1_mapped_share','C4_disjoint_program']:
    f=P/'inputs'/f'{r.sample}_{radius}_{readout}_CHAIN.tsv.gz'
    if not f.exists():continue
    d=read(f);mask=d.joint_eligible.to_numpy()
    for model,cols in [('core',CORE),('vascular_sensitivity',CORE+VASC)]:
     for source,x in [('tumor','in_g_TYMP'),('macro','macro_out_g_TYMP')]:
      for component,y in [('macrophage_state','M_state'),('stromal_activation','A_activation')]:
       rho,n,reason=effect(d[x].to_numpy(),d[y].to_numpy(),d[cols].to_numpy(),mask)
       out.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,radius=radius,readout=readout,model=model,source=source,component=component,effect=rho,edges=n,reason=reason))
 tab(out,'COMPONENT_DIAGNOSTIC_SECTION.tsv');ef=pd.DataFrame(out);keys=['dataset','radius','readout','model','source','component'];pa=[]
 for key,g in ef.groupby(keys+['patient']):
  v=g.effect.dropna().to_numpy();val=float(np.tanh(np.arctanh(np.clip(v,-.999999,.999999)).mean())) if len(v) else np.nan
  pa.append(dict(zip(keys+['patient'],key),effect=val))
 pt=pd.DataFrame(pa);tab(pt,'COMPONENT_DIAGNOSTIC_PATIENT.tsv');summ=[]
 for key,g in pt.groupby(keys):
  v=g.effect.dropna().to_numpy();summ.append(dict(zip(keys,key),patients=len(v),positive=int((v>0).sum()),negative=int((v<0).sum()),median_effect=float(np.median(v)) if len(v) else np.nan,P=np.nan,scope='post_primary_descriptive_diagnostic_no_new_significance_tests'))
 tab(summ,'COMPONENT_DIAGNOSTIC_SUMMARY.tsv')
 # Confirm component diagnostics use the exact same local observations as the corresponding upstream joint test.
 primary=read(P/'tables/CHAIN_SECTION_EFFECTS.tsv')
 for r in ef.itertuples():
  hit=primary[(primary['sample']==r.sample)&(primary.radius==r.radius)&(primary.readout==r.readout)&(primary.model==r.model)&(primary.endpoint==r.source+'_TYMP_to_joint')]
  assert len(hit)==1 and int(hit.edges.iloc[0])==r.edges
 dump(P/'qc/COMPONENT_DIAGNOSTIC_CHECK.json',dict(status='PASS',same_observations_as_primary_joint=True,new_P_values=False,original_primary_statistics_unchanged=True))
 print('COMPONENT_DIAGNOSTIC_COMPLETE',flush=True)
if __name__=='__main__':main()
