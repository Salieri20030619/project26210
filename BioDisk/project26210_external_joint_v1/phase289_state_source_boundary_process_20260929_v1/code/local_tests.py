from common import *
from scipy import stats
def partial(x,y,covs,mask):
 good=mask&np.isfinite(x)&np.isfinite(y)&np.isfinite(covs).all(1)
 n=int(good.sum())
 if n<20:return np.nan,n,'fewer_than_20_edges'
 xx=stats.rankdata(x[good]);yy=stats.rankdata(y[good]);cc=stats.rankdata(covs[good],axis=0)
 cc=cc[:,cc.std(0)>1e-9];cc=(cc-cc.mean(0))/cc.std(0) if cc.shape[1] else cc
 D=np.column_stack([np.ones(n),cc]);rank=np.linalg.matrix_rank(D)
 if n-rank<5:return np.nan,n,'residual_df_lt5'
 xx=xx-D@np.linalg.lstsq(D,xx,rcond=None)[0];yy=yy-D@np.linalg.lstsq(D,yy,rcond=None)[0]
 if np.std(xx)<1e-8 or np.std(yy)<1e-8:return np.nan,n,'constant_residual'
 return float(np.corrcoef(xx,yy)[0,1]),n,'OK'
def main():
 gate('boundary_associations')
 q=read(P/'tables/BOUNDARY_QUALIFICATION.tsv');effects=[]
 core=['in_log1p_UMI','out_log1p_UMI','in_total_macro_fraction','out_total_macro_fraction','in_malignant_fraction','out_malignant_fraction','log_side_count_ratio']
 extra=['in_endothelial_fraction','out_endothelial_fraction','in_pericyte_fraction','out_pericyte_fraction']
 for r in q.itertuples():
  if r.definition=='none':continue
  for radius in [200,400]:
   f=P/'inputs'/f'{r.sample}_{r.definition}_{radius}_EDGES.tsv.gz'
   if not f.exists():continue
   d=read(f);mask=d.geometry_eligible.to_numpy()&bool(r.eligible)
   for model,cols in [('core',core),('vascular_sensitivity',core+extra)]:
    for endpoint,x,y,family in [('tumor_side_TYMP_CAF','in_g_TYMP','out_p_CAF_activation','organization'),('macro_region_TYMP_CAF','macro_out_g_TYMP','out_p_CAF_activation','organization'),('macro_region_SPP1_CAF','macro_out_g_SPP1','out_p_CAF_activation','organization'),('joint_CAF_SPP1_inner_CD8','joint_rank','in_CD8_expression','immune_distribution')]:
     covs=cols+(['out_CD8_expression','out_Tcell_fraction'] if family=='immune_distribution' else [])
     e,n,reason=partial(d[x].to_numpy(),d[y].to_numpy(),d[covs].to_numpy(),mask)
     effects.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,definition=r.definition,radius=radius,model=model,family=family,endpoint=endpoint,effect=e,edges=n,reason=reason))
 tab(effects,'LOCAL_SECTION_EFFECTS.tsv');e=pd.DataFrame(effects)
 keys=['dataset','definition','radius','model','family','endpoint'];patients=[]
 for key,g in e.groupby(keys+['patient']):
  v=g.effect.dropna().to_numpy();eff=float(np.tanh(np.mean(np.arctanh(np.clip(v,-.999999,.999999))))) if len(v) else np.nan
  patients.append(dict(zip(keys+['patient'],key),effect=eff,sections=';'.join(g.loc[g.effect.notna(),'sample']),n_sections=len(v)))
 pt=pd.DataFrame(patients);tab(pt,'LOCAL_PATIENT_EFFECTS.tsv')
 out=[]
 for key,g in pt.groupby(keys):
  v=g.effect.dropna().to_numpy();n=len(v);p=float(stats.wilcoxon(v).pvalue) if n>=10 and np.any(v!=0) else (1. if n>=10 else np.nan)
  out.append(dict(zip(keys,key),patients=n,positive=int((v>0).sum()),negative=int((v<0).sum()),median_effect=float(np.median(v)) if n else np.nan,P=p))
 su=pd.DataFrame(out);su['FDR']=np.nan
 for _,g in su.groupby(['dataset','definition','radius','model','family']):
  z=g[g.P.notna()]
  if len(z):su.loc[z.index,'FDR']=stats.false_discovery_control(z.P.to_numpy(),method='bh')
 tab(su,'LOCAL_COHORT_SUMMARY.tsv')
 loo=[]
 for key,g in pt.groupby(keys):
  g=g[g.effect.notna()]
  if len(g)<10:continue
  for pat in g.patient:
   v=g.loc[g.patient!=pat,'effect'];loo.append(dict(zip(keys,key),omitted_patient=pat,patients=len(v),positive=int((v>0).sum()),median_effect=float(v.median())))
 tab(loo,'LOCAL_LEAVE_PATIENT.tsv')
 dump(P/'qc/LOCAL_ASSOCIATIONS_COMPLETE.json',dict(status='COMPUTED',section_effects=len(e),patient_effects=len(pt),summary_rows=len(su),new_CNV_inference=False,source_labels='compartment_enrichment_not_single_cell_attribution',new_figures=0))
 print('LOCAL_ASSOCIATIONS_COMPLETE',flush=True)
if __name__=='__main__':main()
