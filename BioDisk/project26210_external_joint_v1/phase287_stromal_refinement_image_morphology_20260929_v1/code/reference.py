from common import *
from scipy.optimize import nnls
def main():
 gate('reference_check')
 audit=R/'phase286_CAF_identity_boundary_audit_20260929_v1'
 a=read(audit/'inputs/CELL_IDENTITY_AUDIT.tsv.gz')
 m=read(P/'tables/STROMAL_CLUSTER_MEMBERSHIP.tsv.gz')
 cnv=read(P/'inputs/CANDIDATE_CNV_FLAGS.tsv.gz')
 assert m.cell_id.tolist()==cnv.cell_id.tolist()
 m['malignant_conservative_v1']=cnv.malignant_conservative_v1
 labels={i:'Pericyte_mural' for i in m.cluster_primary.unique()}
 labels.update({9:'Fibroblast_like_with_mural_overlap',15:'Collagen_rich_fibroblast',1:'SMC_enriched_mural',3:'SMC_enriched_mural',7:'Endothelial_mural_mixed',11:'Epithelial_stromal_mixed',5:'Mural_with_immune_signal',6:'Mural_with_immune_signal'})
 m['stromal_label_v1']=m.cluster_primary.map(labels)
 bad=m.competing_immune|m.competing_endo|m.competing_epi|m.strong_mural
 ismal=m.malignant_conservative_v1.astype(str).str.lower().isin(['true','1'])
 core=m.cluster_primary.isin([9,15])&m.cell_type_broad_v1.isin(['Fibroblast','Pericyte/SMC'])&~bad&~ismal&(m.fibro_candidate|m.cell_type_broad_v1.eq('Fibroblast'))
 m['conservative_fibro_reference']=core
 m.to_csv(P/'tables/STROMAL_ANNOTATION_V1.tsv.gz',sep='\t',index=False)
 rows=[]
 for (study,patient,tissue),g in m.groupby(['study','curated_patient_id','source_tissue']):
  rows.append(dict(study=study,patient=patient,tissue=tissue,selected_cells=len(g),fibrolike_cluster9=int(g.cluster_primary.eq(9).sum()),collagen_cluster15=int(g.cluster_primary.eq(15).sum()),conservative_fibro=int(g.conservative_fibro_reference.sum())))
 tab(rows,'PATIENT_REFINED_FIBRO_CAPACITY.tsv')
 tab([dict(cluster=k,label=v) for k,v in sorted(labels.items())],'CLUSTER_IDENTITY_ADJUDICATION.tsv')
 old=m[m.cluster_primary==9].cluster_original_scvi.value_counts()
 sens=m[m.cluster_primary==9].cluster_sensitive.value_counts()
 dump(P/'qc/ANNOTATION_FREEZE.json',dict(status='FROZEN_BEFORE_ASSOCIATIONS',membership_sha256=sha(P/'tables/STROMAL_ANNOTATION_V1.tsv.gz'),cluster9=len(m[m.cluster_primary==9]),cluster9_original_scvi_counts=old.to_dict(),cluster9_sensitive_counts=sens.to_dict(),conservative_fibro=int(core.sum()),TYMP_SPP1_used_for_labels=False))
 z=np.load(audit/'inputs/IDENTITY_PANEL_COUNTS.npz');assert z['cell_id'].tolist()==a.cell_id.tolist()
 panels=json.loads((audit/'resources/IDENTITY_PANELS_FROZEN.json').read_text())
 features=sorted(set(sum([v for k,v in panels.items() if k!='activation'],[])));ix=[z['genes'].tolist().index(g) for g in features]
 count=z['counts'][:,ix];total=z['total_counts']
 group=np.full(len(a),'',dtype=object)
 mapping={'Endothelial':'Endothelial','Epithelial cell':'Epithelial','Macrophage':'Macrophage','T cell':'Tcell'}
 for old,new in mapping.items():group[a.cell_type_broad_v1.eq(old)&a.qc_pass]=new
 mural=m.cell_type_broad_v1.eq('Pericyte/SMC')&m.stromal_label_v1.isin(['Pericyte_mural','SMC_enriched_mural'])
 lookup=pd.Index(a.cell_id);group[lookup.get_indexer(m.loc[mural,'cell_id'])]='Mural'
 group[lookup.get_indexer(m.loc[core,'cell_id'])]='Fibroblast_like'
 pb=[];vectors=[]
 for (study,patient,label),gg in a.assign(reference=group).groupby(['study','curated_patient_id','reference']):
  if not label or len(gg)<20:continue
  ii=gg.index.to_numpy();vec=count[ii].sum(0)/total[ii].sum()*1e6
  pb.append(dict(study=study,patient=patient,label=label,cells=len(ii),tumor_cells=int(gg.tumor.sum())));vectors.append(vec)
 pb=pd.DataFrame(pb);v=np.array(vectors);pb.to_csv(P/'tables/REFERENCE_PSEUDOBULK_REGISTER.tsv',sep='\t',index=False)
 pd.DataFrame(v,columns=features).to_csv(P/'inputs/REFERENCE_PSEUDOBULK_CPM.tsv.gz',sep='\t',index=False)
 classes=['Fibroblast_like','Mural','Endothelial','Epithelial','Macrophage','Tcell']
 pred=[];mixtures=[];elig=[]
 for study in sorted(pb.study.unique()):
  train=pb.study.ne(study);test=pb.study.eq(study);cent=[]
  ftrain=pb[train&pb.label.eq('Fibroblast_like')]
  allpresent=all((train&pb.label.eq(x)).any() for x in classes)
  eligible=allpresent and len(ftrain)>=2 and ftrain.study.nunique()>=2
  elig.append(dict(held_out_study=study,training_fibro_patients=len(ftrain),training_fibro_studies=ftrain.study.nunique(),eligible=eligible))
  if not eligible:continue
  for cl in classes:
   studymeans=[v[train&pb.label.eq(cl)&pb.study.eq(st)].mean(0) for st in sorted(pb.loc[train&pb.label.eq(cl),'study'].unique())]
   cent.append(np.mean(studymeans,axis=0))
  mat=np.array(cent).T;scale=np.sqrt(np.maximum(mat.mean(1),1.0));mat=mat/scale[:,None]
  def predict(y):
   coef,res=nnls(mat,y/scale)
   coef=coef/coef.sum() if coef.sum()>0 else coef
   return coef,res
  for i in np.flatnonzero(test):
   coef,res=predict(v[i]);pred.append(dict(study=study,patient=pb.patient.iloc[i],truth=pb.label.iloc[i],predicted=classes[int(coef.argmax())],fibro_coefficient=float(coef[0]),residual=res,cells=int(pb.cells.iloc[i])))
  backgrounds=[v[test&pb.label.eq(c)].mean(0) for c in classes[1:] if (test&pb.label.eq(c)).any()]
  if len(backgrounds)<5:continue
  bg=np.mean(backgrounds,axis=0)
  for i in np.flatnonzero(test&pb.label.eq('Fibroblast_like')):
   for proportion in [0.,.1,.25,.5]:
    coef,res=predict(v[i]*proportion+bg*(1-proportion))
    mixtures.append(dict(study=study,patient=pb.patient.iloc[i],true_fibro_RNA_weight=proportion,estimated_fibro_weight=float(coef[0]),absolute_error=abs(float(coef[0])-proportion),residual=res))
 tab(elig,'REFERENCE_LOSO_ELIGIBILITY.tsv');pr=pd.DataFrame(pred);mi=pd.DataFrame(mixtures)
 tab(pr,'REFERENCE_LOSO_PREDICTIONS.tsv');tab(mi,'REFERENCE_MIXTURE_RECOVERY.tsv')
 fib=pb[pb.label.eq('Fibroblast_like')]
 truth=pr.truth.eq('Fibroblast_like')
 recall=float(pr.loc[truth,'predicted'].eq('Fibroblast_like').mean()) if truth.any() else np.nan
 fp=float(pr.loc[~truth,'predicted'].eq('Fibroblast_like').mean()) if (~truth).any() else np.nan
 mae=float(mi.absolute_error.mean()) if len(mi) else np.nan
 null95=float(mi.loc[mi.true_fibro_RNA_weight==0,'estimated_fibro_weight'].quantile(.95)) if len(mi) else np.nan
 passed=len(fib)>=3 and fib.study.nunique()>=2 and recall>=.8 and fp<=.1 and mae<=.15 and null95<=.15
 dump(P/'qc/REFERENCE_TECHNICAL_GATE.json',dict(status='PASS' if passed else 'NOT_PASSED',conservative_cells=int(core.sum()),conservative_tumor_cells=int((core&m.tumor).sum()),fibro_pseudobulk_patients=len(fib),fibro_studies=fib.study.nunique(),held_out_fibro_evaluated=int(truth.sum()),pure_fibro_recall=recall,other_false_positive_rate=fp,mixture_MAE=mae,null_fibro_weight_q95=null95,biological_P_values=False,technical_RNA_weights_not_cell_fractions=True))
 print('REFERENCE_COMPLETE',len(fib),recall,fp,mae,null95,passed,flush=True)
if __name__=='__main__':main()
