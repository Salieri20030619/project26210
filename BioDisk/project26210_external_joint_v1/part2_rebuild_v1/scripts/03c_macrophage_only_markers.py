from pathlib import Path
import json,datetime,math,warnings,gc
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from scipy.stats import wilcoxon
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/part2_rebuild_v1');C=R/(R/'provenance/CURRENT_CLUSTERING_RUN.txt').read_text().strip();O=R/'04_markers_annotation/run_v3_macrophage_only';V=R/'05_patient_study_validation/run_v3_macrophage_only'
assert json.loads((C/'STATUS.json').read_text())['status']in ['PASS_UNSUPERVISED_CLUSTERING_REVIEW','APPROVED_EXPLORATORY_MARKER_IDENTIFICATION']
O.mkdir();V.mkdir();(O/'full_marker_tables').mkdir();(O/'paired_patient_effects').mkdir()
def log(s):print(datetime.datetime.now().isoformat(),s,flush=True)
def bh(p):
 p=np.asarray(p,float);out=np.full(p.shape,np.nan);valid=np.isfinite(p);a=p[valid];order=np.argsort(a);vals=np.minimum.accumulate((a[order]*len(a)/np.arange(1,len(a)+1))[::-1])[::-1];rev=np.empty(len(a));rev[order]=np.minimum(vals,1);out[valid]=rev;return out
log('load raw expression by exact selected cell identity')
z=ad.read_h5ad(C/'CLUSTERING_GRAPH.h5ad');z=z[z.obs.cluster.astype(str).isin(['0','3','4','5','6','7','8','9','10','11','12','15'])].copy();assert z.n_obs==22444;allmy=ad.read_h5ad(R/'01_input_audit/run_v1/myeloid_raw_counts_v1.h5ad');a=allmy[z.obs_names].copy();mono=allmy[allmy.obs.cell_type_broad_v1.astype(str)=='Monocyte'].copy();del allmy;gc.collect()
x=a.X.tocsr();genes=a.var_names.astype(str);obs=z.obs.copy();obs['cluster']=obs.cluster.astype(str);patients=sorted(obs.curated_patient_id.astype(str).unique());clusters=sorted(obs.cluster.unique(),key=int);pstudy=obs.groupby('curated_patient_id',observed=True).study.first().astype(str).to_dict();assert len(patients)==40
pi={p:i for i,p in enumerate(patients)};ci={c:i for i,c in enumerate(clusters)}
row=np.array([pi[str(p)]*len(clusters)+ci[c] for p,c in zip(obs.curated_patient_id,obs.cluster)]);mat=sparse.csr_matrix((np.ones(len(obs)),(row,np.arange(len(obs)))),shape=(len(patients)*len(clusters),len(obs)));pb=(mat@x).tocsr();nc=np.asarray(mat.sum(axis=1)).ravel().astype(int);pbt=pb.toarray().reshape(len(patients),len(clusters),len(genes));pt=pbt.sum(axis=1);nct=nc.reshape(len(patients),len(clusters));sparse.save_npz(O/'PATIENT_CLUSTER_PSEUDOBULK_COUNTS.npz',pb);pd.DataFrame([(p,c,pstudy[p],int(nct[pi[p],ci[c]])) for p in patients for c in clusters],columns=['patient','cluster','study','n_cells']).to_csv(O/'PSEUDOBULK_ROW_INDEX.tsv',sep='\t',index=False);pd.Series(genes,name='gene').to_csv(O/'PSEUDOBULK_GENE_INDEX.tsv',sep='\t',index=False)
norm=x.astype(np.float64).multiply(10000/np.asarray(x.sum(axis=1))).tocsr();norm.data=np.log1p(norm.data);global_sum=np.asarray(norm.sum(axis=0)).ravel();global_npos=np.asarray((x>0).sum(axis=0)).ravel();global_raw=np.asarray(x.sum(axis=0)).ravel();global_total=global_raw.sum();metrics=[];top=[];canonical=[];study_effects=[]
canonical_genes=['SPP1','GPNMB','FABP5','CD9','LPL','LGALS3','C1QA','C1QB','C1QC','APOE','APOC1','TREM2','FOLR2','LYVE1','MRC1','CD163','IL1B','CXCL8','CCL3','CCL4','FCN1','S100A8','S100A9','VCAN','HLA-DRA','CD74','HLA-DPB1','ISG15','IFIT1','IFIT3','MX1','MKI67','TOP2A','STMN1','LST1','TYROBP','CSF1R','FCER1G','CD68','CTSB','EPCAM','KRT8','KRT18','CD3D','CD3E','TRAC','NKG7','GNLY','KLRD1','MS4A1','CD79A','PECAM1','VWF','TPSAB1']
for c in clusters:
 log('patient-paired complete marker family cluster '+c);k=ci[c];mask=obs.cluster.to_numpy()==c;n=int(mask.sum());inside=np.asarray(norm[mask].sum(axis=0)).ravel();npos=np.asarray((x[mask]>0).sum(axis=0)).ravel();raw=np.asarray(x[mask].sum(axis=0)).ravel();pooled_cpm=raw/raw.sum()*1e6;rest_cpm=(global_raw-raw)/(global_total-raw.sum())*1e6
 eligible=(nct[:,k]>=20)&((nct.sum(axis=1)-nct[:,k])>=20);ii=np.flatnonzero(eligible);npat=len(ii);bx=pbt[ii,k,:];br=pt[ii]-bx
 if npat:
  logc=np.log2(bx/bx.sum(axis=1)[:,None]*1e6+1);logr=np.log2(br/br.sum(axis=1)[:,None]*1e6+1);diff=logc-logr;med=np.median(diff,axis=0);pos=(diff>0).sum(axis=0);neg=(diff<0).sum(axis=0)
 else:diff=np.zeros((0,len(genes)));med=np.full(len(genes),np.nan);pos=neg=np.zeros(len(genes),int)
 test=(np.sum((bx>=10)|(br>=10),axis=0)>=max(3,math.ceil(.2*npat)))&np.any(diff!=0,axis=0) if npat>=3 else np.zeros(len(genes),bool)
 pv=np.full(len(genes),np.nan)
 if test.any():
  with warnings.catch_warnings():
   warnings.simplefilter('ignore');pv[test]=wilcoxon(diff[:,test],axis=0,zero_method='wilcox',correction=True,alternative='two-sided',method='approx').pvalue
 q=bh(pv);pctin=npos/n;pctout=(global_npos-npos)/(len(obs)-n);pooledlf=np.log2(pooled_cpm+1)-np.log2(rest_cpm+1);priority=np.nan_to_num(med,nan=0)*np.maximum(pctin-pctout,0)
 tab=pd.DataFrame({'cluster':c,'gene':genes,'n_cells':n,'pct_in':pctin,'pct_out':pctout,'mean_log1pCP10K_in':inside/n,'mean_log1pCP10K_out':(global_sum-inside)/(len(obs)-n),'log2CPM_difference_pooled':pooledlf,'eligible_patients':npat,'median_patient_log2CPM_difference':med,'n_patients_positive':pos,'n_patients_negative':neg,'patient_wilcoxon_p':pv,'BH_FDR_cluster_gene_family':q,'test_eligible':test,'reporting_tier':'exploratory_n_lt10' if npat<10 else 'patient_FDR_conditional_marker_discovery','annotation_priority':priority})
 tab.to_csv(O/'full_marker_tables'/('cluster_'+c+'_all_genes.tsv.gz'),sep='\t',index=False,compression='gzip')
 np.savez_compressed(O/'paired_patient_effects'/('cluster_'+c+'.npz'),effects=diff,patient=np.array(patients)[ii],study=np.array([pstudy[patients[i]] for i in ii]),genes=np.array(genes),tested=test)
 ranked=tab[(tab.pct_in>=.1)&(tab.pct_in>tab.pct_out)&(tab.median_patient_log2CPM_difference>0)].sort_values(['annotation_priority','median_patient_log2CPM_difference'],ascending=False).head(40);top.append(ranked)
 can=tab[tab.gene.isin(canonical_genes)];canonical.append(can)
 usedgenes=list(dict.fromkeys(ranked.gene.tolist()+[g for g in canonical_genes if g in genes]));gi=genes.get_indexer(usedgenes)
 for st in sorted(set(pstudy.values())):
  kk=np.array([pstudy[patients[i]]==st for i in ii]);vals=diff[kk][:,gi]
  for j,g in enumerate(usedgenes):study_effects.append({'cluster':c,'study':st,'gene':g,'eligible_patients':int(kk.sum()),'median_patient_effect':float(np.median(vals[:,j])) if kk.any() else np.nan,'positive_patients':int((vals[:,j]>0).sum())})
 cellstudy=obs.loc[mask,'study'].astype(str).value_counts();cellpatient=obs.loc[mask,'curated_patient_id'].astype(str).value_counts();qualified=[]
 for st in sorted(obs.study.astype(str).unique()):
  g=obs.loc[mask&(obs.study.astype(str).to_numpy()==st)]
  if len(g)>=20 and g.curated_patient_id.nunique()>=2:qualified.append(st)
 dom=float(cellstudy.iloc[0]/n);robust=len(qualified)>=3 and len(cellpatient)>=5 and dom<=.8
 metrics.append({'cluster':c,'n_cells':n,'patients_any_cells':len(cellpatient),'patients_paired_marker_test':npat,'studies_any_cells':len(cellstudy),'qualified_studies':';'.join(qualified),'n_qualified_studies':len(qualified),'dominant_study':cellstudy.index[0],'dominant_study_fraction':dom,'dominant_patient_fraction':float(cellpatient.iloc[0]/n),'lineage_review_flag_cells':int(obs.loc[mask,'lineage_review_flag'].sum()),'lineage_review_flag_fraction':float(obs.loc[mask,'lineage_review_flag'].mean()),'median_total_counts':float(obs.loc[mask,'total_counts_raw'].median()),'median_n_genes':float(obs.loc[mask,'n_genes_raw'].median()),'median_pct_mito':float(obs.loc[mask,'pct_counts_mt_raw'].median()),'tested_genes':int(test.sum()),'positive_patient_FDR_genes':int(((q<.05)&(med>0)).sum()),'broad_study_distribution_gate':robust})
pd.crosstab(obs.cluster,obs.source_tissue).to_csv(V/'CLUSTER_BY_TISSUE.tsv',sep='\t')
pd.concat(top).to_csv(O/'TOP40_POSITIVE_PATIENT_MARKERS.tsv',sep='\t',index=False);pd.concat(canonical).to_csv(O/'CANONICAL_PANEL_ALL_CLUSTERS.tsv',sep='\t',index=False);pd.DataFrame(study_effects).to_csv(V/'TOP_AND_CANONICAL_MARKER_STUDY_EFFECTS.tsv',sep='\t',index=False);pd.DataFrame(metrics).to_csv(V/'CLUSTER_COMPOSITION_QC.tsv',sep='\t',index=False)
# Monocyte reference is descriptive and outside clustering selection.
use=[g for g in canonical_genes if g in mono.var_names];ii=mono.var_names.get_indexer(use);mx=mono.X.tocsr();total=np.asarray(mx.sum(axis=1)).ravel();sm=mx[:,ii].astype(float).multiply(10000/total[:,None]).tocsr();sm.data=np.log1p(sm.data);pd.DataFrame({'gene':use,'monocyte_mean_log1pCP10K':np.asarray(sm.mean(axis=0)).ravel(),'monocyte_pct_positive':np.asarray((mx[:,ii]>0).mean(axis=0)).ravel()}).to_csv(O/'MONOCYTE_CONTEXT_MARKERS.tsv',sep='\t',index=False)
status={'status':'MARKER_AND_COMPOSITION_TABLES_COMPLETE_PENDING_ANNOTATION','clusters':len(clusters),'genes_per_cluster':len(genes),'patients':len(patients),'cell_level_P_values_generated':False,'marker_tests_are_conditional_not_independent_validation':True};(O/'STATUS.json').write_text(json.dumps(status,indent=2));print(json.dumps(status),flush=True)
print(pd.DataFrame(metrics).to_string(index=False),flush=True)
for c in clusters:print('TOP',c,pd.concat(top).query('cluster==@c').head(15)[['gene','median_patient_log2CPM_difference','BH_FDR_cluster_gene_family']].to_dict('records'),flush=True)
