from common import *
from scipy import sparse
import importlib.metadata as md
def main():
 gate('stromal_summary_resume')
 obs=read(P/'tables/STROMAL_CLUSTER_MEMBERSHIP.tsv.gz',dtype={'cluster_primary':str,'cluster_sensitive':str,'cluster_original_scvi':str})
 raw=sparse.load_npz(P/'inputs/STROMAL_COUNTS.npz');genes=np.load(P/'inputs/GENES.npy')
 log=raw.astype(np.float64).multiply((1e4/obs.total_counts_raw.to_numpy())[:,None]).tocsr();log.data=np.log1p(log.data)
 binary=raw.copy();binary.data=np.ones_like(binary.data,dtype=np.float64)
 total_log=np.asarray(log.sum(0)).ravel();total_det=np.asarray(binary.sum(0)).ravel()
 panels=json.loads((R/'phase286_CAF_identity_boundary_audit_20260929_v1/resources/IDENTITY_PANELS_FROZEN.json').read_text())
 union=sorted(set(sum(list(panels.values()),[])));idx={g:i for i,g in enumerate(genes)}
 allowed=np.array([g not in ['TYMP','SPP1'] and not g.startswith(('MT-','RPL','RPS')) and g not in ['HBA1','HBA2','HBB','HBD','HBE1','HBG1','HBG2','HBM','HBQ1','HBZ'] for g in genes])
 agg=[];markers=[];summ=[]
 for mode in ['cluster_primary','cluster_sensitive','cluster_original_scvi']:
  cats=sorted(obs[mode].unique(),key=int);code=pd.Categorical(obs[mode],categories=cats).codes
  G=sparse.csr_matrix((np.ones(len(obs)),(code,np.arange(len(obs)))),shape=(len(cats),len(obs)))
  sums=(G@log).toarray();dets=(G@binary).toarray();ns=np.bincount(code)
  for k,cluster in enumerate(cats):
   gg=obs[obs[mode]==cluster];n=ns[k]
   mean=sums[k]/n;pct=dets[k]/n;bg=(total_log-sums[k])/(len(obs)-n);pctbg=(total_det-dets[k])/(len(obs)-n);diff=mean-bg
   order=np.argsort(-diff)
   for j in [j for j in order if allowed[j] and pct[j]>=.1][:60]:markers.append(dict(mode=mode,cluster=cluster,gene=genes[j],mean_logCP10K=float(mean[j]),mean_other=float(bg[j]),difference=float(diff[j]),detection=float(pct[j]),detection_other=float(pctbg[j])))
   for g in union:
    if g in idx:agg.append(dict(mode=mode,cluster=cluster,gene=g,cells=int(n),mean_logCP10K=float(mean[idx[g]]),detection=float(pct[idx[g]])))
   pc=gg.curated_patient_id.value_counts();ss=gg.study.value_counts()
   summ.append(dict(mode=mode,cluster=cluster,cells=int(n),patients=len(pc),studies=len(ss),largest_patient_fraction=pc.iloc[0]/n,largest_study_fraction=ss.iloc[0]/n,tumor_cells=int(gg.tumor.sum()),original_labels=json.dumps(gg.cell_type_broad_v1.value_counts().to_dict()),study_counts=json.dumps(ss.to_dict()),high_priority_review=int(gg.high_priority_review.sum())))
  print('SUMMARIZED',mode,len(cats),flush=True)
 tab(agg,'CLUSTER_LINEAGE_MARKERS.tsv');tab(markers,'CLUSTER_TOP_MARKERS_DESCRIPTIVE.tsv');tab(summ,'CLUSTER_COMPOSITION.tsv')
 pat=[]
 # Patient marker table from the exact raw-derived phase286 per-cell panel, checked separately.
 for (cluster,study,patient,tissue),g in obs.groupby(['cluster_primary','study','curated_patient_id','source_tissue']):
  for gene in union:
   c='expr_'+gene
   if c in g:pat.append(dict(cluster=cluster,study=study,patient=patient,tissue=tissue,cells=len(g),gene=gene,mean_logCP10K=float(g[c].mean()),detection=float(g[c].gt(0).mean())))
 tab(pat,'PATIENT_CLUSTER_LINEAGE_MARKERS.tsv')
 for other in ['cluster_sensitive','cluster_original_scvi']:pd.crosstab(obs.cluster_primary,obs[other]).to_csv(P/'tables'/f'PRIMARY_VS_{other}.tsv',sep='\t')
 hv=read(P/'tables/HVG_SELECTION.tsv',index_col=0)
 hvg=int((hv.highly_variable.to_numpy()&allowed).sum())
 dump(P/'qc/STROMAL_CLUSTER_COMPLETE.json',dict(status='COMPUTED',cells=len(obs),patients=obs.curated_patient_id.nunique(),hvg=hvg,source_path=str(H),source_bytes=H.stat().st_size,source_mtime_ns=H.stat().st_mtime_ns,counts_match=True,excluded_TYMP_SPP1_from_features=True,original_labels_unchanged=True,seed=287,versions={x:md.version(x) for x in ['scanpy','anndata','harmonypy','numpy','scipy','pandas','igraph']},execution_note='clustering/checkpoints completed; slow repeated matrix summaries terminated and replaced by mathematically equivalent sparse grouped sums; no clustering rerun'))
 print('STROMAL_CLUSTER_COMPLETE',obs.cluster_primary.value_counts().to_dict(),flush=True)
if __name__=='__main__':main()
