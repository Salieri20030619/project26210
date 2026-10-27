from pathlib import Path
import pandas as pd,numpy as np,json
from scipy.stats import spearmanr
w=Path(__file__).resolve().parents[1];r=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/05_cnmf')
c=pd.read_csv(w/'data/rank7_run_component_loadings.tsv.gz',sep='\t').set_index('gene');g=pd.read_csv(w/'data/rank7_component_labels.tsv',sep='\t').set_index('component_id').loc[c.columns];target=pd.read_csv(r/'full_ranks4_8_v1/consensus_loading_rank7_v1.tsv.gz',sep='\t').set_index('gene');reconstructed={}
for k in range(1,8):
 selected=g[g.consensus_cluster==k];byrun=[]
 for run,rows in selected.groupby('run'):byrun.append(c[rows.index].mean(axis=1).to_numpy())
 vals=np.mean(byrun,axis=0);reconstructed['F'+str(k)]=vals/vals.sum()
reconstructed=pd.DataFrame(reconstructed,index=c.index);err=float(np.max(np.abs(reconstructed.loc[target.index].values-target.values)));assert err<1e-12,err
u=pd.read_csv(r/'full_ranks4_8_v1/consensus_usage_rank7_v1.tsv.gz',sep='\t').set_index('exact_cell_id');m=pd.read_csv(r/'input_v1/cnmf_cell_metadata_v1.tsv.gz',sep='\t').set_index('exact_cell_id').loc[u.index]
pc=pd.read_csv(r/'tymp_rank_evidence_v1/factor_patient_stratified_correlations_v1.tsv.gz',sep='\t');pc=pc[(pc['rank']==7)&(pc.factor=='F1')&(pc.outcome=='TYMP_log1p_CP10K_v1')]
rechecks=[]
for row in pc.itertuples():
 ii=m.curated_patient_id==row.curated_patient_id;rho=float(spearmanr(u.loc[ii,'F1'],m.loc[ii,'TYMP_log1p_CP10K_v1']).statistic);rechecks.append({'patient':row.curated_patient_id,'source_rho':row.rho,'recomputed_rho':rho,'absolute_delta':abs(rho-row.rho)})
df=pd.DataFrame(rechecks);df.to_csv(w/'qc/INDEPENDENT_PATIENT_RECHECK.tsv',sep='\t',index=False)
maxerr=float(df.absolute_delta.max());assert maxerr<1e-5,maxerr
result={'status':'PASS','rank7_consensus_loading_reconstruction_max_abs_error':err,'components':140,'original_runs':20,'patient_rho_rechecks':len(df),'patient_rho_max_abs_error':maxerr,'positive_patients':int((df.recomputed_rho>0).sum()),'inference_recalculated':False}
(w/'qc/INDEPENDENT_VALIDATION.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))
