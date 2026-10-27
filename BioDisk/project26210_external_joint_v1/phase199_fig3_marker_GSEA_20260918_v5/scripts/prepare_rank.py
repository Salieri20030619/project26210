from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,shutil
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=r/'phase199_fig3_marker_GSEA_20260918_v5'
s=r/'part2_rebuild_v1/04_markers_annotation/run_v4_macrophage_only/full_marker_tables/cluster_4_all_genes.tsv.gz'
a=pd.read_csv(s,sep='\t')
good=a.test_eligible.astype(str).str.lower().eq('true')&a.patient_wilcoxon_p.notna()&np.isfinite(a.median_patient_log2CPM_difference)
z=a.loc[good].copy()
assert z.gene.is_unique and z.eligible_patients.eq(26).all()
assert z.patient_wilcoxon_p.between(0,1).all()
z['ranking_score']=np.sign(z.median_patient_log2CPM_difference)*(-np.log10(z.patient_wilcoxon_p.clip(lower=np.finfo(float).tiny)))
z=z.sort_values(['ranking_score','median_patient_log2CPM_difference','gene'],ascending=[False,False,True],kind='stable')
z.to_csv(d/'tables/RANKED_PATIENT_GENES.tsv',sep='\t',index=False)
npz=np.load(r/'part2_rebuild_v1/04_markers_annotation/run_v4_macrophage_only/paired_patient_effects/cluster_4.npz',allow_pickle=True)
print('NPZ',npz.files, {k:npz[k].shape for k in npz.files})
assert npz['effects'].shape[0]==26
idxmap={str(g):i for i,g in enumerate(npz['genes'])}
effects=npz['effects']
for i,row in enumerate(z.itertuples()):
 idx=idxmap[row.gene]
 assert np.isclose(np.median(effects[:,idx]),row.median_patient_log2CPM_difference,atol=1e-6)
gmt=r/'part2_rebuild_v1/06_functional_programs/run_v1/resources/KEGG_2021_Human.gmt'
shutil.copy2(gmt,d/'inputs'/gmt.name)
(d/'qc/RANK_AUDIT.json').write_text(json.dumps({'source':str(s),'source_sha256':hashlib.sha256(s.read_bytes()).hexdigest(),'total_genes':len(a),'eligible_ranked_genes':len(z),'patients':26,'ties':int(z.ranking_score.duplicated().sum()),'zero_scores':int((z.ranking_score==0).sum()),'median_effects_verified':len(z)},indent=2))
print('RANK_READY',len(z),flush=True)
