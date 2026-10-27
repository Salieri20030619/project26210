from pathlib import Path
import json,hashlib,csv,gzip
import pandas as pd,numpy as np
from scipy.stats import wilcoxon,false_discovery_control
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');R=P/'part2_rebuild_v1';O=P/'phase149_part2_scope_initial_comparison_20260915_v1';E=O/'evidence'
oldpath=P/'phase34_part2_y03_myeloid_pyscenic_v1/pyscenic_core_v1/auc_grnboost2_all_myeloid_v1.csv';newpath=R/'08_regulons/network_v2_20260915/AUC_ALL_CELLS.csv'
a=pd.read_csv(oldpath,index_col=0);b=pd.read_csv(newpath,index_col=0)
m=pd.read_csv(R/'07_trajectory/input_v1/CELL_METADATA.tsv.gz',sep='\t',index_col=0,dtype={'candidate':str},low_memory=False)
oldmeta=pd.read_csv(P/'phase34_part2_y03_myeloid_pyscenic_v1/input_v1/cell_metadata_v1.tsv.gz',sep='\t',index_col=0,low_memory=False)
assert m.index.is_unique and a.index.is_unique and b.index.is_unique and m.index.isin(a.index).all() and set(m.index)==set(b.index)
tf=['CEBPD(+)','JUN(+)','JUNB(+)'];assert set(tf)<=set(a.columns)&set(b.columns)
m['old_state']=oldmeta.loc[m.index,'state_v2'];rows=[]
for pat,z in m.groupby('curated_patient_id'):
 mono=z.index[z.candidate=='Monocyte'];oldsp=z.index[z.old_state=='SPP1-state Macrophage'];c4=z.index[z.candidate=='4']
 if min(len(mono),len(oldsp),len(c4))<20:continue
 for net,auc in [('old_network',a),('new_network',b)]:
  for label,cells in [('old_SPP1_state',oldsp),('new_C4',c4)]:
   effects=auc.loc[cells,tf].mean()-auc.loc[mono,tf].mean()
   for reg,v in effects.items():rows.append(dict(patient=pat,study=z.study.iloc[0],network=net,label=label,regulon=reg,cells=len(cells),monocytes=len(mono),effect=v))
p=pd.DataFrame(rows);p.to_csv(E/'MATCHED_UNIVERSE_2X2_PATIENT_EFFECTS.tsv',sep='\t',index=False)
summ=[]
for keys,g in p.groupby(['network','label','regulon']):
 v=g.effect.to_numpy();pv=wilcoxon(v,method='approx',correction=True).pvalue if len(v)>=3 and np.any(v) else 1. if len(v)>=3 else np.nan
 summ.append(dict(zip(['network','label','regulon'],keys),n_patients=g.patient.nunique(),median_effect=np.median(v),positive=int((v>0).sum()),negative=int((v<0).sum()),P=pv))
d=pd.DataFrame(summ);d['FDR_posthoc_threeTF']=np.nan
for _,idx in d.groupby(['network','label']).groups.items():d.loc[idx,'FDR_posthoc_threeTF']=false_discovery_control(d.loc[idx,'P'].to_numpy())
d.to_csv(E/'MATCHED_UNIVERSE_2X2_SUMMARY.tsv',sep='\t',index=False)
# Independent stdlib CSV check of the July overlap, avoiding pandas joins.
with gzip.open(R/'09_downstream_crosswalk/C4_CANDIDATE_CELL_IDS.tsv.gz','rt') as f:ids={row['cell_id'] for row in csv.DictReader(f,delimiter='\t')}
with gzip.open(P/'phase145_macrophage_provenance_trace_remote_v1/legacy_object/metadata.tsv.gz','rt') as f:
 oldids={row['cell_id'] for row in csv.DictReader(f,delimiter='\t') if row['macro_subtype']=='SPP1+Macro'}
assert len(ids&oldids)==36
# No absolute cross-network score subtraction or causal attribution.
meta=dict(status='PASS',posthoc=True,comparison='2 labels x 2 networks; same cells and same eligible patients; each contrast vs Monocyte',eligible_patients=p.patient.nunique(),universe_cells=len(m),no_model_rerun=True,no_figures=True,July_overlap_independent_check=36)
(E/'MATCHED_UNIVERSE_STATUS.json').write_text(json.dumps(meta,indent=2))
with (O/'SOURCE_MANIFEST.tsv').open('a') as f:
 for path in [oldpath,newpath]:
  h=hashlib.sha256()
  with path.open('rb') as g:
   for block in iter(lambda:g.read(8*1024*1024),b''):h.update(block)
  f.write(str(path)+'\t'+str(path.stat().st_size)+'\t'+h.hexdigest()+'\n')
print(d.to_string(index=False))

