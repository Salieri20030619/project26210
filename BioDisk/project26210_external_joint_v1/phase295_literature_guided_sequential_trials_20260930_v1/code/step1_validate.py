import pandas as pd,numpy as np,json,subprocess,hashlib
from pathlib import Path
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');P=R/'phase295_literature_guided_sequential_trials_20260930_v1';O=P/'step1_receiver_network'
with open(P/'logs/memory_gate.log','a') as f:subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],stdout=f,check=True)
q=pd.read_csv(O/'PATIENT_LIBRARY_QUALIFICATION.tsv',sep='\t');sel=q[q.selected]
pb=pd.read_csv(O/'RECEIVER_PSEUDOBULK.tsv.gz',sep='\t');g=pd.read_csv(O/'RECEIVER_GENE_INPUT.tsv',sep='\t');d=pd.read_csv(O/'PATIENT_CONTRASTS.tsv.gz',sep='\t',index_col=0)
checks=[]
def check(n,b):
 assert b,n
 checks.append(n)
check('independent_patients',len(sel.patient.unique())==len(sel))
check('eligibility',((sel.activated>=20)&(sel.reference>=20)).all())
for r in sel.itertuples():
 z=pb[pb.patient.eq(r.patient)]
 a=z[z.state.eq('activated')].set_index('gene').loc[d.columns]
 b=z[z.state.eq('reference')].set_index('gene').loc[d.columns]
 for label,t in [('activated',a),('reference',b)]:
  check(r.patient+label+'CPM',np.allclose(t.cpm,1e6*t.counts/t.counts.sum(),atol=1e-10))
 check(r.patient+'contrast',np.allclose(d.loc[r.patient].values,np.log2(a.cpm+1)-np.log2(b.cpm+1),atol=1e-12))
check('median_effect',np.allclose(g.median_log2_cpm_delta,d.median().reindex(g.gene)))
check('target_rule',np.array_equal(g.target,g.background&(g.median_log2_cpm_delta>=.25)&(g.positive_patient_fraction>=.7)))
check('definition_exclusion',not (g.target & g.excluded_definition_gene).any())
rank=pd.read_csv(O/'LIGAND_ACTIVITIES_ALL.tsv',sep='\t')
check('complete_ranking',np.array_equal(rank.rank_all,np.arange(1,len(rank)+1)))
check('rank_order',np.all(np.diff(rank.pearson)<=1e-12))
check('finite_scores',np.isfinite(rank.pearson).all())
check('network_validation_present',(O/'NETWORK_VALIDATION.tsv').exists())
check('R_exit',int((P/'logs/step1_network.exit').read_text())==0)
check('official_prior_checksum',all(x['verified'] for x in json.loads((P/'resources/OFFICIAL_PRIOR_WEB_CHECKSUM_VALIDATION.json').read_text())['files'].values()))
summary=dict(status='COMPUTED_AND_VALIDATED_DESCRIPTIVE',checks_passed=len(checks),checks=checks,patients=len(sel),targets=len((O/'TARGETS_FROZEN.txt').read_text().splitlines()),all_ligands=len(rank),macrophage_ligands=int(rank.macrophage_expressed.sum()),inference='model ranking only; no patient P/FDR because n<10',top_macrophage=rank[rank.macrophage_expressed].head(10).replace({np.nan:None}).to_dict('records'),named_candidates=rank[rank.test_ligand.isin(['SPP1','TGFB1','IL1A','IL1B'])].replace({np.nan:None}).to_dict('records'))
(O/'COMPLETION.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
