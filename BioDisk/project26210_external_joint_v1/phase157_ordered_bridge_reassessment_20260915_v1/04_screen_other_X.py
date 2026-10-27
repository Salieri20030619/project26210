from pathlib import Path
import json,subprocess,re
import numpy as np,pandas as pd
from scipy.stats import false_discovery_control
from bridge_stats import bridge,bootstrap_ab
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');R=P/'part2_rebuild_v1';O=P/'phase157_ordered_bridge_reassessment_20260915_v1';I=O/'02_inputs';D=O/'04_other_X_screen';D.mkdir(exist_ok=False)
assert json.loads((O/'03_MTA_to_new_C4/STATUS.json').read_text())['status']=='COMPLETE_MTA_PROXY_REASSESSMENT_NEW_C4'
with (D/'memory.log').open('w') as f:subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],stdout=f,stderr=subprocess.STDOUT,check=True)
npz=np.load(I/'GROUP_PSEUDOBULK.npz',allow_pickle=True);counts=npz['counts'];det=npz['detected'];genes=npz['genes'].astype(str);gi={g:i for i,g in enumerate(genes)}
group=pd.read_csv(I/'GROUP_INDEX.tsv',sep='\t');base=pd.read_csv(P/'phase154_same_library_TYMP_C4_test_20260915_v1/SAME_LIBRARY_PATIENT_METRICS.tsv',sep='\t')
mal=pd.read_csv(I/'MATCHED_MALIGNANT_LOG2CPM.tsv.gz',sep='\t',index_col=0).loc[base.patient]
mc=np.vstack([counts[group.index[(group.sample_id==s)&(group.role=='Malignant')][0]] for s in base.sample_id])
covered=set(genes[(mc>=10).sum(axis=0)>=12]);cand={};registry=[];excluded=[]
lr=pd.read_csv(P/'phase35_part1_c01_c04_cell_communication_v1/resource_lock_v1/cellchat_interactions_v1.tsv',sep='\t',low_memory=False)
secret=lr[lr.annotation=='Secreted Signaling'].copy()
for label,g in secret.groupby('ligand'):
 syms=set()
 for s in g['ligand.symbol'].dropna().astype(str):syms.update(x.strip() for x in s.split(','))
 reason=None
 if len(syms)!=1:reason='multigene_or_unmapped_ligand'
 else:
  gene=next(iter(syms))
  if gene not in gi or gene!=str(label):reason='not_direct_single_gene_ligand'
  elif gene in ['TYMP','SPP1']:reason='exposure_or_outcome_gene'
  elif gene not in covered:reason='fewer_than_12_libraries_with_10_counts'
 if reason:excluded.append(dict(label=label,reason=reason));continue
 name='ligand:'+gene;cand[name]=mal[gene].to_numpy();registry.append(dict(candidate=name,family='secreted_ligand',genes=gene,n_genes=1))
meta_genes=['TK1','SLC29A1','SLC29A2','ODC1','SRM','SMS','MTAP','CHKA','PCYT1A','PCYT2','CEPT1','AMD1','MAT2A']
for gene in meta_genes:
 if gene in covered:
  name='metabolism:'+gene;cand[name]=mal[gene].to_numpy();registry.append(dict(candidate=name,family='metabolism_transport_gene',genes=gene,n_genes=1))
 else:excluded.append(dict(label='metabolism:'+gene,reason='coverage'))
gmt=R/'06_functional_programs/run_v1/resources/MSigDB_Hallmark_2020.gmt'
for line in gmt.read_text().splitlines():
 p=line.split('\t');gs=sorted((set(p[2:])&covered)-{'TYMP','SPP1'})
 if len(gs)<10:excluded.append(dict(label='program:'+p[0],reason='fewer_than_10_covered_genes'));continue
 x=mal[gs].to_numpy();sd=x.std(axis=0,ddof=1);good=sd>0;gs=list(np.asarray(gs)[good]);x=x[:,good];sd=sd[good]
 name='program:'+p[0];cand[name]=((x-x.mean(axis=0))/sd).mean(axis=1);registry.append(dict(candidate=name,family='Hallmark_program',genes=';'.join(gs),n_genes=len(gs)))
reg=pd.DataFrame(registry);reg.to_csv(D/'CANDIDATE_REGISTRY.tsv',sep='\t',index=False);pd.DataFrame(excluded).to_csv(D/'CANDIDATE_EXCLUSIONS.tsv',sep='\t',index=False)
X=pd.DataFrame(cand,index=base.patient);X.to_csv(D/'PATIENT_CANDIDATE_VALUES.tsv',sep='\t',index_label='patient')
rows=[]
for i,r in reg.iterrows():
 name=r.candidate;m=cand[name]
 for adj in ['raw','study_adjusted']:
  v=bridge(base.TYMP_pseudobulk,m,base.C4_fraction,base.study,adj)
  v.update(bootstrap_ab(base.TYMP_pseudobulk,m,base.C4_fraction,base.study,adj))
  v.update(candidate=name,family=r.family,adjustment=adj);rows.append(v)
out=pd.DataFrame(rows)
for p in ['a_P','XY_P','b_P']:
 out[p+'_FDR_family']=np.nan
 for _,ix in out.groupby(['family','adjustment']).groups.items():
  val=out.loc[ix,p].dropna()
  if len(val):out.loc[val.index,p+'_FDR_family']=false_discovery_control(val.to_numpy())
out.to_csv(D/'ALL_CANDIDATE_BRIDGE_RESULTS.tsv',sep='\t',index=False)
# C4 receiver compatibility uses explicit public receptor.symbol subunits, not inferred labels.
c4expr=[];otherexpr=[];c4det=[]
gg=group[group.source_tissue.str.lower()=='tumor']
for patient,g in gg.groupby('patient'):
 a=g[g.role=='C4'];b=g[g.role=='OtherMacrophage'];na=a.cells.sum();nb=b.cells.sum()
 if min(na,nb)<20:continue
 ca=counts[a.index].sum(axis=0);cb=counts[b.index].sum(axis=0)
 c4expr.append(np.log2(ca/ca.sum()*1e6+1));otherexpr.append(np.log2(cb/cb.sum()*1e6+1));c4det.append(det[a.index].sum(axis=0)/na)
c4expr=np.asarray(c4expr);otherexpr=np.asarray(otherexpr);c4det=np.asarray(c4det)
receivers=[]
for name in reg.loc[reg.family=='secreted_ligand','candidate']:
 gene=name.split(':',1)[1];options=[]
 for _,row in secret[secret.ligand==gene].iterrows():
  if pd.isna(row['receptor.symbol']):continue
  subs=sorted(set(x.strip() for x in str(row['receptor.symbol']).split(',')))
  if not subs or not all(x in gi for x in subs):continue
  ids=[gi[x] for x in subs];minimum=c4det[:,ids].min(axis=1);delta=c4expr[:,ids].min(axis=1)-otherexpr[:,ids].min(axis=1)
  options.append(dict(candidate=name,receptor_label=row.receptor,receptor_genes=';'.join(subs),n_C4_patients=len(minimum),patients_all_subunits_detected_10pct=int((minimum>=.1).sum()),median_min_subunit_detection=float(np.median(minimum)),median_limiting_subunit_C4_minus_other=float(np.median(delta))))
 receivers.extend(options)
rc=pd.DataFrame(receivers);rc.to_csv(D/'ALL_RECEPTOR_COMPATIBILITY.tsv',sep='\t',index=False)
best=rc.sort_values(['patients_all_subunits_detected_10pct','median_min_subunit_detection'],ascending=False).drop_duplicates('candidate') if len(rc) else pd.DataFrame(columns=['candidate'])
raw=out[out.adjustment=='raw'].copy();adj=out[out.adjustment=='study_adjusted'].set_index('candidate')
raw['study_a_rho']=raw.candidate.map(adj.a_rho);raw['study_a_P']=raw.candidate.map(adj.a_P);raw['study_b_rho']=raw.candidate.map(adj.b_partial_rho);raw['study_b_P']=raw.candidate.map(adj.b_P);raw['study_ab']=raw.candidate.map(adj.indirect_ab);raw['study_ab_low']=raw.candidate.map(adj.ab_low);raw['study_ab_high']=raw.candidate.map(adj.ab_high)
raw['same_path_directions']=(np.sign(raw.a_rho)==np.sign(raw.study_a_rho))&(np.sign(raw.b_partial_rho)==np.sign(raw.study_b_rho))
raw['both_legs_nominal']=(raw.a_P<.05)&(raw.b_P<.05)
raw['positive_indirect']=(raw.indirect_ab>0)&(raw.study_ab>0)
raw=raw.merge(best,on='candidate',how='left')
raw['receptor_supported']=raw.patients_all_subunits_detected_10pct>=raw.n_C4_patients*.5
raw['weakest_raw_leg']=np.minimum(raw.a_rho.abs(),raw.b_partial_rho.abs())
raw['weakest_study_leg']=np.minimum(raw.study_a_rho.abs(),raw.study_b_rho.abs())
raw['support_sort']=raw.same_path_directions.astype(int)+raw.both_legs_nominal.astype(int)+raw.positive_indirect.astype(int)
raw=raw.sort_values(['support_sort','weakest_study_leg','weakest_raw_leg'],ascending=False)
# Patient omission stability for every candidate; no p-value selection.
stability=[]
for _,row in raw.iterrows():
 name=row.candidate;vals=[];aa=[];bb=[];sv=[]
 for i in range(len(base)):
  keep=np.arange(len(base))!=i
  v=bridge(base.TYMP_pseudobulk[keep],cand[name][keep],base.C4_fraction[keep],base.study[keep],'raw')
  vals.append(v['indirect_ab']);aa.append(v['a_rho']);bb.append(v['b_partial_rho'])
 for study in base.study.unique():
  keep=base.study!=study;v=bridge(base.TYMP_pseudobulk[keep],cand[name][keep],base.C4_fraction[keep],base.study[keep],'raw');sv.append(v['indirect_ab'])
 stability.append(dict(candidate=name,LOPO_ab_min=min(vals),LOPO_ab_max=max(vals),LOPO_ab_positive=sum(x>0 for x in vals),LOPO_a_sign_agreement=sum(np.sign(x)==np.sign(row.a_rho) for x in aa),LOPO_b_sign_agreement=sum(np.sign(x)==np.sign(row.b_partial_rho) for x in bb),LOSO_ab_positive=sum(x>0 for x in sv),LOSO_ab_min=min(sv),LOSO_ab_max=max(sv)))
stab=pd.DataFrame(stability);stab.to_csv(D/'ALL_CANDIDATE_OMISSION_STABILITY.tsv',sep='\t',index=False)
raw=raw.merge(stab,on='candidate');raw.to_csv(D/'CANDIDATE_EVIDENCE_RANKING.tsv',sep='\t',index=False)
status={'status':'COMPLETE_OTHER_X_DISCOVERY_SCREEN','families':reg.family.value_counts().to_dict(),'candidate_count':len(reg),'public_prior_not_proof':True,'no_FDR_gate':True,'bootstrap_intervals':'exploratory standardized-rank path product; fixed-rank study-stratified patient bootstrap 2000, not causal identification or multiplicity-corrected significance','figures_created':False,'models_retrained':False}
(D/'STATUS.json').write_text(json.dumps(status,indent=2))
cols=['candidate','a_rho','a_P','b_partial_rho','b_P','indirect_ab','study_a_rho','study_b_rho','study_a_P','study_b_P','same_path_directions','both_legs_nominal','receptor_supported','LOPO_ab_positive']
print('COUNTS',status)
for fam in raw.family.unique():print('TOP',fam,raw[raw.family==fam][cols].head(4).to_dict('records'))
