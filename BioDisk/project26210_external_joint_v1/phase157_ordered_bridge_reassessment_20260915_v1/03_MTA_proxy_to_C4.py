from pathlib import Path
import json,subprocess
import numpy as np,pandas as pd
from scipy.stats import wilcoxon,false_discovery_control
from bridge_stats import bridge,bootstrap_ab
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');O=P/'phase157_ordered_bridge_reassessment_20260915_v1';I=O/'02_inputs';D=O/'03_MTA_to_new_C4';D.mkdir(exist_ok=False)
assert json.loads((I/'STATUS.json').read_text())['status']=='INPUTS_COMPLETE'
with (D/'memory.log').open('w') as f:subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],stdout=f,stderr=subprocess.STDOUT,check=True)
z=np.load(I/'GROUP_PSEUDOBULK.npz',allow_pickle=True);counts=z['counts'];genes=z['genes'].astype(str);gi={g:i for i,g in enumerate(genes)}
groups=pd.read_csv(I/'GROUP_INDEX.tsv',sep='\t');cell=pd.read_csv(I/'CELL_MTA_RESPONSE_SCORES.tsv.gz',sep='\t')
base=pd.read_csv(P/'phase154_same_library_TYMP_C4_test_20260915_v1/SAME_LIBRARY_PATIENT_METRICS.tsv',sep='\t')
expr=[];detect=[];metrics=base.copy()
for _,row in base.iterrows():
 ix=groups.index[(groups.sample_id==row.sample_id)&(groups.role=='Malignant')];assert len(ix)==1
 x=counts[ix[0]];expr.append(np.log2(x/x.sum()*1e6+1));detect.append(z['detected'][ix[0]]/groups.loc[ix[0],'cells'])
mal=pd.DataFrame(expr,index=base.patient,columns=genes);mal.to_csv(I/'MATCHED_MALIGNANT_LOG2CPM.tsv.gz',sep='\t',index_label='patient')
pd.DataFrame(detect,index=base.patient,columns=genes).to_csv(I/'MATCHED_MALIGNANT_DETECTION.tsv.gz',sep='\t',index_label='patient')
assert np.allclose(mal.TYMP.to_numpy(),base.TYMP_pseudobulk)
scoremap=cell[cell.role=='Malignant'].groupby('sample_id').MTA_response_transcript_score.mean()
metrics['direct_MTA_response_RNA_proxy']=metrics.sample_id.map(scoremap)
for name,gs in [('polyamine_synthesis_core_RNA',['ODC1','SRM','SMS']),('polyamine_SAM_supply_extended_RNA',['ODC1','SRM','SMS','AMD1','MAT2A'])]:
 x=mal[gs].to_numpy();sd=x.std(axis=0,ddof=1);assert (sd>0).all()
 metrics[name]=((x-x.mean(axis=0))/sd).mean(axis=1)
metrics['MTAP_clearance_RNA']=mal.MTAP.to_numpy()
metrics.to_csv(D/'MATCHED_MTA_PROXY_PATIENTS.tsv',sep='\t',index=False)
names=['direct_MTA_response_RNA_proxy','polyamine_synthesis_core_RNA','polyamine_SAM_supply_extended_RNA','MTAP_clearance_RNA'];rows=[]
for name in names:
 for adj in ['raw','study_adjusted']:
  r=bridge(metrics.TYMP_pseudobulk,metrics[name],metrics.C4_fraction,metrics.study,adj)
  r.update(bootstrap_ab(metrics.TYMP_pseudobulk,metrics[name],metrics.C4_fraction,metrics.study,adj))
  r.update(metric=name,adjustment=adj);rows.append(r)
out=pd.DataFrame(rows)
for p in ['a_P','XY_P','b_P']:
 out[p+'_FDR_MTA_four']=np.nan
 for _,ix in out.groupby('adjustment').groups.items():out.loc[ix,p+'_FDR_MTA_four']=false_discovery_control(out.loc[ix,p].to_numpy())
out.to_csv(D/'MTA_PROXY_C4_BRIDGE_RESULTS.tsv',sep='\t',index=False)
pairs=[];summ=[]
for scope in ['tumor_only','all_tissues']:
 m=cell[cell.role!='Malignant']
 if scope=='tumor_only':m=m[m.source_tissue.str.lower()=='tumor']
 for patient,g in m.groupby('curated_patient_id'):
  a=g[g.role=='C4'];b=g[g.role=='OtherMacrophage']
  if min(len(a),len(b))<20:continue
  pairs.append(dict(scope=scope,patient=patient,study=g.study.iloc[0],C4_cells=len(a),other_cells=len(b),effect=float(a.MTA_response_transcript_score.mean()-b.MTA_response_transcript_score.mean())))
 for_scope=pd.DataFrame(pairs);v=for_scope[for_scope.scope==scope].effect.to_numpy()
 p=float(wilcoxon(v,method='approx',correction=True).pvalue) if np.any(v) else 1.
 summ.append(dict(scope=scope,n_patients=len(v),median_effect=float(np.median(v)),positive=int((v>0).sum()),negative=int((v<0).sum()),P=p))
s=pd.DataFrame(summ);s['FDR_two_scopes']=false_discovery_control(s.P.to_numpy())
pd.DataFrame(pairs).to_csv(D/'C4_MTA_RESPONSE_PATIENT_CONTRASTS.tsv',sep='\t',index=False);s.to_csv(D/'C4_MTA_RESPONSE_SUMMARY.tsv',sep='\t',index=False)
status={'status':'COMPLETE_MTA_PROXY_REASSESSMENT_NEW_C4','measured_MTA_C4_direct_test':'NOT_ESTIMABLE_NO_MATCHED_MTA_AND_C4','public_signature':'human fibroblast direct MTA exposure; 84/89 covered, 68 UP/16 DOWN; exploratory cross-cell-context response proxy','new_C4_test_completed':True,'client_wetlab_used':False,'figures_created':False}
(D/'STATUS.json').write_text(json.dumps(status,indent=2))
print('MTA_STATE',s.to_dict('records'))
print('MTA_BRIDGE',out[['metric','adjustment','a_rho','a_P','XY_rho','XY_P','b_partial_rho','b_P','indirect_ab','ab_low','ab_high']].to_dict('records'))
