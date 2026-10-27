from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='4'
import json,subprocess,gc
import numpy as np,pandas as pd,anndata as ad
P=Path(__file__).resolve().parents[1];R=P.parent
rule={'OBRADOVIC':{'macrophage':[2,3,5,6,8,9,10,15,17],'T_candidates':[4,5,6,7,8,9,10],'tumour_like':[13,14,15,16,22,23,24,26],'stroma':[21],'primary_immune_fraction':'CD45+','primary_tumour_stroma_fraction':'CD45-'},'HCRN':{'macrophage':[1,2,4,5,9,13,16,17],'T_candidates':[0,1,2,5,6,7,8,16],'tumour_like':[4,11,17,19,22,23],'stroma':[13]},'macrophage_gate':'Reviewed refined clusters; >=2 macrophage support genes; exclude >=2 T lineage OR >=2 epithelial OR >=2 endothelial markers. DC/monocyte/osteoclast/clear mixed clusters excluded.','CD8_gate':'Reviewed T/NK candidate clusters; >=2 of CD3D CD3E TRAC and >=1 CD8A CD8B; exclude >=2 epithelial/endothelial and >=2 myeloid markers. Treg clusters excluded.','n_gates':{'Macrophage':100,'CD8':20,'Tumour_like':50,'Fibro_strict_candidate':20,'Stromal_comarker':20,'Mural':50},'patient_selection':'OBRADOVIC same-patient tumour CD45+ fractions aggregated; companion CD45- tumour/stroma fractions matched at patient level. HCRN baseline/progression kept separate; donor_first_available selects baseline if present irrespective of cell or target eligibility; no duplicate donor in pooled summaries.','association_tests_run':False}
f=P/'resources/LINEAGE_IDENTITY_FROZEN.json'
if not f.exists():f.write_text(json.dumps(rule,indent=2))
else:assert json.loads(f.read_text())==rule
pan=json.loads((R/'phase286_CAF_identity_boundary_audit_20260929_v1/resources/IDENTITY_PANELS_FROZEN.json').read_text());rows=[]
for ds in ['OBRADOVIC','HCRN']:
 subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],check=True,stdout=subprocess.DEVNULL)
 a=ad.read_h5ad(P/f'inputs/{ds}_counts_qc.h5ad');m=pd.read_csv(P/f'tables/{ds}_MYELOID_METADATA.tsv.gz',sep='\t',index_col=0);r=rule[ds];x=a.X.tocsr()
 def ng(gs):return np.asarray((x[:,a.var_names.isin(gs)]>0).sum(axis=1)).ravel()
 a.obs['myeloid_cluster']=m.myeloid_cluster.astype(str).reindex(a.obs_names).fillna('not_myeloid').to_numpy()
 cl=a.obs.cluster.astype(str);mc=a.obs.myeloid_cluster.astype(str)
 epi=ng(pan['epithelial'])>=2;endo=ng(pan['endothelial'])>=2;tl=ng(['CD3D','CD3E','TRAC']);my=ng(['LST1','TYROBP','FCER1G']);macsupport=ng(['C1QA','C1QB','C1QC','CSF1R','CD163','MRC1','CD68','APOE','TREM2'])
 immune_fraction=a.obs.selection.eq('CD45+').to_numpy() if ds=='OBRADOVIC' else np.ones(a.n_obs,dtype=bool)
 nonimmune_fraction=a.obs.selection.eq('CD45-').to_numpy() if ds=='OBRADOVIC' else np.ones(a.n_obs,dtype=bool)
 masks={'Macrophage':mc.isin([str(c) for c in r['macrophage']]).to_numpy()&(macsupport>=2)&(tl<2)&~epi&~endo&immune_fraction,'CD8':cl.isin([str(c) for c in r['T_candidates']]).to_numpy()&(tl>=2)&(ng(['CD8A','CD8B'])>=1)&(my<2)&~epi&~endo&immune_fraction,'Tumour_like':cl.isin([str(c) for c in r['tumour_like']]).to_numpy()&nonimmune_fraction}
 stromal=cl.isin([str(c) for c in r['stroma']]).to_numpy()&nonimmune_fraction;fib=(ng(pan['fibro_core'])>=2)&(ng(pan['fibro_support'])>=1);contam=(ng(pan['immune'])>=2)|epi|endo;mural=(ng(pan['pericyte'])>=2)|((ng(['MYH11','CNN1'])>=1)&(ng(['ACTA2','TAGLN'])>=1))
 masks['Fibro_strict_candidate']=stromal&fib&~contam&~mural;masks['Stromal_comarker']=stromal&fib&~contam&mural;masks['Mural']=stromal&mural&~contam
 for k,v in masks.items():a.obs['lineage_'+k]=v
 if ds=='OBRADOVIC':a.obs['stratum']='Untreated';a.obs['technical_batch']=a.obs.patient.str.contains(r'Patient[A-C]$').map({True:'batch1_v2',False:'batch2_v3'});a.obs['analysis_unit']=a.obs.patient.astype(str)
 else:a.obs['stratum']=a.obs.timepoint.astype(str);a.obs['technical_batch']='HCRN';a.obs['analysis_unit']=a.obs.patient.astype(str)+'::'+a.obs.timepoint.astype(str)
 first={pt:('Baseline' if 'Baseline' in g.stratum.astype(str).tolist() else str(g.stratum.iloc[0])) for pt,g in a.obs.groupby('patient',observed=True)}
 a.obs['donor_first_available']=[str(s)==first[str(pt)] for pt,s in zip(a.obs.patient,a.obs.stratum)]
 a.obs.to_csv(P/f'tables/{ds}_QUALIFIED_CELL_METADATA.tsv.gz',sep='\t')
 for unit,g in a.obs.groupby('analysis_unit',observed=True):
  ix=a.obs.analysis_unit.eq(unit).to_numpy();rec={'dataset':ds,'analysis_unit':unit,'patient':str(g.patient.iloc[0]),'stratum':str(g.stratum.iloc[0]),'technical_batch':str(g.technical_batch.iloc[0]),'donor_first_available':bool(g.donor_first_available.iloc[0]),'libraries':';'.join(sorted(g.library.unique().astype(str))), 'total_cells':len(g)}
  for k,v in masks.items():rec[k]=int((ix&v).sum())
  rows.append(rec)
 print('QUALIFIED',ds,flush=True);del a,x;gc.collect()
out=pd.DataFrame(rows);out.to_csv(P/'tables/LINEAGE_COUNTS_BY_PATIENT.tsv',sep='\t',index=False);print(out[['dataset','patient','stratum','Macrophage','CD8','Tumour_like','Fibro_strict_candidate','Stromal_comarker']].to_string(index=False))
(P/'qc/LINEAGE_INPUTS_COMPLETE.json').write_text(json.dumps({'status':'COMPLETE','association_tests_run':False},indent=2))
