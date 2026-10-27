from pathlib import Path
import json,hashlib,subprocess,math,datetime,gc
import numpy as np,pandas as pd,h5py
from scipy import sparse,stats
from anndata.io import read_elem
P=Path(__file__).resolve().parents[1];R=P.parent
reg=pd.read_csv(P/'resources/EXTERNAL_SECTION_REGISTER.tsv',sep='\t');sec=pd.read_csv(P/'tables/SECTION_EFFECTS_ALL.tsv',sep='\t',float_precision='round_trip');pa=pd.read_csv(P/'tables/PATIENT_EFFECTS_ALL.tsv',sep='\t',float_precision='round_trip');su=pd.read_csv(P/'tables/COHORT_SUMMARY_ALL.tsv',sep='\t',float_precision='round_trip');checks=[]
def ck(n,b,d=None):checks.append({'check':n,'pass':bool(b),'detail':d})
def sha(f):
 h=hashlib.sha256()
 with open(f,'rb') as x:
  for b in iter(lambda:x.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def bo(s):return s.astype(str).str.lower().isin(['true','1'])
ck('frozen_protocol_sha',sha(P/'resources/EXTERNAL_PROTOCOL_FROZEN.json')==(P/'resources/EXTERNAL_PROTOCOL_FROZEN.sha256').read_text())
features=json.loads((P/'resources/FEATURES_FROM_PAPER.json').read_text());ck('disjoint_state_features',not(set(features['FAP_CAF25_disjoint'])&set(features['SPP1_MAC25_disjoint'])))
for n in ['ECM_org_heldout','Collagen_org_heldout','TGFb_response_heldout']:ck(n+'_excludes_state_features',not(set(features[n])&(set(features['FAP_CAF25'])|set(features['SPP1_MAC25']))))
maxscore=0;maxeffect=0;complete_source=[]
for a in reg.itertuples():
 subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],check=True)
 d=P/'inputs'/a.sample;info=json.loads((d/'INPUT_QC.json').read_text());m=pd.read_csv(d/'METADATA.tsv',sep='\t',float_precision='round_trip');sc=pd.read_csv(d/'MODULE_SCORES.tsv',sep='\t',float_precision='round_trip');member=json.loads((d/'SCORING_MEMBERSHIP.json').read_text());gn=pd.read_csv(d/'GENE_NAME_MAPPING.tsv',sep='\t')
 ck(a.sample+'_source_sha',sha(a.input_path)==info['input_sha256'])
 old=R/'phase289_state_source_boundary_process_20260929_v1/inputs'/f'{a.sample}_SPOT_PROGRAMS.tsv.gz';ck(a.sample+'_old_labels_sha',sha(old)==info['old_metadata_sha256'])
 with h5py.File(a.input_path,'r') as h:
  gg=np.asarray(read_elem(h['var/gene_symbol_original'])) if 'gene_symbol_original' in h['var'] else np.asarray(read_elem(h['var'][h['var'].attrs['_index']]));gg=gg.astype(str)
  ids=np.asarray(read_elem(h['obs'][h['obs'].attrs['_index']])).astype(str);idx=pd.Index(ids).get_indexer(m.spot_id);assert (idx>=0).all()
  hh=h[info['layer']];X=sparse.csr_matrix((hh['data'][:],hh['indices'][:],hh['indptr'][:]),shape=tuple(hh.attrs['shape']))[idx].astype(float)
 if pd.Index(gg).has_duplicates:
  ug,inv=np.unique(gg,return_inverse=True);X=(X@sparse.csr_matrix((np.ones(len(gg)),(np.arange(len(gg)),inv)),shape=(len(gg),len(ug)))).tocsr();gg=ug
 ck(a.sample+'_matrix_genes_match',np.array_equal(gg,gn.input_gene.to_numpy()))
 totals=np.asarray(X.sum(1)).ravel();ck(a.sample+'_library_totals',np.array_equal(totals,m.library_total.to_numpy()))
 norm=X.multiply((10000/totals)[:,None]).tocsr();norm.data=np.log1p(norm.data);ix={g:i for i,g in enumerate(gn.scoring_gene)}
 recomputed={}
 for name in member['features']:
  fg=member['features'][name];cg=member['controls'][name];ck(a.sample+'_'+name+'_coverage',len(fg)>=max(5,math.ceil(.8*len(features[name]))))
  ck(a.sample+'_'+name+'_unique_controls',len(cg)==len(set(cg)))
  v=np.asarray(norm[:,[ix[g] for g in fg]].mean(1)).ravel()-np.asarray(norm[:,[ix[g] for g in cg]].mean(1)).ravel();err=float(np.max(abs(v-sc[name])));maxscore=max(maxscore,err);ck(a.sample+'_'+name+'_score_from_original',err<1e-10,err);recomputed[name]=v
 z=m.copy()
 for name,v in recomputed.items():z[name]=v
 for gene in ['FAP','SPP1','CD3D','CD8A']:
  if gene in gg:
   raw_gene=np.asarray(X[:,np.flatnonzero(gg==gene)[0]].toarray()).ravel();vals=np.log1p(raw_gene/totals*10000);ck(a.sample+'_'+gene+'_expression',np.max(abs(vals-m['g_'+gene]))<1e-12);z['g_'+gene]=vals;z['det_'+gene]=vals>0
 origmeta=pd.read_csv(old,sep='\t',float_precision='round_trip').set_index('spot_id').loc[m.spot_id]
 for n in ['boundary_input_valid','inside_morphology_supported','inside_expression_sensitivity']:ck(a.sample+'_'+n+'_unchanged',np.array_equal(bo(m[n]),bo(origmeta[n])))
 ck(a.sample+'_mapped_share_original',np.allclose(m.mapped_SPP1_share,origmeta.SPP1_share,equal_nan=True))
 for row in sec[sec['sample'].eq(a.sample)].itertuples():
  if pd.isna(row.effect):continue
  if row.family in ['coupling','process']:
   xname,yname=row.endpoint.split('__');x=np.asarray(z[xname]);y=np.asarray(z[yname]);valid=np.isfinite(x)&np.isfinite(y)
   if row.model=='depth_adjusted':valid&=np.isfinite(z.log1p_UMI)
   xv=stats.rankdata(x[valid]);yv=stats.rankdata(y[valid])
   if row.model=='depth_adjusted':
    dep=stats.rankdata(z.log1p_UMI.to_numpy()[valid]);dep=(dep-dep.mean())/dep.std();Q,_=np.linalg.qr(np.column_stack([np.ones(len(dep)),dep]));xv=xv-Q@(Q.T@xv);yv=yv-Q@(Q.T@yv)
   ef=float(np.corrcoef(xv,yv)[0,1])
  else:
   gene,typ=row.endpoint.split('__');valid=bo(z.boundary_input_valid).to_numpy();inn=valid&bo(z['inside_'+row.model]).to_numpy();out=valid&~inn;v=np.asarray(z[('g_' if typ=='delta_mean_log_expression' else 'det_')+gene],float)
   ef=math.fsum(v[inn])/int(inn.sum())-math.fsum(v[out])/int(out.sum())
  err=abs(ef-row.effect);maxeffect=max(maxeffect,err);ck(a.sample+'_'+row.family+'_'+row.endpoint+'_'+row.model+'_effect',err<1e-9,err)
 complete_source.append(a.sample);print('VALIDATED_SOURCE',a.sample,flush=True);del X,norm,z;gc.collect()
subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],check=True)
keys=['dataset','stratum','family','endpoint','model','effect_type','expected_direction']
for a in pa.itertuples():
 q=sec[sec.dataset.eq(a.dataset)&sec.patient.eq(a.patient)&sec.family.eq(a.family)&sec.endpoint.eq(a.endpoint)&sec.model.eq(a.model)]
 if a.stratum!='overall':q=q[q.preparation.eq(a.stratum.removeprefix('preparation_'))]
 v=q.effect.dropna().to_numpy();expected=math.tanh(math.fsum(math.atanh(float(np.clip(x,-1+1e-12,1-1e-12))) for x in v)/len(v)) if len(v) and a.effect_type=='spearman_rho' else math.fsum(v)/len(v) if len(v) else np.nan
 ck('patient_'+str(a.Index), (pd.isna(expected) and pd.isna(a.effect)) or abs(expected-a.effect)<1e-12)
for a in su.itertuples():
 q=pa[(pa.dataset==a.dataset)&(pa.stratum==a.stratum)&(pa.family==a.family)&(pa.endpoint==a.endpoint)&(pa.model==a.model)];v=q.effect.dropna().to_numpy();pos=sum(v>0);neg=sum(v<0);n=len(v);nt=pos+neg
 expected=min(1.,2*sum(math.comb(nt,k) for k in range(min(pos,neg)+1))/(2**nt)) if n>=10 and nt else 1. if n>=10 else np.nan
 ck('summary_n_'+str(a.Index),n==a.n_patients and pos==a.n_positive and neg==a.n_negative)
 ck('summary_P_'+str(a.Index),(pd.isna(expected) and pd.isna(a.p_raw)) or abs(expected-a.p_raw)<1e-14)
 ck('summary_median_'+str(a.Index),(not n and pd.isna(a.median_effect)) or (n and abs(float(np.median(v))-a.median_effect)<1e-12))
for kk,grp in su.groupby(['dataset','stratum','family']):
 vv=grp.p_raw.fillna(1).to_numpy();expected=stats.false_discovery_control(vv,method='bh');ok=grp.p_raw.notna().to_numpy();ck('BH_'+str(kk),np.allclose(expected[ok],grp.FDR.to_numpy()[ok],atol=1e-14,rtol=0))
ck('n_lt10_no_inference',su.loc[su.n_patients<10,['p_raw','FDR']].isna().all().all())
ck('section_rows_38x28',len(sec)==38*28);ck('unique_patients_37',reg[['dataset','patient']].drop_duplicates().shape[0]==37)
ref=pd.read_csv(P/'qc/REFERENCE_CALIBRATION.tsv',sep='\t');ck('failed_reference_gate_preserved',not ref['pass'].any())
pd.DataFrame(checks).to_csv(P/'qc/FINAL_VALIDATION_CHECKS.tsv',sep='\t',index=False)
out={'status':'PASS' if all(c['pass'] for c in checks) else 'FAIL','checks':len(checks),'failed':[c for c in checks if not c['pass']],'source_sections_validated':len(complete_source),'max_score_error':maxscore,'max_effect_error':maxeffect,'time':datetime.datetime.now().astimezone().isoformat(),'meaning':'Technical reproduction and complete reporting validated; biological validation assessed endpoint-wise; exact original score gate failed'}
(P/'qc/FINAL_VALIDATION.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2),flush=True)
if out['status']!='PASS':raise SystemExit(1)
