from pathlib import Path
import importlib.util,json,datetime,hashlib,gc,traceback
import numpy as np,pandas as pd,h5py
from scipy import sparse
from anndata.io import read_elem
H=Path(__file__).resolve().parents[1];P=H.parent;R=P.parent
spec=importlib.util.spec_from_file_location('base',P/'code/run_tests.py');b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b);b.P=H
plan=json.loads((H/'ADDENDUM_BEFORE_HARMONIZED_MODELS.json').read_text());panels={k:plan[k] for k in ['C4_platform17','C4_platform15_noAPOEshared']}
read=lambda f:pd.read_csv(f,sep='\t',float_precision='round_trip')
def main():
 b.gate('harmonized_start');reg=read(P/'tables/INPUT_SECTION_REGISTER.tsv');effects=[];elig=[];design=[];ident=[]
 zz=np.load(R/'phase200_union_C4_C10_C15_reassessment_20260918_v1/core/contrasts/tumor_C4.npz',allow_pickle=True);gg=list(zz['genes'])
 for name,genes in panels.items():
  jj=[gg.index(g) for g in genes];a=zz['target_log2CPM'][:,jj].mean(1);c=zz['reference_log2CPM'][:,jj].mean(1)
  for i in range(len(a)):ident.append(dict(patient=str(zz['patients'][i]),study=str(zz['study'][i]),program=name,C4_score=float(a[i]),other_macro_score=float(c[i]),difference=float(a[i]-c[i])))
 b.tab(ident,'PLATFORM_SIGNATURE_PAIRED_DESCRIPTIVE.tsv')
 ref=R/'phase219_spatial_reference_identifiability_20260920_v1/reference';rm=read(ref/'PATIENT_STATE_METADATA.tsv');rx=np.load(ref/'PATIENT_STATE_COUNTS.npz');gg=list(rx['genes']);rl=np.log2(rx['counts']/rx['counts'].sum(1)[:,None]*1e6+1)
 rr=[]
 for name,genes in panels.items():
  v=rl[:,[gg.index(g) for g in genes]].mean(1)
  for i,row in rm.iterrows():rr.append(dict(program=name,**row.to_dict(),score=float(v[i])))
 rd=pd.DataFrame(rr);b.tab(rd,'REFERENCE_LINEAGE_DESCRIPTIVE.tsv');b.tab(rd.groupby(['program','study','label']).agg(patients=('patient','nunique'),median_score=('score','median')).reset_index(),'REFERENCE_LINEAGE_SUMMARY.tsv')
 for i,r in enumerate(reg.itertuples()):
  b.gate(r.sample);b.status('RUNNING',current=r.sample,completed=i,total=38)
  z=read(P/'inputs'/f'{r.sample}_MODEL_INPUT.tsv.gz');a=read(P/'inputs'/f'{r.sample}_SOURCE_MEMBERSHIP.tsv.gz')
  with h5py.File(r.input_path,'r') as f:
   var=read_elem(f['var']);obs=read_elem(f['obs']);genes=panels['C4_platform17'];assert all(g in var.index for g in genes)
   X=read_elem(f['layers/counts'] if 'counts' in f.get('layers',{}) else f['X']).astype(float);total=np.asarray(X.sum(1)).ravel();raw=X[:,var.index.get_indexer(genes)].toarray()
  scores={}
  for name,panel in panels.items():
   ix=[genes.index(g) for g in panel];v=np.log1p(raw[:,ix]/total[:,None]*10000);sd=v.std(0);sc=((v-v.mean(0))/np.where(sd>0,sd,1)).mean(1);det=(raw[:,ix]>0).sum(1)>=3;ok=det.sum()>=20 and np.std(sc)>0
   scores[name]=sc if ok else np.full(len(obs),np.nan);elig.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,program=name,genes=';'.join(panel),genes_present=len(panel),detected_spots=int(det.sum()),eligible=ok))
   z[name]=np.nan
   for sid,g in a.groupby('segment'):z.at[sid,name]=float(scores[name][g.spot_index.to_numpy(int)].mean())
  pd.DataFrame(dict(spot_id=obs.index,**scores)).to_csv(H/'inputs'/f'{r.sample}_PROGRAMS.tsv.gz',sep='\t',index=False);z.to_csv(H/'inputs'/f'{r.sample}_MODEL_INPUT.tsv.gz',sep='\t',index=False)
  def add(family,ep,model,e,n,why='',typ='rho'):
   effects.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,family=family,endpoint=ep,model=model,effect_type=typ,effect=e,eligible_segments=n,reason=why))
  for source in ['cholesterol_disjoint','polyamine_MTA_context','TYMP_RNA']:
   for model in ['unadjusted','core_adjusted']:
    e,n,why=b.fit(z[source].to_numpy(),z.C4_platform17.to_numpy(),cov=None if model=='unadjusted' else z[b.CORE].to_numpy(),mask=z.source_eligible.to_numpy(bool),min_n=50)
    add('metabolic_context',source+'__C4_platform17',model,e,n,why)
  for M,family in [('C4_platform17','C4_primary_interaction'),('C4_platform15_noAPOEshared','C4_disjoint_sensitivity')]:
   for F in ['POSTNCAF','ECM4']:
    for model,cov in [('core_adjusted',b.CORE),('state_adjusted',b.CORE+['SPP1_share','APOEMph_noC4','source_endothelial'])]:
     for target in ['CD8_RNA','T_abundance_log2']:
      o=b.interaction(z,M,F,target+'_delta',cov);ep=M+'x'+F+'__'+target;add(family,ep,model,o['effect'],o['n'],o['reason'],'interaction_beta')
      design.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,family=family,endpoint=ep,model=model,macrophage=M,stroma=F,target=target,covariates=';'.join(cov),**o))
  print('HARMONIZED_SECTION',i+1,r.sample,flush=True);del X;gc.collect()
 b.tab(elig,'PROGRAM_ELIGIBILITY.tsv');b.tab(design,'INTERACTION_MODEL_DETAILS.tsv');sec,pat,su=b.summarize(effects)
 b.dump(H/'qc/COMPUTATION_COMPLETED.json',dict(at=b.now(),section_effects=len(sec),patient_effects=len(pat),cohort_endpoints=len(su),addendum_sha256=b.sha(H/'ADDENDUM_BEFORE_HARMONIZED_MODELS.json'),figures=0))
 b.status('COMPUTED_AWAITING_VALIDATION',completed=38,total=38);print('HARMONIZED_COMPLETE',len(su),flush=True)
if __name__=='__main__':
 try:main()
 except Exception:b.status('FAILED',traceback=traceback.format_exc());raise
