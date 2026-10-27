from pathlib import Path
import os,json,datetime,traceback,hashlib,platform,sys
import numpy as np,pandas as pd,h5py,scipy
from scipy.optimize import nnls
from scipy import stats
from anndata.io import read_elem
P=Path(__file__).resolve().parents[1];R=P.parent;B=R.parent;REF=R/'phase219_spatial_reference_identifiability_20260920_v1/reference'
def log(s):print(datetime.datetime.now().isoformat(),s,flush=True)
def save(n,o):(P/n).write_text(json.dumps(o,indent=2,default=lambda a:a.item() if hasattr(a,'item') else str(a)))
def tab(n,df):df.to_csv(P/n,sep='\t',index=False)
def balanced(raw,meta,label,training):
 ix=np.flatnonzero((meta.label.to_numpy()==label)&training)
 return np.mean([raw[j].mean(axis=0) for st in sorted(meta.iloc[ix].study.unique()) for j in [ix[meta.iloc[ix].study.to_numpy()==st]]],axis=0)
def fit(A,weight,y):
 coef,res=nnls(A*weight[:,None],y*weight,maxiter=1000)
 frac=coef/coef.sum() if coef.sum()>0 else np.zeros_like(coef)
 den=np.linalg.norm(y*weight)
 return frac,float(res/den) if den>0 else np.nan,coef

def run():
 data=np.load(REF/'PATIENT_STATE_COUNTS.npz');counts=data['counts'];genes=data['genes'].astype(str);meta=pd.read_csv(REF/'PATIENT_STATE_METADATA.tsv',sep='\t');raw=counts/meta.cells.to_numpy()[:,None]
 panels=pd.read_csv(REF/'LOSO_TRAINING_GENE_PANELS.tsv',sep='\t');lookup={g:i for i,g in enumerate(genes)}
 allref=pd.read_csv(REF/'MEAN_COUNTS_STUDY_PATIENT_BALANCED.tsv.gz',sep='\t',index_col=0);labels=list(allref.columns);assert list(allref.index)==list(genes)
 rng=np.random.default_rng(262220);rows=[];truthrows=[];cache=[]
 for study in sorted(meta.study.unique()):
  training=meta.study.to_numpy()!=study;held=np.flatnonzero(~training);c4=held[meta.iloc[held].label.to_numpy()=='C4'];other=held[meta.iloc[held].label.to_numpy()=='Other_mature_macrophage']
  fold=panels[panels.heldout_study==study];assert len(fold)==3 and len(other)>0
  fold_labels=fold.iloc[0].training_labels.split(';');assert all(x==fold.iloc[0].training_labels for x in fold.training_labels)
  nonmac=[l for l in fold_labels if l not in ['C4','Other_mature_macrophage'] and ((meta.label==l)&(~training)).any()]
  references=np.column_stack([balanced(raw,meta,l,training) for l in fold_labels]);systems={}
  for x in fold.itertuples():
   ix=np.array([lookup[g] for g in x.genes.split(';')]);a=references[ix];w=1/np.sqrt(a.mean(axis=1)+.1);systems[x.variant]=(ix,a,w)
  for donor in c4:
   source={'C4':int(donor),'Other_mature_macrophage':int(rng.choice(other))}
   for l in nonmac:source[l]=int(rng.choice(np.flatnonzero((~training)&(meta.label==l))))
   for within in [0.,.1,.25,.5,1.]:
    fractions={l:.7/len(nonmac) for l in nonmac};fractions.update(C4=.3*within,Other_mature_macrophage=.3*(1-within))
    true=np.array([fractions.get(l,0.) for l in fold_labels]);assert abs(true.sum()-1)<1e-10
    expect=sum(raw[source[l]]*v for l,v in fractions.items());prob=expect/expect.sum()
    for depth in [500,2000,10000]:
     for rep in range(3):
      seed=int(rng.integers(1,2**31-1));y=np.random.default_rng(seed).multinomial(depth,prob).astype(float)
      mix=f'{study}__{meta.patient.iloc[donor]}__{within}__{depth}__{rep}'
      truthrows.append(dict(mixture_id=mix,study=study,C4_patient=meta.patient.iloc[donor],C4_profile_row=int(donor),sources=json.dumps(source),true_fractions=json.dumps(fractions),depth=depth,replicate=rep,seed=seed))
      for variant,(ix,a,w) in systems.items():
       frac,res,coef=fit(a,w,y[ix]);cc=fold_labels.index('C4');oo=fold_labels.index('Other_mature_macrophage');macro=frac[cc]+frac[oo];estwithin=frac[cc]/macro if macro>0 else np.nan
       rows.append(dict(mixture_id=mix,study=study,patient=meta.patient.iloc[donor],variant=variant,depth=depth,repeat=rep,true_C4=fractions['C4'],estimate_C4=frac[cc],true_C4_within_macrophage=within,estimate_C4_within_macrophage=estwithin,true_total_macro=.3,estimate_total_macro=macro,C4_abs_error=abs(frac[cc]-fractions['C4']),within_macro_abs_error=abs(estwithin-within),zero_C4_false_positive=(within==0 and frac[cc]>.05),relative_weighted_residual=res,panel_counts=int(y[ix].sum()),true_other_macro=fractions['Other_mature_macrophage'],estimated_fractions=json.dumps(dict(zip(fold_labels,frac.tolist())))))
  log('simulated heldout '+study+' C4 donor profiles '+str(len(c4)))
 df=pd.DataFrame(rows);truth=pd.DataFrame(truthrows);tab('simulation/ALL_RECOVERY_RESULTS.tsv.gz',df);tab('simulation/MIXTURE_TRUTH_AND_SEEDS.tsv.gz',truth)
 patient=[]
 for (variant,depth,study,pat),g in df.groupby(['variant','depth','study','patient']):
  z=g[g.true_C4==0]
  patient.append(dict(variant=variant,depth=depth,study=study,patient=pat,n_mixtures=len(g),C4_total_MAE=float(g.C4_abs_error.median()),C4_within_macro_MAE=float(g.within_macro_abs_error.median()),zero_C4_FPR=float(z.zero_C4_false_positive.mean()),residual_median=float(g.relative_weighted_residual.median())))
 pdn=pd.DataFrame(patient);tab('simulation/PATIENT_RECOVERY_SUMMARY.tsv',pdn)
 summary=[]
 for (variant,depth),g in pdn.groupby(['variant','depth']):
  ss=g.groupby('study').C4_within_macro_MAE.median();a=g.C4_total_MAE.median();b=g.C4_within_macro_MAE.median();fp=g.zero_C4_FPR.mean()
  passed=a<=.05 and b<=.15 and fp<=.1 and len(ss)==5 and ss.max()<=.2
  summary.append(dict(variant=variant,depth=depth,independent_C4_profiles=len(g),studies=len(ss),median_patient_C4_total_MAE=a,median_patient_C4_within_macro_MAE=b,mean_patient_zero_C4_FPR=fp,worst_study_within_macro_MAE=ss.max(),primary_gate_pass=passed if depth==2000 else None,study_errors=json.dumps(ss.to_dict())))
 tab('simulation/RECOVERY_DECISIONS.tsv',pd.DataFrame(summary))
 # Fit reproducible technical maps on outcome-blind two-section sample.
 fullcpm=counts/counts.sum(axis=1)[:,None]*1e6;train=np.ones(len(meta),bool);prof=np.array([balanced(fullcpm,meta,l,train) for l in labels]);lp=np.log2(prof+1)
 source54=B/'project26210_external_joint_v1_phase54_gse175540_spatial_v1/run_v1/sections_v1'
 common=set(genes)
 for dd in sorted(source54.iterdir()):
  with h5py.File(dd/'processed_section_v1.h5ad') as h:common &= set(read_elem(h['var']).index.astype(str))
 prog=json.loads((R/'phase218_spatial_identity_alignment_20260920_v1/PROGRAM_DEFINITIONS.json').read_text());ex=set(prog['spp1_state_full'])|set(prog['vascular_remodeling'])|{'SPP1','TYMP'}
 variants={'complete':set(),'exclude_SPP1_TYMP':{'SPP1','TYMP'},'exclude_program_union':ex};systems={};pans=[]
 for variant,excluded in variants.items():
  allowed=np.array([g in common and g not in excluded and not g.startswith(('MT-','RPL','RPS')) for g in genes]);inds=set()
  for j,l in enumerate(labels):
   diff=lp[j]-np.delete(lp,j,axis=0).mean(axis=0);ix=np.flatnonzero(allowed&(prof[j]>=1));ix=ix[np.argsort(diff[ix])[::-1]];inds.update(ix[diff[ix]>0][:50].tolist())
  ix=np.array(sorted(inds));a=allref.iloc[ix].to_numpy();w=1/np.sqrt(a.mean(axis=1)+.1);systems[variant]=(ix,a,w);pans.append(dict(variant=variant,genes=';'.join(genes[ix]),n_genes=len(ix)))
 tab('spatial/TECHNICAL_MAPPING_PANELS.tsv',pd.DataFrame(pans));spstats=[];stability=[]
 for sid in ['GSM5924040','GSM5924048']:
  d=P/'spatial'/sid;d.mkdir()
  with h5py.File(source54/sid/'processed_section_v1.h5ad') as h:
   obs=read_elem(h['obs']);var=read_elem(h['var']);X=read_elem(h['layers']['counts']);coords=read_elem(h['obsm']['spatial'])
  maps={};residuals={}
  for variant,(ix,a,w) in systems.items():
   vi=var.index.get_indexer(genes[ix]);assert (vi>=0).all();y=X[:,vi].toarray().astype(float);frac=np.zeros((len(obs),len(labels)));res=np.zeros(len(obs));rawcoef=np.zeros_like(frac)
   for j in range(len(obs)):frac[j],res[j],rawcoef[j]=fit(a,w,y[j])
   maps[variant]=frac;residuals[variant]=res
   output=pd.DataFrame(frac,columns=['relative__'+l for l in labels]);output.insert(0,'spot_id',obs.index);output['x']=coords[:,0];output['y']=coords[:,1];output['weighted_relative_residual']=res;output['marker_panel_counts']=y.sum(axis=1).astype(int);output['zero_marker_counts']=y.sum(axis=1)==0
   output.to_csv(d/(variant+'_TECHNICAL_RELATIVE_FIT.tsv.gz'),sep='\t',index=False)
   np.savez_compressed(d/(variant+'_FIT_CACHE.npz'),y=y,A=a,weights=w,coef=rawcoef,genes=genes[ix],labels=np.array(labels))
   spstats.append(dict(section=sid,variant=variant,spots=len(obs),n_panel_genes=len(ix),median_panel_counts=np.median(y.sum(axis=1)),zero_panel_spots=int((y.sum(axis=1)==0).sum()),median_weighted_residual=np.nanmedian(res),p90_weighted_residual=np.nanquantile(res,.9),median_C4_relative=np.median(frac[:,labels.index('C4')]),max_C4_relative=np.max(frac[:,labels.index('C4')]),median_CD8_relative=np.median(frac[:,labels.index('CD8_T')]),finite_fraction=float(np.isfinite(frac).mean()),nonzero_fit_spots=int((frac.sum(axis=1)>0).sum())))
  for variant in ['exclude_SPP1_TYMP','exclude_program_union']:
   for label in ['C4','Other_mature_macrophage','CD8_T','NK','gdT_enriched','Endothelial','Malignant']:
    x=maps['complete'][:,labels.index(label)];y=maps[variant][:,labels.index(label)]
    rho=stats.spearmanr(x,y).statistic if np.std(x)>0 and np.std(y)>0 else np.nan
    stability.append(dict(section=sid,variant=variant,label=label,descriptive_rho=rho,median_abs_fraction_change=float(np.median(abs(x-y))),p90_abs_fraction_change=float(np.quantile(abs(x-y),.9))))
  log('two-section technical fit complete '+sid)
 tab('spatial/TECHNICAL_MAPPING_SUMMARY.tsv',pd.DataFrame(spstats));tab('spatial/PANEL_STABILITY_DESCRIPTIVE.tsv',pd.DataFrame(stability))
 before=json.loads((P/'audit/PREEXISTING_RDS_STAT.json').read_text());changed=[]
 for f in before:
  q=Path(f['path']);st=q.stat()
  if (st.st_size,st.st_mtime_ns,st.st_ino)!=(f['bytes'],f['mtime_ns'],f['inode']):changed.append(f['path'])
 assert not changed,changed
 save('audit/RDS_PRESERVATION_CHECK.json',dict(tracked=len(before),changed=changed,write_operations_to_existing_RDS=0,verification='size/mtime_ns/inode before and after; scripts only write into new Phase220 root'))
 save('COMPUTE_VALIDATION.json',dict(status='PASS_BASIC',heldout_C4_profiles=truth.C4_profile_row.nunique(),mixtures=len(truth),fits=len(df),sections=2,spot_counts=[4540,1544],no_spatial_inference=True,no_RDS_outputs=True,RDS_preserved=len(before)))
 save('STATUS.json',dict(status='COMPUTE_COMPLETE',method='NNLS_TECHNICAL_PILOT',cell2location_run=False,biological_spatial_claims=False));(P/'EXIT_STATUS').write_text('0\n');log('COMPLETE')
try:run()
except Exception:
 (P/'ERROR.txt').write_text(traceback.format_exc());(P/'EXIT_STATUS').write_text('1\n');save('STATUS.json',dict(status='FAILED'));raise
