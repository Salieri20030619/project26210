from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
from scipy.stats import rankdata,t,norm
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');O=P/'phase162_MTA_polyamine_part1_reassessment_20260915_v1'
R=pd.read_csv(O/'PUBLIC_COHORT_AXIS_ASSOCIATIONS.tsv',sep='\t');M=pd.read_csv(O/'PUBLIC_META_AXIS_ASSOCIATIONS.tsv',sep='\t');checks=[]
for cohort,g in R.groupby('cohort'):
 f=pd.read_csv(O/('PUBLIC_PATIENTS_'+cohort+'.tsv'),sep='\t',index_col=0)
 for i,r in g.iterrows():
  cov=[] if r.adjustment=='raw' else ['MYELOID_PROXY']
  if r.edge=='a':cols=['TYMP',r.gene]+cov
  elif r.edge=='XY':cols=[r.gene,r.outcome]+cov
  else:cols=[r.gene,r.outcome,'TYMP']+(['MTAP'] if r.edge=='b_adjust_MTAP' else [])+cov
  d=f[cols].dropna();mat=np.array([rankdata(d[c]) for c in cols]);prec=np.linalg.inv(np.corrcoef(mat));rho=-prec[0,1]/np.sqrt(prec[0,0]*prec[1,1]);k=len(cols)-2;p=2*t.sf(abs(rho)*np.sqrt((len(d)-k-2)/(1-rho*rho)),len(d)-k-2)
  checks.append(dict(check='precision_'+str(i),passed=abs(rho-r.rho)<1e-9 and abs(p-r.P)<1e-9))
for i,r in M.iterrows():
 g=R[(R.gene==r.gene)&(R.adjustment==r.adjustment)&(R.outcome==r.outcome)&(R.edge==r.edge)]
 z=np.arctanh(g.rho.to_numpy());v=1/(g.n.to_numpy()-g.k.to_numpy()-3);w=1/v;mu=np.average(z,weights=w);q=np.sum(w*(z-mu)**2);tau=max(0,(q-len(g)+1)/(w.sum()-(w*w).sum()/w.sum()));wr=1/(v+tau);m=np.average(z,weights=wr)
 checks.append(dict(check='meta_'+str(i),passed=abs(np.tanh(m)-r.random_rho)<1e-10))
for _,ix in M.groupby(['family','adjustment','outcome','edge']).groups.items():
 q=multipletests(M.loc[ix,'random_P'],method='fdr_bh')[1];checks.append(dict(check='meta_bh_'+str(min(ix)),passed=np.allclose(q,M.loc[ix,'random_FDR'])))
V=pd.DataFrame(checks);assert V.passed.all();V.to_csv(O/'INDEPENDENT_VALIDATION.tsv',sep='\t',index=False)
status=json.loads((O/'STATUS.json').read_text());status.update(independent_checks=len(V),independent_passed=int(V.passed.sum()),total_checks=status['checks']+len(V));(O/'STATUS.json').write_text(json.dumps(status,indent=2))
# Copy relevant first-part and historical results for exact provenance in the version directory.
import shutil
Q=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1');old=O/'FIRST_PART_AND_LEGACY_EVIDENCE';old.mkdir(exist_ok=True)
files=[Q/'09_metabolism_trajectory/metabolism_v1/tables/four_pathway_AUCell_paired_statistics_v1.tsv',Q/'09_metabolism_trajectory/metabolism_v1/tables/scFEA_four_axis_paired_statistics_v1.tsv',Q/'09_metabolism_trajectory/metabolism_v1/tables/scFEA_scope_limitations_v1.tsv',Path('/mnt/BioDisk/wangmingyuan/project26210_phase79_mtap_natural_perturbation_v1/results_v1/tcga_kirc_interaction_primary_v1.tsv')]
for p in files:shutil.copy2(p,old/p.name)
sources=files+[Q/'07_pyscenic_v1/input/raw_counts_genes_by_cells_v1.mtx.gz',Q/'07_pyscenic_v1/input/genes_v1.tsv',Q/'07_pyscenic_v1/input/cells_v1.tsv',Q/'09_metabolism_trajectory/metabolism_v1/source_data/four_pathway_AUCell_cell_scores_v1.tsv.gz']
sources+=list(Path('/mnt/BioDisk/wangmingyuan/project26210_phase98_public_paired_metabolomics_v1/downloads_v1').glob('matched_tpm_*.csv'))
(O/'SOURCE_MANIFEST.json').write_text(json.dumps([dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sources],indent=2))
print('MTA RAW\n'+M[(M.outcome=='MET:5-methylthioadenosine (MTA)')&(M.adjustment=='raw')][['gene','edge','random_rho','random_P','random_FDR','I2','positive_cohorts']].to_string(index=False))
print('SPP1 MYELOID\n'+M[(M.outcome=='SPP1')&(M.adjustment=='myeloid')&(M.edge!='a')][['gene','edge','random_rho','random_P','random_FDR','positive_cohorts']].to_string(index=False))
print('STATUS',json.dumps(status))

