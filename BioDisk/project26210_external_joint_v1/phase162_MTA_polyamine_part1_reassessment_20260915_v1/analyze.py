from pathlib import Path
import pandas as pd,numpy as np,json,gzip,sys,subprocess,hashlib,datetime,ast
from scipy.io import mmread
from scipy.stats import wilcoxon,spearmanr,rankdata,t,norm,chi2
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');Q=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1');O=P/'phase162_MTA_polyamine_part1_reassessment_20260915_v1';O.mkdir(exist_ok=True)
gate=subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True);(O/'MEMORY_GATE.log').write_text(gate.stdout+gate.stderr);assert gate.returncode==0
sys.path.insert(0,str(P/'phase157_ordered_bridge_reassessment_20260915_v1'));from bridge_stats import rankz,corrp
axis=['ODC1','SRM','SMS','MTAP','AMD1','MAT2A'];target=axis+['TYMP','SLC16A3']
plan=dict(firstpart='Frozen 11950 malignant cells; original TYMP-positive/negative and F1-high-clone assignments; exact prior eligible patient lists',axis=axis,upstream_supplement='SLC16A3 full Part1 cohort and F1 context, no outcome-selected subgroup',effects='patient pseudobulk differences; within-patient rank correlations raw/depth/depth+G2M, signed rank across patients; >=50 cells and >=10 detected per gene for correlations',public='4 matched observed MTA cohorts; separate production and clearance; raw/myeloid covariate; gene to MTA/SPP1 and continuous TYMP x low-MTAP interaction on MTA',MTAP_genomic='Separate deep -2 from shallow -1 and neutral 0; no calling scRNA zero expression genomic deletion',no_figures=True,no_client_data_in_models=True,fold44='Client-mentioned magnitude not independently derivable from provided figure; not a public-data target or selection rule')
(O/'PLAN_BEFORE_FIT.json').write_text(json.dumps(plan,indent=2))
inp=Q/'07_pyscenic_v1/input';genes=pd.read_csv(inp/'genes_v1.tsv',sep='\t').gene.astype(str).tolist();ids=pd.read_csv(inp/'cells_v1.tsv',sep='\t').exact_cell_id.astype(str)
with gzip.open(inp/'raw_counts_genes_by_cells_v1.mtx.gz','rb') as f:mat=mmread(f).tocsr()
assert mat.shape==(len(genes),11950)
tot=np.asarray(mat.sum(0)).ravel()
meta=pd.read_csv(Q/'09_metabolism_trajectory/metabolism_v1/source_data/four_pathway_AUCell_cell_scores_v1.tsv.gz',sep='\t').set_index('exact_cell_id').loc[ids].reset_index()
labels=pd.read_csv(Q/'05_cnmf/input_v1/cnmf_cell_metadata_v1.tsv.gz',sep='\t').set_index('exact_cell_id').loc[ids]
meta['sample_id']=labels.sample_id.to_numpy();meta['total_counts']=tot
raw=mat[[genes.index(g) for g in target]].toarray().T
E=np.log1p(raw/tot[:,None]*1e4)
assert np.allclose(raw[:,target.index('TYMP')],meta.TYMP_raw_count_v1)
gmt=P/'part2_rebuild_v1/06_functional_programs/run_v1/resources/MSigDB_Hallmark_2020.gmt';cycle=[]
for line in gmt.read_text().splitlines():
 p=line.split('\t')
 if p[0]=='G2-M Checkpoint':cycle=[g for g in p[2:] if g in genes and g not in target]
assert len(cycle)>10
CE=np.log1p(mat[[genes.index(g) for g in cycle]].toarray().T/tot[:,None]*1e4)
for j,g in enumerate(target):meta['count:'+g]=raw[:,j];meta['expr:'+g]=E[:,j]
meta.to_csv(O/'PART1_TARGET_CELL_METADATA.tsv.gz',sep='\t',index=False,compression='gzip')
prior=pd.read_csv(Q/'09_metabolism_trajectory/metabolism_v1/tables/four_pathway_AUCell_patient_deltas_v1.tsv',sep='\t')
deltas=[];coverage=[]
for contrast,group,high,low in [('TYMP-positive_vs_negative','TYMP_group','TYMP_positive','TYMP_negative'),('F1-high_clone_vs_other','clone_group','F1_high_primary_clone','other_primary_clones')]:
 eligible=prior.loc[prior.contrast==contrast,'curated_patient_id'].unique()
 for pat in eligible:
  d=meta[meta.curated_patient_id==pat];hi=d[d[group]==high];lo=d[d[group]==low]
  assert len(hi)>0 and len(lo)>0
  for gene in axis+['SLC16A3']:
   hv=np.log2(hi['count:'+gene].sum()/hi.total_counts.sum()*1e6+1);lv=np.log2(lo['count:'+gene].sum()/lo.total_counts.sum()*1e6+1)
   deltas.append(dict(patient=pat,contrast=contrast,gene=gene,high_cells=len(hi),low_cells=len(lo),high_log2CPM=hv,low_log2CPM=lv,delta=hv-lv))
D=pd.DataFrame(deltas);D.to_csv(O/'PART1_FROZEN_PATIENT_GENE_CONTRASTS.tsv',sep='\t',index=False)
ds=[]
for (contrast,gene),d in D.groupby(['contrast','gene']):
 v=d.delta;ds.append(dict(contrast=contrast,gene=gene,n_patients=len(v),positive=int((v>0).sum()),negative=int((v<0).sum()),median_delta=v.median(),P=wilcoxon(v).pvalue))
DS=pd.DataFrame(ds)
for contrast,ix in DS.groupby('contrast').groups.items():
 ia=[i for i in ix if DS.loc[i,'gene']!='SLC16A3'];DS.loc[ia,'FDR_axis6']=multipletests(DS.loc[ia,'P'],method='fdr_bh')[1]
DS.to_csv(O/'PART1_GENE_CONTRAST_SUMMARY.tsv',sep='\t',index=False)
cors=[];patient=[]
for pat,d in meta.groupby('curated_patient_id'):
 for gene in axis+['SLC16A3']:
  patient.append(dict(patient=pat,study=d.study.iloc[0],gene=gene,cells=len(d),counts=d['count:'+gene].sum(),detected=int((d['count:'+gene]>0).sum()),detection=(d['count:'+gene]>0).mean(),log2CPM=np.log2(d['count:'+gene].sum()/d.total_counts.sum()*1e6+1)))
  nx=int((d['count:TYMP']>0).sum());ny=int((d['count:'+gene]>0).sum());ok=len(d)>=50 and nx>=10 and ny>=10
  coverage.append(dict(patient=pat,gene=gene,cells=len(d),TYMP_detected=nx,gene_detected=ny,eligible=ok))
  if not ok:continue
  cy=CE[d.index];sd=cy.std(0,ddof=1);v=sd>0;cycle_score=((cy[:,v]-cy[:,v].mean(0))/sd[v]).mean(1)
  A=rankz(d['expr:TYMP']);Y=rankz(d['expr:'+gene])
  for adj in ['raw','depth','depth_cycle']:
   C=np.ones((len(d),1))
   if adj!='raw':C=np.column_stack([C,rankz(np.log1p(d.total_counts))])
   if adj=='depth_cycle':C=np.column_stack([C,rankz(cycle_score)])
   rho,_=corrp(A,Y,C);cors.append(dict(patient=pat,gene=gene,adjustment=adj,cells=len(d),rho=rho))
PC=pd.DataFrame(cors);PC.to_csv(O/'PART1_PATIENT_TYMP_GENE_CORRELATIONS.tsv',sep='\t',index=False)
pd.DataFrame(coverage).to_csv(O/'PART1_GENE_COVERAGE.tsv',sep='\t',index=False);PAT=pd.DataFrame(patient);PAT.to_csv(O/'PART1_PATIENT_TARGET_EXPRESSION.tsv',sep='\t',index=False)
ps=[]
for (gene,adj),d in PC.groupby(['gene','adjustment']):
 ps.append(dict(gene=gene,adjustment=adj,n_patients=len(d),positive=int((d.rho>0).sum()),negative=int((d.rho<0).sum()),median_rho=d.rho.median(),P=wilcoxon(d.rho).pvalue))
PS=pd.DataFrame(ps)
for adj,ix in PS.groupby('adjustment').groups.items():
 ia=[i for i in ix if PS.loc[i,'gene']!='SLC16A3'];PS.loc[ia,'FDR_axis6']=multipletests(PS.loc[ia,'P'],method='fdr_bh')[1]
PS.to_csv(O/'PART1_TYMP_GENE_CORRELATION_SUMMARY.tsv',sep='\t',index=False)
# Public paired MTA gene decomposition, using frozen Phase160 patients and covariate.
PUB=Path('/mnt/BioDisk/wangmingyuan/project26210_phase98_public_paired_metabolomics_v1')
mapping=pd.read_csv(PUB/'MasterMapping_ccRCC.csv',dtype=str).fillna('')
rows=[];interaction=[]
for cohort in ['RC18','RC20','CPTAC','CPTAC_val']:
 f=pd.read_csv(P/'phase160_external_bridge_investigation_20260915_v1'/('PATIENT_VALUES_'+cohort+'.tsv'),sep='\t',index_col=0)
 rr=pd.read_csv(PUB/'downloads_v1'/('matched_tpm_'+cohort+'.csv'),index_col=0)
 mm=mapping[(mapping.Dataset==cohort)&(mapping.TN=='Tumor')].copy()
 mm['patient']=mm.CommonID.where(mm.CommonID!='',mm.RNAID.str.replace(r'-T$','',regex=True)).str.replace(r'^RC(?:18|20):','',regex=True)
 ee=np.log2(rr.loc[mm.RNAID,axis]+1);ee.index=mm.patient.to_numpy();f=f.join(ee)
 z=(f[axis]-f[axis].mean())/f[axis].std(ddof=1)
 f['production_core']=z[['ODC1','SRM','SMS']].mean(1);f['SAM_supply']=z[['AMD1','MAT2A']].mean(1)
 f.to_csv(O/('PUBLIC_PATIENTS_'+cohort+'.tsv'),sep='\t')
 for gene in axis+['production_core','SAM_supply']:
  family='gene6' if gene in axis else 'module2'
  for adj in ['raw','myeloid']:
   cov=[] if adj=='raw' else ['MYELOID_PROXY']
   for outcome in ['MET:5-methylthioadenosine (MTA)','SPP1']:
    d=f[['TYMP',gene,outcome]+cov].dropna();n=len(d);C=np.column_stack([np.ones(n)]+[rankz(d[x]) for x in cov]);A=rankz(d.TYMP);M=rankz(d[gene]);Y=rankz(d[outcome])
    for edge,xx,yy,cc in [('a',A,M,C),('XY',M,Y,C),('b',M,Y,np.column_stack([C,A]))]:
     r,p=corrp(xx,yy,cc);rows.append(dict(cohort=cohort,gene=gene,family=family,adjustment=adj,outcome=outcome,edge=edge,n=n,k=np.linalg.matrix_rank(cc)-1,rho=r,P=p))
   if gene!='MTAP':
    d=f[['TYMP',gene,'MET:5-methylthioadenosine (MTA)','MTAP']+cov].dropna();C=np.column_stack([np.ones(len(d)),rankz(d.TYMP),rankz(d.MTAP)]+[rankz(d[x]) for x in cov])
    r,p=corrp(rankz(d[gene]),rankz(d['MET:5-methylthioadenosine (MTA)']),C)
    rows.append(dict(cohort=cohort,gene=gene,family=family,adjustment=adj,outcome='MET:5-methylthioadenosine (MTA)',edge='b_adjust_MTAP',n=len(d),k=np.linalg.matrix_rank(C)-1,rho=r,P=p))
 for adj in ['raw','myeloid']:
  cov=[] if adj=='raw' else ['MYELOID_PROXY'];d=f[['TYMP','MTAP','MET:5-methylthioadenosine (MTA)']+cov].dropna()
  a=rankz(d.TYMP);low=-rankz(d.MTAP);y=rankz(d['MET:5-methylthioadenosine (MTA)']);design=np.column_stack([np.ones(len(d)),a,low,a*low]+[rankz(d[x]) for x in cov])
  fit=sm.OLS(y,design).fit(cov_type='HC3',use_t=True)
  interaction.append(dict(cohort=cohort,adjustment=adj,n=len(d),beta=float(fit.params[3]),SE=float(fit.bse[3]),P=float(fit.pvalues[3])))
R=pd.DataFrame(rows)
for _,ix in R.groupby(['cohort','family','adjustment','outcome','edge']).groups.items():R.loc[ix,'FDR']=multipletests(R.loc[ix,'P'],method='fdr_bh')[1]
R.to_csv(O/'PUBLIC_COHORT_AXIS_ASSOCIATIONS.tsv',sep='\t',index=False)
def meta_corr(d):
 z=np.arctanh(d.rho.to_numpy());w=(d.n-d.k-3).to_numpy();mu=np.average(z,weights=w);q=np.sum(w*(z-mu)**2);df=len(d)-1;tau=max(0,(q-df)/(w.sum()-(w*w).sum()/w.sum()));wr=1/(1/w+tau);m=np.average(z,weights=wr);se=1/np.sqrt(wr.sum())
 return dict(cohorts=len(d),total_n=int(d.n.sum()),positive_cohorts=int((d.rho>0).sum()),random_rho=np.tanh(m),random_low=np.tanh(m-1.96*se),random_high=np.tanh(m+1.96*se),random_P=2*norm.sf(abs(m/se)),I2=max(0,(q-df)/q)*100 if q else 0)
mr=[]
for (gene,fam,adj,out,e),d in R.groupby(['gene','family','adjustment','outcome','edge']):mr.append(dict(gene=gene,family=fam,adjustment=adj,outcome=out,edge=e,**meta_corr(d)))
MR=pd.DataFrame(mr)
for _,ix in MR.groupby(['family','adjustment','outcome','edge']).groups.items():MR.loc[ix,'random_FDR']=multipletests(MR.loc[ix,'random_P'],method='fdr_bh')[1]
MR.to_csv(O/'PUBLIC_META_AXIS_ASSOCIATIONS.tsv',sep='\t',index=False)
IT=pd.DataFrame(interaction);IT.to_csv(O/'PUBLIC_TYMP_LOW_MTAP_INTERACTION.tsv',sep='\t',index=False)
imt=[]
for adj,d in IT.groupby('adjustment'):
 z=d.beta.to_numpy();v=d.SE.to_numpy()**2;w=1/v;mu=np.average(z,weights=w);q=np.sum(w*(z-mu)**2);tau=max(0,(q-len(d)+1)/(w.sum()-(w*w).sum()/w.sum()));wr=1/(v+tau);m=np.average(z,weights=wr);se=1/np.sqrt(wr.sum())
 imt.append(dict(adjustment=adj,cohorts=4,n=int(d.n.sum()),positive_cohorts=int((z>0).sum()),beta=m,low=m-1.96*se,high=m+1.96*se,P=2*norm.sf(abs(m/se))))
pd.DataFrame(imt).to_csv(O/'PUBLIC_META_MTAP_INTERACTION.tsv',sep='\t',index=False)
# Genetic distinction, no new FDR gate and no transfer of old endpoint to C4.
tc=pd.read_csv('/mnt/BioDisk/wangmingyuan/project26210_phase79_mtap_natural_perturbation_v1/results_v1/tcga_kirc_natural_perturbation_cohort_v1.tsv',sep='\t')
tcrows=[]
for val,d in tc.groupby('mtap_cna'):
 tcrows.append(dict(MTAP_GISTIC_call=val,patients=len(d),fraction=len(d)/len(tc),median_MTAP_expression_z=d.mtap_expr_z.median(),median_TYMP_z=d.tymp_z.median()))
TC=pd.DataFrame(tcrows);TC.to_csv(O/'TCGA_MTAP_SHALLOW_DEEP_DISTINCTION.tsv',sep='\t',index=False)
checks=[]
checks.append(dict(check='raw_cell_number',passed=len(meta)==11950));checks.append(dict(check='raw_TYMP_identity',passed=np.allclose(raw[:,target.index('TYMP')],meta.TYMP_raw_count_v1)))
for _,d in D.groupby(['contrast','gene']):checks.append(dict(check='frozen_pairs_'+d.contrast.iloc[0]+d.gene.iloc[0],passed=len(d)==14))
for i,r in PC[PC.adjustment=='raw'].iterrows():
 d=meta[meta.curated_patient_id==r.patient];rr=spearmanr(d['expr:TYMP'],d['expr:'+r.gene]).statistic
 checks.append(dict(check='scipy_'+str(i),passed=abs(rr-r.rho)<1e-10))
checks.append(dict(check='TCGA_deep_shallow_partition',passed=int((tc.mtap_cna==-2).sum())==6 and int((tc.mtap_cna==-1).sum())==142))
V=pd.DataFrame(checks);assert V.passed.all();V.to_csv(O/'VALIDATION.tsv',sep='\t',index=False)
status=dict(status='COMPUTATION_COMPLETE_REPORT_PENDING',part1_cells=11950,part1_patients=26,frozen_paired_patients=14,public_patients=255,TCGA_patients=len(tc),checks=len(V),passed=int(V.passed.sum()),no_figures=True,no_client_data_in_models=True,time=datetime.datetime.now(datetime.timezone.utc).isoformat())
(O/'STATUS.json').write_text(json.dumps(status,indent=2))
print('STATUS',json.dumps(status));print('PART1 CONTRAST\n'+DS.to_string(index=False));print('CORRELATIONS\n'+PS.to_string(index=False))
print('MTAP EXPRESSION\n'+PAT[(PAT.gene=='MTAP')&(PAT.cells>=50)].to_string(index=False))
print('PUBLIC MTA\n'+MR[(MR.outcome=='MET:5-methylthioadenosine (MTA)')&(MR.adjustment=='myeloid')][['gene','edge','positive_cohorts','random_rho','random_P','random_FDR','I2']].to_string(index=False))
print('PUBLIC INTERACTION\n'+pd.DataFrame(imt).to_string(index=False));print('TCGA\n'+TC.to_string(index=False))

