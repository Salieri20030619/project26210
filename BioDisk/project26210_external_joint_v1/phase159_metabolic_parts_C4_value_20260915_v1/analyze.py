from pathlib import Path
import pandas as pd,numpy as np,sys,json,subprocess,hashlib,datetime
from scipy.stats import spearmanr,wilcoxon,rankdata
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');S=P/'phase157_ordered_bridge_reassessment_20260915_v1';O=P/'phase159_metabolic_parts_C4_value_20260915_v1';O.mkdir(exist_ok=True)
gate=subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True);(O/'MEMORY_GATE.log').write_text(gate.stdout+gate.stderr);assert gate.returncode==0
sys.path.insert(0,str(S));from bridge_stats import bridge,rankz,corrp
blocks={
'cholesterol_mevalonate':['HMGCS1','HMGCR','MVK','PMVK','MVD','IDI1','FDPS'],
'cholesterol_sterol_synthesis':['FDFT1','SQLE','LSS','CYP51A1','MSMO1','NSDHL','HSD17B7','EBP','SC5D','DHCR7','DHCR24'],
'cholesterol_oxysterol_generation':['CYP27A1','CH25H'],
'cholesterol_oxysterol_clearance':['CYP7B1','HSD3B7'],
'cholesterol_efflux_carrier':['ABCA1','ABCG1','APOE'],
'glycolysis_glucose_entry':['SLC2A1','HK2'],
'glycolysis_core':['PFKP','ALDOA','GAPDH','PGK1','ENO1','PKM'],
'glycolysis_lactate_conversion':['LDHA','LDHB'],
'glycolysis_monocarboxylate_transport':['SLC16A3','SLC16A1','BSG']}
receivers={'cholesterol':['LDLR','SCARB1','NR1H3','NR1H2','RXRA','ABCA1','ABCG1','CYP27A1'],'glycolysis':['SLC16A1','SLC16A3','BSG','HCAR1','GPR132','LDHA','LDHB']}
sources={'cholesterol_synthesis':'https://reactome.org/content/detail/R-HSA-191273','CYP27A1':'https://reactome.org/content/detail/R-HSA-193807','MCT4':'https://www.uniprot.org/uniprotkb/O15427/entry','MCT1':'https://www.uniprot.org/uniprotkb/P53985/entry','GPR132':'https://pmc.ncbi.nlm.nih.gov/articles/PMC5255630/'}
plan=dict(blocks=blocks,receiver_genes=receivers,sources=sources,primary_outcomes=['C4_fraction','SPP1_pseudobulk','SPP1_detection_fraction'],QC='sender gene >=10 pseudobulk counts in >=12/16 libraries; modules >=2 and >=50% genes eligible; recipient expression described without censoring by association',statistics='patient units; families pathway x gene/module x endpoint x edge x adjustment; selected secondary cohort comparison; rank partials t approximation; sample paired differences; stratified paired patient bootstrap rho-difference reranks each resample',limits='post-selection decomposition, no causal mediator claim, no metabolite concentration or flux claims, LDHA/LDHB and MCT1/MCT4 separate directions in interpretation; no figures; no client data')
(O/'PLAN_BEFORE_FIT.json').write_text(json.dumps(plan,indent=2))
B=pd.read_csv(P/'phase154_same_library_TYMP_C4_test_20260915_v1/SAME_LIBRARY_PATIENT_METRICS.tsv',sep='\t').set_index('patient')
G=pd.read_csv(S/'02_inputs/GROUP_INDEX.tsv',sep='\t')
Z=np.load(S/'02_inputs/GROUP_PSEUDOBULK.npz',allow_pickle=True);counts=Z['counts'];det=Z['detected'];genes=list(Z['genes']);lookup={g:i for i,g in enumerate(genes)}
GI=G.set_index('group_index');assert list(G.group_index)==list(range(len(G)))
malids=[]
for pat,b in B.iterrows():
 g=G[(G.patient==pat)&(G.sample_id==b.sample_id)&(G.role=='Malignant')]
 assert len(g)==1;malids.append(int(g.iloc[0].group_index))
MC=counts[malids];MX=pd.DataFrame(np.log2(MC/MC.sum(1)[:,None]*1e6+1),index=B.index,columns=genes)
values={};reg=[];coverage=[]
for block,gg in blocks.items():
 keep=[]
 for gene in gg:
  idx=lookup.get(gene);covered=0 if idx is None else int((MC[:,idx]>=10).sum())
  coverage.append(dict(block=block,gene=gene,pseudobulk_libraries_ge10=covered,eligible=covered>=12))
  if covered<12:continue
  keep.append(gene);name='gene:'+gene;values[name]=MX[gene]
  reg.append(dict(candidate=name,family=block.split('_')[0],level='gene',block=block,genes=gene))
 if len(keep)>=2 and len(keep)>=len(gg)*.5:
  z=MX[keep];z=(z-z.mean())/z.std(ddof=1)
  name='module:'+block;values[name]=z.mean(axis=1)
  reg.append(dict(candidate=name,family=block.split('_')[0],level='module',block=block,genes=','.join(keep)))
X=pd.DataFrame(values,index=B.index);REG=pd.DataFrame(reg).drop_duplicates('candidate')
pd.DataFrame(coverage).to_csv(O/'SENDER_COVERAGE_ALL_GENES.tsv',sep='\t',index=False);REG.to_csv(O/'CANDIDATE_REGISTRY.tsv',sep='\t',index=False);X.to_csv(O/'PATIENT_SENDER_VALUES.tsv',sep='\t')
rows=[]
for _,r in REG.iterrows():
 for endpoint in plan['primary_outcomes']:
  for adj in ['raw','study']:
   rows.append(dict(**r.to_dict(),endpoint=endpoint,adjustment=adj,**bridge(B.TYMP_pseudobulk,X[r.candidate],B[endpoint],B.study,adj)))
R=pd.DataFrame(rows)
for _,idx in R.groupby(['family','level','endpoint','adjustment']).groups.items():
 for edge in ['a','XY','b']:R.loc[idx,edge+'_FDR']=multipletests(R.loc[idx,edge+'_P'],method='fdr_bh')[1]
R.to_csv(O/'ALL_COMPONENT_BRIDGE_RESULTS.tsv',sep='\t',index=False)
# Same-library endpoint audit with exact raw RNA and detection contributions.
spp=lookup['SPP1'];crows=[]
for pat,b in B.iterrows():
 q=G[(G.patient==pat)&(G.sample_id==b.sample_id)&G.role.isin(['C4','OtherMacrophage'])]
 ix=q.group_index.to_numpy(int); c=q[q.role=='C4'].group_index.to_numpy(int);o=q[q.role=='OtherMacrophage'].group_index.to_numpy(int)
 assert int(q.cells.sum())==b.macrophage_cells
 assert abs(counts[ix,spp].sum()-b.macrophage_SPP1_counts)<1e-7
 ccounts=counts[c].sum(0);ocounts=counts[o].sum(0);nc=int(q[q.role=='C4'].cells.sum());no=int(q[q.role=='OtherMacrophage'].cells.sum())
 crows.append(dict(patient=pat,C4_cells=nc,Other_cells=no,C4_fraction=nc/(nc+no),C4_SPP1_RNA_share=ccounts[spp]/(ccounts[spp]+ocounts[spp]),C4_SPP1_positive_cell_share=det[c,spp].sum()/det[ix,spp].sum(),C4_SPP1_log2CPM=np.log2(ccounts[spp]/ccounts.sum()*1e6+1) if nc else np.nan,Other_SPP1_log2CPM=np.log2(ocounts[spp]/ocounts.sum()*1e6+1),C4_SPP1_detection=det[c,spp].sum()/nc if nc else np.nan,Other_SPP1_detection=det[o,spp].sum()/no))
D=pd.DataFrame(crows).set_index('patient').loc[B.index];D.to_csv(O/'SAME_LIBRARY_C4_CONTRIBUTION.tsv',sep='\t')
# Tumor patient aggregation of C4 and Other, >=20 each for paired state biology.
panel=sorted(set([g for gg in blocks.values() for g in gg]+[g for gg in receivers.values() for g in gg]+['SPP1','GPNMB','CCL18']))
paired=[]
for pat,gg in G[G.source_tissue.str.lower()=='tumor'].groupby('patient'):
 c=gg[gg.role=='C4'];o=gg[gg.role=='OtherMacrophage']
 nc=int(c.cells.sum());no=int(o.cells.sum())
 if nc<20 or no<20:continue
 ci=c.group_index.to_numpy(int);oi=o.group_index.to_numpy(int);cc=counts[ci].sum(0);oc=counts[oi].sum(0)
 for gene in panel:
  j=lookup.get(gene)
  if j is None:continue
  cexp=np.log2(cc[j]/cc.sum()*1e6+1);oexp=np.log2(oc[j]/oc.sum()*1e6+1)
  family='marker' if gene in ['SPP1','GPNMB','CCL18'] else ('cholesterol' if any(gene in v for k,v in blocks.items() if k.startswith('cholesterol')) or gene in receivers['cholesterol'] else 'glycolysis')
  paired.append(dict(patient=pat,gene=gene,family=family,C4_cells=nc,Other_cells=no,C4_log2CPM=cexp,Other_log2CPM=oexp,difference=cexp-oexp,C4_detection=det[ci,j].sum()/nc,Other_detection=det[oi,j].sum()/no))
PD=pd.DataFrame(paired);assert len(PD)>0, G.source_tissue.unique()
PD.to_csv(O/'C4_OTHER_PAIRED_GENE_PATIENTS.tsv',sep='\t',index=False)
pr=[]
for (gene,family),g in PD.groupby(['gene','family']):
 d=g.difference.to_numpy();p=float(wilcoxon(d).pvalue) if np.any(d) else 1.
 pr.append(dict(gene=gene,family=family,n_patients=len(g),positive=int((d>0).sum()),negative=int((d<0).sum()),median_difference=np.median(d),P=p,C4_detection_median=g.C4_detection.median(),C4_receiver_ge10pct=int((g.C4_detection>=.1).sum()),Other_detection_median=g.Other_detection.median()))
PR=pd.DataFrame(pr)
for _,idx in PR.groupby('family').groups.items():PR.loc[idx,'FDR']=multipletests(PR.loc[idx,'P'],method='fdr_bh')[1]
PR.to_csv(O/'C4_OTHER_GENE_SUMMARY.tsv',sep='\t',index=False)
# Fair matched endpoint comparison, with full and common 12-patient subsets.
oldX=pd.read_csv(S/'04_other_X_screen/PATIENT_CANDIDATE_VALUES.tsv',sep='\t',index_col=0).loc[B.index]
predictors={'TYMP':B.TYMP_pseudobulk,'IGFBP3':oldX['ligand:IGFBP3'],'Cholesterol Homeostasis':oldX['program:Cholesterol Homeostasis'],'Glycolysis':oldX['program:Glycolysis']}
audit=[];differences=[];rng=np.random.default_rng(262159)
for scope,idx in [('all16',B.index),('C4_ge20_common',D.index[D.C4_cells>=20])]:
 ys={'C4_fraction':B.loc[idx,'C4_fraction'],'SPP1_all':B.loc[idx,'SPP1_pseudobulk'],'SPP1_detection_all':B.loc[idx,'SPP1_detection_fraction']}
 if scope!='all16':ys.update(SPP1_C4=D.loc[idx,'C4_SPP1_log2CPM'],SPP1_other=D.loc[idx,'Other_SPP1_log2CPM'])
 studies=B.loc[idx,'study'];ng=len(idx);Cstudy=pd.get_dummies(studies,dtype=float).to_numpy()
 for name,pred in predictors.items():
  xx=pred.loc[idx]
  for endpoint,yy in ys.items():
   for adj in ['raw','study']:
    C=np.ones((ng,1)) if adj=='raw' else Cstudy
    rho,p=corrp(rankz(xx),rankz(yy),C)
    audit.append(dict(scope=scope,predictor=name,endpoint=endpoint,adjustment=adj,n_patients=ng,rho=rho,P=p))
  contrasts=[('C4_fraction','SPP1_all'),('C4_fraction','SPP1_detection_all')]
  if scope!='all16':contrasts+=[('SPP1_C4','SPP1_all'),('SPP1_C4','SPP1_other')]
  groups=[np.flatnonzero(studies.to_numpy()==s) for s in studies.unique()]
  ii=np.tile(np.arange(ng),(2000,1))
  for ids in groups:ii[:,ids]=rng.choice(ids,(2000,len(ids)),replace=True)
  def bcorr(a,b):
   aa=rankdata(np.asarray(a)[ii],axis=1);bb=rankdata(np.asarray(b)[ii],axis=1)
   aa-=aa.mean(axis=1)[:,None];bb-=bb.mean(axis=1)[:,None]
   denom=np.sqrt((aa*aa).sum(1)*(bb*bb).sum(1))
   return np.divide((aa*bb).sum(1),denom,out=np.full(2000,np.nan),where=denom>0)
  for e1,e2 in contrasts:
   delta=bcorr(xx,ys[e1])-bcorr(xx,ys[e2]);valid=delta[np.isfinite(delta)]
   differences.append(dict(scope=scope,predictor=name,endpoint1=e1,endpoint2=e2,n=ng,rho_difference=spearmanr(xx,ys[e1]).statistic-spearmanr(xx,ys[e2]).statistic,bootstrap_low=np.quantile(valid,.025),bootstrap_high=np.quantile(valid,.975),valid=len(valid)))
A=pd.DataFrame(audit)
for _,ix in A.groupby(['scope','predictor','adjustment']).groups.items():A.loc[ix,'FDR_endpoint_audit']=multipletests(A.loc[ix,'P'],method='fdr_bh')[1]
A.to_csv(O/'FAIR_ENDPOINT_COMPARISON.tsv',sep='\t',index=False);pd.DataFrame(differences).to_csv(O/'PAIRED_CORRELATION_DIFFERENCE_BOOTSTRAP.tsv',sep='\t',index=False)
# Component stability for all, avoiding choosing a subset by P.
stab=[]
for _,r in REG.iterrows():
 for endpoint in plan['primary_outcomes']:
  for kind in ['patient','study']:
   groups=B.index if kind=='patient' else B.study.unique()
   rr=[]
   for omit in groups:
    mask=B.index!=omit if kind=='patient' else B.study!=omit
    fit=bridge(B.loc[mask,'TYMP_pseudobulk'],X.loc[mask,r.candidate],B.loc[mask,endpoint],B.loc[mask,'study'],'raw')
    rr.append(dict(candidate=r.candidate,endpoint=endpoint,kind=kind,omitted=omit,**fit))
   stab.extend(rr)
ST=pd.DataFrame(stab);ST.to_csv(O/'ALL_COMPONENT_LEAVE_ONE.tsv',sep='\t',index=False)
v=[]
for i,r in R.iterrows():
 v.append(dict(check='identity_'+str(i),passed=abs(r.identity_error)<1e-10))
 if r.adjustment=='raw':
  rho,p=spearmanr(X[r.candidate],B[r.endpoint])
  v.append(dict(check='independent_scipy_'+str(i),passed=abs(rho-r.XY_rho)<1e-10 and abs(p-r.XY_P)<1e-10))
v.append(dict(check='C4_counts_identity',passed=np.allclose(D.C4_fraction,B.C4_fraction)))
v.append(dict(check='C4_share_bounded',passed=bool(D.C4_SPP1_RNA_share.between(0,1).all())))
V=pd.DataFrame(v);V.to_csv(O/'VALIDATION.tsv',sep='\t',index=False);assert V.passed.all()
status=dict(status='COMPUTATION_COMPLETE_REPORT_PENDING',candidates=len(REG),genes=int((REG.level=='gene').sum()),modules=int((REG.level=='module').sum()),paired_state_patients=int(PR.n_patients.max()),matched_C4_expression_patients=int((D.C4_cells>=20).sum()),checks=len(V),passed=int(V.passed.sum()),time=datetime.datetime.now(datetime.timezone.utc).isoformat(),figures_created=False,client_data_used=False)
(O/'STATUS.json').write_text(json.dumps(status,indent=2))
files=[S/'02_inputs/GROUP_PSEUDOBULK.npz',S/'02_inputs/GROUP_INDEX.tsv',P/'phase154_same_library_TYMP_C4_test_20260915_v1/SAME_LIBRARY_PATIENT_METRICS.tsv',S/'04_other_X_screen/PATIENT_CANDIDATE_VALUES.tsv',S/'bridge_stats.py']
(O/'SOURCE_MANIFEST.json').write_text(json.dumps([dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in files],indent=2))
print(json.dumps(status))
print('RAW BOTH EDGES POSITIVE P<.1')
q=R[(R.adjustment=='raw')&(R.a_rho>0)&(R.b_partial_rho>0)&(R.a_P<.1)&(R.b_P<.1)]
print(q[['candidate','endpoint','a_rho','a_P','a_FDR','XY_rho','XY_P','b_partial_rho','b_P','b_FDR']].to_string(index=False))
print('ALL MODULES')
print(R[R.level=='module'][['candidate','endpoint','adjustment','a_rho','a_P','XY_rho','XY_P','b_partial_rho','b_P']].to_string(index=False))
print('STATE')
print(PR[PR.gene.isin(['SPP1','GPNMB','CCL18','CYP27A1','CH25H','ABCA1','ABCG1','LDLR','SCARB1','NR1H3','NR1H2','RXRA','SLC16A1','SLC16A3','BSG','HCAR1','GPR132','LDHA','LDHB'])].to_string(index=False))
print('CONTRIBUTION MEDIAN',D[['C4_fraction','C4_SPP1_RNA_share','C4_SPP1_positive_cell_share']].median().to_dict())
print('FAIR RAW\n'+A[A.adjustment=='raw'].to_string(index=False))
print('DIFF\n'+pd.DataFrame(differences).to_string(index=False))

