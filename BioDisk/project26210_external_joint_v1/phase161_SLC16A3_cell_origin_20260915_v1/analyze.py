from pathlib import Path
import pandas as pd,numpy as np,h5py,json,sys,subprocess,hashlib,datetime
from scipy.sparse import csr_matrix
from scipy.stats import spearmanr,wilcoxon
from anndata.io import read_elem
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');S=P/'phase157_ordered_bridge_reassessment_20260915_v1';O=P/'phase161_SLC16A3_cell_origin_20260915_v1';O.mkdir(exist_ok=True)
g=subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True);(O/'MEMORY_GATE.log').write_text(g.stdout+g.stderr);assert g.returncode==0
sys.path.insert(0,str(S));from bridge_stats import rankz,corrp
B=pd.read_csv(P/'phase154_same_library_TYMP_C4_test_20260915_v1/SAME_LIBRARY_PATIENT_METRICS.tsv',sep='\t').set_index('patient')
C=pd.read_csv(S/'02_inputs/CELL_IDENTITY.tsv.gz',sep='\t')
C=C[C.curated_patient_id.isin(B.index)].copy()
C=C[C.apply(lambda r:r.sample_id==B.loc[r.curated_patient_id,'sample_id'],axis=1)].reset_index(drop=True)
assert len(C)==23018
targets=['TYMP','SLC16A3','LDHA','MSMO1','SPP1','CYP27A1','BSG','SLC16A1']
edges={'Malignant':[('TYMP','SLC16A3'),('TYMP','MSMO1')],'Macrophage':[('SLC16A3','SPP1'),('LDHA','SPP1')],'C4':[('SLC16A3','SPP1'),('LDHA','SPP1')],'OtherMacrophage':[('SLC16A3','SPP1'),('LDHA','SPP1')]}
plan=dict(matched_patients=16,cells_expected=23018,targets=targets,edges=edges,minimum_cells=50,minimum_detected_each_gene=10,primary_adjustment='depth',sensitivity=['raw','depth_hypoxia'],unit='Independent patient; within-patient cell correlations are descriptive effects, no cell-level inferential P exported; signed rank test of patient rho versus zero',C4='frozen members; secondary subset, same QC',source_scope='Only malignant and mature macrophage compartments available here; not an exhaustive whole-tissue origin comparison',hypoxia='Hallmark gene z-expression mean within each patient and compartment; exclude all target genes',no_figures=True,no_client_data_in_models=True)
(O/'PLAN_BEFORE_FIT.json').write_text(json.dumps(plan,indent=2))
gmt=P/'part2_rebuild_v1/06_functional_programs/run_v1/resources/MSigDB_Hallmark_2020.gmt'
hyp=[]
for line in gmt.read_text().splitlines():
 parts=line.split('\t')
 if parts[0]=='Hypoxia':hyp=[x for x in parts[2:] if x not in targets]
assert len(hyp)>10
raw='/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/03_malignancy_cnv/cnv_annotated_h5ad_v1/part1_scvi_integrated_cell_annotated_cnv_v1.h5ad'
with h5py.File(raw,'r') as f:
 ids=np.asarray(read_elem(f['obs/cell_id'])).astype(str);genes=list(np.asarray(read_elem(f['var/gene'])).astype(str))
 ix=pd.Index(ids).get_indexer(C.cell_id);assert (ix>=0).all() and len(set(ix))==len(ix)
 ptr=f['layers/counts/indptr'][:];idx=f['layers/counts/indices'][:];dat=f['layers/counts/data'][:]
 mat=csr_matrix((dat,idx,ptr),shape=(len(ids),len(genes)),copy=False)
 selected=mat[ix];del mat,dat,idx,ptr
 totals=np.asarray(selected.sum(1)).ravel();assert (totals>0).all()
 hs=[x for x in hyp if x in genes];gg=targets+hs;rawsmall=selected[:,[genes.index(x) for x in gg]].toarray();del selected
E=np.log1p(rawsmall/totals[:,None]*1e4);expr=pd.DataFrame(E,columns=gg)
C['total_counts']=totals
for gene in targets:C['counts:'+gene]=rawsmall[:,gg.index(gene)];C['expr:'+gene]=expr[gene].to_numpy()
C.to_csv(O/'CELL_TARGET_EXPRESSION.tsv.gz',sep='\t',index=False,compression='gzip')
np.savez_compressed(O/'HYPoxia_CELL_INPUT.npz',expression=E[:,len(targets):],genes=np.array(hs),cell_id=C.cell_id.to_numpy())
rows=[];qc=[];source=[]
for pat in B.index:
 for compartment in ['Malignant','Macrophage']:
  mask=(C.curated_patient_id==pat)&(C.role=='Malignant' if compartment=='Malignant' else C.role.isin(['C4','OtherMacrophage']))
  d=C[mask]
  for gene in targets:
   source.append(dict(patient=pat,compartment=compartment,gene=gene,cells=len(d),counts=float(d['counts:'+gene].sum()),detection=float((d['counts:'+gene]>0).mean()),log2CPM=float(np.log2(d['counts:'+gene].sum()/d.total_counts.sum()*1e6+1))))
for compartment,pairs in edges.items():
 for pat in B.index:
  mask=(C.curated_patient_id==pat)&(C.role.isin(['C4','OtherMacrophage']) if compartment=='Macrophage' else C.role==compartment)
  d=C[mask];n=len(d)
  for x,y in pairs:
   ndx=int((d['counts:'+x]>0).sum());ndy=int((d['counts:'+y]>0).sum());ok=n>=50 and ndx>=10 and ndy>=10
   qc.append(dict(compartment=compartment,patient=pat,x=x,y=y,cells=n,detected_x=ndx,detected_y=ndy,eligible=ok))
   if not ok:continue
   H=E[d.index,len(targets):];sd=H.std(0,ddof=1);valid=sd>0;h=((H[:,valid]-H[:,valid].mean(0))/sd[valid]).mean(1)
   A=rankz(d['expr:'+x]);Y=rankz(d['expr:'+y])
   for adj in ['raw','depth','depth_hypoxia']:
    cov=np.ones((n,1))
    if adj!='raw':cov=np.column_stack([cov,rankz(np.log1p(d.total_counts))])
    if adj=='depth_hypoxia':cov=np.column_stack([cov,rankz(h)])
    rho,_,_= (None,None,None) if False else (0,0,0)
    rho,_=corrp(A,Y,cov)
    rows.append(dict(compartment=compartment,patient=pat,study=B.loc[pat,'study'],x=x,y=y,adjustment=adj,cells=n,rho=rho))
R=pd.DataFrame(rows);R.to_csv(O/'PATIENT_WITHIN_CELL_CORRELATIONS.tsv',sep='\t',index=False)
pd.DataFrame(qc).to_csv(O/'PATIENT_COVERAGE_QC.tsv',sep='\t',index=False)
source=pd.DataFrame(source);source.to_csv(O/'PAIRED_COMPARTMENT_SOURCE.tsv',sep='\t',index=False)
ss=[]
for (comp,x,y,adj),d in R.groupby(['compartment','x','y','adjustment']):
 rho=d.rho.to_numpy();p=float(wilcoxon(rho).pvalue) if len(rho)>=3 else np.nan
 ss.append(dict(compartment=comp,x=x,y=y,adjustment=adj,n_patients=len(rho),positive=int((rho>0).sum()),negative=int((rho<0).sum()),median_rho=np.median(rho),P=p,exploratory=len(rho)<10))
SS=pd.DataFrame(ss)
for _,ix in SS.groupby(['compartment','adjustment']).groups.items():
 vals=SS.loc[ix,'P'];ok=vals.notna();SS.loc[vals.index[ok],'FDR']=multipletests(vals[ok],method='fdr_bh')[1]
SS.to_csv(O/'PATIENT_LEVEL_CORRELATION_SUMMARY.tsv',sep='\t',index=False)
sr=[]
for gene,d in source.groupby('gene'):
 wide=d.pivot(index='patient',columns='compartment',values='log2CPM');delta=wide.Malignant-wide.Macrophage
 det=d.pivot(index='patient',columns='compartment',values='detection')
 sr.append(dict(gene=gene,n_patients=len(delta),malignant_higher=int((delta>0).sum()),macrophage_higher=int((delta<0).sum()),median_malignant_minus_macro_log2CPM=delta.median(),P=wilcoxon(delta).pvalue,malignant_detection_median=det.Malignant.median(),macro_detection_median=det.Macrophage.median()))
SR=pd.DataFrame(sr);SR['FDR']=multipletests(SR.P,method='fdr_bh')[1];SR.to_csv(O/'COMPARTMENT_SOURCE_SUMMARY.tsv',sep='\t',index=False)
# Same-patient across-compartment coordination is kept separate from within-cell effects.
cross=[];wide=source.pivot(index='patient',columns=['compartment','gene'],values='log2CPM').loc[B.index]
for xname,x,yname,y in [('malignant_MCT4',wide[('Malignant','SLC16A3')],'macro_MCT4',wide[('Macrophage','SLC16A3')]),('malignant_TYMP',wide[('Malignant','TYMP')],'malignant_MCT4',wide[('Malignant','SLC16A3')]),('malignant_MCT4',wide[('Malignant','SLC16A3')],'macro_SPP1',wide[('Macrophage','SPP1')])]:
 for adj in ['raw','study']:
  cov=np.ones((16,1)) if adj=='raw' else pd.get_dummies(B.study,dtype=float).to_numpy()
  rho,p=corrp(rankz(x),rankz(y),cov);cross.append(dict(x=xname,y=yname,adjustment=adj,n_patients=16,rho=rho,P=p))
CR=pd.DataFrame(cross)
for _,ix in CR.groupby('adjustment').groups.items():CR.loc[ix,'FDR_current_coordination_family']=multipletests(CR.loc[ix,'P'],method='fdr_bh')[1]
CR.to_csv(O/'CROSS_COMPARTMENT_COORDINATION.tsv',sep='\t',index=False)
checks=[]
for pat in B.index:
 b=B.loc[pat];m=source[(source.patient==pat)&(source.compartment=='Malignant')&(source.gene=='TYMP')].iloc[0];s=source[(source.patient==pat)&(source.compartment=='Macrophage')&(source.gene=='SPP1')].iloc[0]
 checks.append(dict(check=pat+'_raw_TYMP',passed=abs(m.log2CPM-b.TYMP_pseudobulk)<1e-6))
 checks.append(dict(check=pat+'_raw_SPP1',passed=abs(s.log2CPM-b.SPP1_pseudobulk)<1e-6))
for i,r in R[R.adjustment=='raw'].iterrows():
 mask=(C.curated_patient_id==r.patient)&(C.role.isin(['C4','OtherMacrophage']) if r.compartment=='Macrophage' else C.role==r.compartment)
 rho=spearmanr(C.loc[mask,'expr:'+r.x],C.loc[mask,'expr:'+r.y]).statistic
 checks.append(dict(check='scipy_'+str(i),passed=abs(rho-r.rho)<1e-10))
V=pd.DataFrame(checks);assert V.passed.all();V.to_csv(O/'VALIDATION.tsv',sep='\t',index=False)
status=dict(status='COMPUTATION_COMPLETE_REPORT_PENDING',matched_patients=16,cells=len(C),checks=len(V),passed=int(V.passed.sum()),figures=False,client_data_in_models=False,time=datetime.datetime.now(datetime.timezone.utc).isoformat())
(O/'STATUS.json').write_text(json.dumps(status,indent=2))
(O/'SOURCE_MANIFEST.json').write_text(json.dumps(dict(raw_h5ad=raw,raw_gene_index=genes.index('SLC16A3'),hypoxia_gene_count=len(hs),files=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in [S/'02_inputs/CELL_IDENTITY.tsv.gz',P/'phase154_same_library_TYMP_C4_test_20260915_v1/SAME_LIBRARY_PATIENT_METRICS.tsv',gmt]]),indent=2))
print('STATUS',json.dumps(status));print('WITHIN\n'+SS.to_string(index=False));print('SOURCE\n'+SR.to_string(index=False));print('CROSS\n'+CR.to_string(index=False))

