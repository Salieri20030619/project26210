from pathlib import Path
import json,hashlib,math,re,subprocess
import pandas as pd,numpy as np
P=Path(__file__).resolve().parents[1];R=P.parent
q=subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True)
(P/'logs/terminal_evidence_memory.log').write_text(q.stdout+q.stderr);assert q.returncode==0
checks=[]
def ck(n,v):checks.append(dict(check=n,passed=bool(v)));assert v,n
e=pd.read_csv(P/'tables/PERTURBATION_ALL_GENES.tsv.gz',sep='\t',float_precision='round_trip')
sm=pd.read_csv(P/'tables/PERTURBATION_SAMPLE_MAP.tsv',sep='\t')
e['measurement_status']='evaluable_descriptive'
e['total_input_count']=np.nan
errors=[]
for file,z in e.groupby('file'):
 c=np.load(P/'inputs'/(Path(file).stem+'.npz'),allow_pickle=True)
 C=c['counts'];cols=c['columns'];genes=c['ensembl']
 design=sm.set_index('column').loc[cols]
 ix=design.RIG.to_numpy().astype(bool)
 ck('mapping_unique_'+file,design.GSM.is_unique)
 total=C.sum(0);uq=np.array([np.quantile(v[v>0],.75) for v in C.T]);uq=uq/np.exp(np.log(uq).mean())*np.exp(np.log(total).mean())
 for norm,den in [('CPM',total),('upper_quartile',uq)]:
  sel=z[z.normalization==norm].set_index('ENSEMBL').loc[genes]
  # Recompute entire output using column-wise normalization and mean differences.
  logs=np.stack([np.log2(C[:,j]*(1e6/den[j])+1) for j in range(C.shape[1])],axis=1)
  ef=logs[:,ix].sum(1)/ix.sum()-logs[:,~ix].sum(1)/(~ix).sum()
  err=float(np.max(np.abs(ef-sel.effect.to_numpy())));errors.append(err)
  ck('all_effects_'+file+'_'+norm,err<1e-12)
  globalidx=z[z.normalization==norm].index
  countmap=dict(zip(genes,C.sum(1)))
  e.loc[globalidx,'total_input_count']=e.loc[globalidx,'ENSEMBL'].map(countmap)
  e.loc[globalidx,'measurement_status']=np.where(e.loc[globalidx,'total_input_count'].eq(0),'all_zero_not_evaluable','evaluable_descriptive')
e.to_csv(P/'tables/PERTURBATION_ALL_GENES_QUALIFIED.tsv.gz',sep='\t',index=False)
targets=['TYMP','SPP1','CD44','ICAM1','VCAM1','CXCL9','CXCL10','CXCL11','CCL5']
e[e.SYMBOL.isin(targets)].to_csv(P/'tables/PERTURBATION_TARGETS_QUALIFIED.tsv',sep='\t',index=False)
# Recompute all program scores from complete counts and one-to-one mapping.
mp=pd.read_csv(R/'phase280_TYMP_source_niche_function_20260928_v1/tables/EXTERNAL_ONE_TO_ONE_MAPPING.tsv',sep='\t').set_index('ENSEMBL').SYMBOL.to_dict()
panels=json.loads((R/'phase292_barrier_execution_and_external_feasibility_20260930_v1/PANELS_FROZEN.json').read_text())
scores=pd.read_csv(P/'tables/PERTURBATION_PROGRAM_SCORES.tsv',sep='\t')
for file,z in scores.groupby('file'):
 c=np.load(P/'inputs'/(Path(file).stem+'.npz'),allow_pickle=True);C=c['counts'];genes=c['ensembl'];cols=list(c['columns']);total=C.sum(0)
 for row in z.itertuples():
  idx=[i for i,g in enumerate(genes) if mp.get(g) in panels[row.program]];j=cols.index(row.column)
  v=math.fsum(math.log1p(C[i,j]/total[j]*1e4) for i in idx)/len(idx)
  ck('program_'+str(len(checks)),abs(v-row.score)<1e-12)
# GSE264586 sample-level GEO identity alignment: normalized title, exact 1:1.
geo=pd.read_csv(P/'tables/GSE264586_GEO_SAMPLE_IDENTITY.tsv',sep='\t').fillna('')
norm=lambda s:re.sub('[^a-z0-9]','',s.lower().replace('cycle2','c2').replace('end_of_treatment','eot'))
geo['key']=geo.title.map(norm)
lib=pd.read_csv(P/'tables/EXTERNAL_LIBRARY_INPUTS.tsv',sep='\t')
samples=lib[lib.dataset=='GSE264586'].library.drop_duplicates()
align=[]
for sample in samples:
 hits=geo[geo.key.eq(norm(sample))]
 ck('geo_sample_'+sample,len(hits)==1)
 r=hits.iloc[0]
 align.append(dict(library=sample,GSM=r.GSM,title=r.title,tissue=r.source_name_ch1,characteristics=r.characteristics_ch1,ccRCC_reported='clear cell renal cell carcinoma' in r.characteristics_ch1))
pd.DataFrame(align).to_csv(P/'tables/GSE264586_SAMPLE_TISSUE_ALIGNMENT.tsv',sep='\t',index=False)
# Freeze archive-derived 10x panel, without assuming newer donor equivalence.
panel=json.loads((P/'inputs/10x_gene_panel.json').read_text())
names=[x['type']['data']['name'] for x in panel['payload']['targets'] if x['type']['descriptor']=='gene']
exp=json.loads((P/'inputs/10x_experiment.xenium').read_text())
ck('10x_405_unique_genes',len(set(names))==405)
ck('10x_old_465534_cells',exp['num_cells']==465534)
rows=[dict(gene=g,present=g in names) for g in ['TYMP','SPP1','CD68','CD163','DCN','LUM','COL1A1','COL3A1','FAP','POSTN','CD3D','CD3E','CD8A','CD8B','EPCAM','PAX8','CA9']]
pd.DataFrame(rows).to_csv(P/'tables/TENX_OLD_RNA_PANEL_QUALIFICATION.tsv',sep='\t',index=False)
pd.DataFrame(dict(gene=sorted(set(names)))).to_csv(P/'resources/TENX_OLD_FULL_GENE_LIST.tsv',sep='\t',index=False)
ck('10x_TYMP_SPP1_unmeasured','TYMP' not in names and 'SPP1' not in names)
tree=json.loads((P/'resources/nicheverse_tree.json').read_text());f=P/'inputs/xenium_rcc_core.h5ad'
entry=[x for x in tree['tree'] if x['path']=='examples/data/xenium_rcc_core.h5ad'][0]
ck('nicheverse_official_blob',hashlib.sha1(b'blob '+str(f.stat().st_size).encode()+b'\0'+f.read_bytes()).hexdigest()==entry['sha'])
pd.DataFrame(checks).to_csv(P/'qc/TERMINAL_EVIDENCE_CHECKS.tsv',sep='\t',index=False)
(P/'qc/TERMINAL_EVIDENCE_VALIDATION.json').write_text(json.dumps(dict(status='PASS',checks=len(checks),all_gene_effect_max_error=max(errors),aligned_GSE264586_samples=len(align),perturbation_all_zero_rows=int(e.measurement_status.eq('all_zero_not_evaluable').sum()),no_inferential_P=True),indent=2))
print(json.dumps(dict(checks=len(checks),effect_max_error=max(errors),tissues=pd.DataFrame(align).tissue.value_counts().to_dict())))
print(e[(e.normalization=='CPM')&e.SYMBOL.eq('SPP1')][['cell_line','TPI','effect','total_input_count','measurement_status']].to_string(index=False))
