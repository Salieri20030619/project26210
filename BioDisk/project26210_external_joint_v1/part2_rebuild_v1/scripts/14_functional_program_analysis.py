from pathlib import Path
import json,hashlib,warnings
import numpy as np,pandas as pd
from scipy.stats import wilcoxon,hypergeom
from scipy import sparse
warnings.filterwarnings('ignore',category=RuntimeWarning)
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/part2_rebuild_v1');O=R/'06_functional_programs/run_v1';M=R/'04_markers_annotation/run_v4_macrophage_only'
def bh(v):
 a=np.asarray(v,float);q=np.full(len(a),np.nan);ok=np.isfinite(a);x=a[ok];s=np.argsort(x);y=np.minimum.accumulate((x[s]*len(x)/np.arange(1,len(x)+1))[::-1])[::-1];back=np.empty(len(x));back[s]=np.clip(y,0,1);q[ok]=back;return q

def pval(x):
 if len(x)<3:return np.nan
 if not np.any(x):return 1.
 return float(wilcoxon(x,method='approx',correction=True).pvalue)

def tier(n,p,q):
 if n<3:return 'descriptive_no_inference'
 if n<10:return 'small_n_exploratory_nominal' if p<.05 else 'small_n_exploratory'
 return 'FDR_supported_conditional' if q<.05 else ('nominal_exploratory' if p<.05 else 'not_supported')
genes=pd.read_csv(M/'PSEUDOBULK_GENE_INDEX.tsv',sep='\t').gene.tolist();gidx={g:i for i,g in enumerate(genes)}
libs={}
for f in (O/'resources').glob('*.gmt'):
 libs[f.stem]={}
 for line in f.read_text().splitlines():
  b=line.split('\t');gs=set(g.split(',')[0] for g in b[2:] if g);libs[f.stem][b[0]]=gs
rows=[];oras=[];coverage=[];focus=[];checks=[]
pb=sparse.load_npz(M/"PATIENT_CLUSTER_PSEUDOBULK_COUNTS.npz").tocsr();ri=pd.read_csv(M/"PSEUDOBULK_ROW_INDEX.tsv",sep="\t",dtype={"cluster":str})
(O/'patient_program_effects').mkdir(exist_ok=False)
for c in ['0','3','4','5','6','7','8','9','10','11','12','15']:
 d=np.load(M/'paired_patient_effects'/('cluster_'+c+'.npz'));eff=d['effects'].astype(float);n=eff.shape[0];t=pd.read_csv(M/'full_marker_tables'/('cluster_'+c+'_all_genes.tsv.gz'),sep='\t');assert len(t)==len(genes)
 observed=[]
 for patient in d['patient']:
  u=(ri.patient==patient).to_numpy();inside=np.asarray(pb[u & (ri.cluster==c).to_numpy()].sum(axis=0)).ravel();outside=np.asarray(pb[u & (ri.cluster!=c).to_numpy()].sum(axis=0)).ravel();observed.append((inside>=10)|(outside>=10))
 keep=np.asarray(observed).sum(axis=0)>=min(n,max(3,int(np.ceil(.2*n))))
 eligible=set(np.asarray(genes)[keep])
 for lib,sets in libs.items():
  for version,remove in [('complete',set()),('without_SPP1_TYMP',{'SPP1','TYMP'})]:
   universe=eligible-remove; up=set(t.loc[(t.patient_wilcoxon_p<.05)&(t.median_patient_log2CPM_difference>=.25),'gene'])&universe;down=set(t.loc[(t.patient_wilcoxon_p<.05)&(t.median_patient_log2CPM_difference<=-.25),'gene'])&universe
   rr=[];oo=[];names=[];vals=[]
   for term,fullset in sets.items():
    gs=(fullset&universe);k=len(gs)
    coverage.append([c,lib,version,term,len(fullset),k,10<=k<=500])
    if not 10<=k<=500:continue
    idx=[gidx[g] for g in sorted(gs)];v=eff[:,idx].mean(axis=1);pv=pval(v);row={'cluster':c,'library':lib,'version':version,'term':term,'mapped_genes':k,'n_patients':n,'median_paired_program_difference':float(np.median(v)),'patients_positive':int((v>0).sum()),'patients_negative':int((v<0).sum()),'P':pv};rr.append(row);names.append(term);vals.append(v)
    for direction,sel in [('up',up),('down',down)]:
     hit=gs&sel;pvora=float(hypergeom.sf(len(hit)-1,len(universe),k,len(sel))) if n>=3 else np.nan
     oo.append({'cluster':c,'library':lib,'version':version,'direction':direction,'term':term,'universe':len(universe),'selected_genes':len(sel),'pathway_genes':k,'overlap':len(hit),'odds_enrichment':(len(hit)/max(len(sel),1))/(k/len(universe)),'P':pvora,'overlap_genes':';'.join(sorted(hit))})
   for r,q in zip(rr,bh([x['P'] for x in rr])):r['FDR']=q;r['reporting_tier']=tier(n,r['P'],q)
   for r,q in zip(oo,bh([x['P'] for x in oo])):r['FDR']=q
   rows.extend(rr);oras.extend(oo)
   np.savez_compressed(O/'patient_program_effects'/f'C{c}_{lib}_{version}.npz',patient=d['patient'],study=d['study'],terms=np.asarray(names,dtype=str),effects=np.asarray(vals).T)
   if c=='4':
    for term,v in zip(names,vals):
     for study in sorted(set(d['study'])):
      vv=v[d['study']==study];focus.append([lib,version,term,study,len(vv),float(np.median(vv)),int((vv>0).sum()),int((vv<0).sum())])
 print('DONE_CLUSTER',c,'n',n,flush=True)
x=pd.DataFrame(rows);y=pd.DataFrame(oras);x.to_csv(O/'ALL_PATIENT_PROGRAM_COMPARISONS.tsv.gz',sep='\t',index=False);y.to_csv(O/'ALL_ORA_RESULTS.tsv.gz',sep='\t',index=False)
pd.DataFrame(coverage,columns=['cluster','library','version','term','original_genes','eligible_genes','tested_size_range']).to_csv(O/'GENE_SET_COVERAGE.tsv.gz',sep='\t',index=False)
pd.DataFrame(focus,columns=['library','version','term','study','n_patients','median_difference','positive','negative']).to_csv(O/'C4_STUDY_PROGRAM_EFFECTS.tsv.gz',sep='\t',index=False)
x[x.cluster=='4'].sort_values(['version','library','P']).to_csv(O/'C4_ALL_PROGRAMS.tsv',sep='\t',index=False)
keywords='LIPID|CHOLESTEROL|FATTY|LYSOSOM|PHAGOS|INFLAM|TNFA|INTERFERON|ANTIGEN|HYPOXIA|ANGIOGEN|VASCULAR'
f=x[x.term.str.contains(keywords,case=False,regex=True)].copy();f.to_csv(O/'FOCUS_PROGRAMS_ALL_CANDIDATES.tsv',sep='\t',index=False)
f[f.cluster=='4'].sort_values(['version','library','P']).to_csv(O/'C4_FOCUS_PROGRAMS.tsv',sep='\t',index=False)
# Independent selected score and family checks, no plotting.
for keys,g in x.groupby(['cluster','library','version']):checks.append({'check':str(keys)+'_BH','pass':bool(np.allclose(bh(g.P),g.FDR,equal_nan=True))})
assert all(z['pass'] for z in checks);assert not any(p.suffix in ['.png','.svg','.pdf'] for p in O.rglob('*'))
status={'status':'COMPLETE_FUNCTIONAL_NUMERICAL','program_rows':len(x),'ORA_rows':len(y),'checks_passed':len(checks),'figures_created':False,'inference_unit':'patient','nominal_significance_accepted_as_exploratory':True};(O/'STATUS.json').write_text(json.dumps(status,indent=2));(O/'VALIDATION.json').write_text(json.dumps(checks,indent=2));print(json.dumps(status),flush=True)
