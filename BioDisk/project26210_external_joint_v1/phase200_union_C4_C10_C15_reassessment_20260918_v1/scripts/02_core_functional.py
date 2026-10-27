from common import *
from scipy import sparse
from scipy.stats import hypergeom
import json,math
m=pd.read_csv(D/'inputs/MACRO_METADATA.tsv.gz',sep='\t',index_col=0,dtype={'candidate':str})
groups=pd.read_csv(D/'inputs/GROUP_INDEX.tsv',sep='\t',dtype={'cluster':str})
genes=pd.read_csv(D/'inputs/GENES.tsv',sep='\t').gene.to_numpy();gidx={g:i for i,g in enumerate(genes)}
counts=sparse.load_npz(D/'inputs/GROUP_COUNTS.npz').toarray()
det=sparse.load_npz(D/'inputs/GROUP_DETECTED.npz').toarray()
logsum=sparse.load_npz(D/'inputs/GROUP_LOGSUM.npz').toarray()
union=['4','10','15'];dom=json.loads((D/'inputs/FROZEN_UNION.json').read_text())['dominant_C10_patient']
allcs=set(groups.cluster);ptstudy=m.groupby('curated_patient_id').study.first().to_dict()
configs=[];store={}
for scope in ['all','tumor']:
 g=groups if scope=='all' else groups[groups.tissue.str.lower().eq('tumor')]
 pats=sorted(g.patient.unique());tot={p:counts[g.index[g.patient==p]].sum(0) for p in pats}
 nc={p:int(g.loc[g.patient==p,'cells'].sum()) for p in pats}
 specs={'U':union,'C4':['4'],'C10':['10'],'C15':['15'],'U_no10':['4','15'],'U_no15':['4','10']}
 for label,cs in specs.items():
  ep=[]
  for p in pats:
   ix=g.index[(g.patient==p)&g.cluster.isin(cs)];ni=int(g.loc[ix,'cells'].sum());nr=nc[p]-ni
   if min(ni,nr)>=20:ep.append(p)
  store[(scope,label)]={'patients':ep,'g':g,'target':cs,'background':list(allcs-set(cs))}
 common=sorted(set(store[(scope,'U')]['patients'])&set(store[(scope,'C4')]['patients']))
 for label in ['U','C4','C10','C15','U_no10','U_no15']:
  z=store[(scope,label)];configs.append(dict(name=scope+'_'+label,scope=scope,label=label,**z))
 for label in ['U','C4']:
  z=store[(scope,label)].copy();z['patients']=common
  configs.append(dict(name=scope+'_'+label+'_common',scope=scope,label=label,**z))
  z=z.copy();z['background']=list(allcs-set(union))
  z['patients']=[p for p in common if int(g.loc[(g.patient==p)&g.cluster.isin(z['background']),'cells'].sum())>=20]
  configs.append(dict(name=scope+'_'+label+'_shared_background',scope=scope,label=label,**z))
 for label in ['U','C4']:
  z=store[(scope,label)].copy();z['patients']=[p for p in z['patients'] if p!=dom]
  configs.append(dict(name=scope+'_'+label+'_without_C10_dominant_patient',scope=scope,label=label,**z))
# Patient composition and SPP1 RNA coverage.
comp=[]
for scope in ['all','tumor']:
 g=groups if scope=='all' else groups[groups.tissue.str.lower().eq('tumor')]
 for p,z in g.groupby('patient'):
  indices=z.index;total=counts[indices].sum(0);n=int(z.cells.sum());spp=total[gidx['SPP1']]
  for label,cs in {'U':union,'C4':['4'],'C10':['10'],'C15':['15'],'nonU':list(allcs-set(union))}.items():
   ix=z.index[z.cluster.isin(cs)];nn=int(z.loc[ix,'cells'].sum());vv=counts[ix].sum(0);dd=det[ix].sum(0)
   comp.append(dict(scope=scope,patient=p,study=z.study.iloc[0],group=label,cells=nn,total_macrophages=n,fraction=nn/n,SPP1_raw=vv[gidx['SPP1']],SPP1_RNA_share=vv[gidx['SPP1']]/spp if spp>0 else np.nan,SPP1_log2CPM=np.log2(vv[gidx['SPP1']]/vv.sum()*1e6+1) if vv.sum()>0 else np.nan,SPP1_detected_fraction=dd[gidx['SPP1']]/nn if nn else np.nan))
pd.DataFrame(comp).to_csv(D/'core/PATIENT_COMPOSITION_SPP1.tsv',sep='\t',index=False)
libs={}
for f in (R/'part2_rebuild_v1/06_functional_programs/run_v1/resources').glob('*.gmt'):
 sets={}
 for ln in f.read_text().splitlines():
  a=ln.split('\t');sets[a[0]]=set(v.split(',')[0] for v in a[2:] if v)
 libs[f.stem]=sets
prows=[];orows=[];reg=[];markerfocus=[];studyrows=[]
focusgenes=['SPP1','TYMP','GPNMB','CCL18','CYP27A1','LPL','APOE','MKI67','TOP2A','USP53','ALOX5AP','JUN','JUNB','CEBPD','STAT3','SOCS3','ADORA2A','ADORA2B','IL6','IL6R','IL6ST','IL10','IL10RA','IL10RB','LGALS3','SLC16A3','LDHA','SRM','SMS','MTAP']
(D/'core/contrasts').mkdir(exist_ok=True);(D/'functional/patient_effects').mkdir(exist_ok=True)
for cfg in configs:
 name=cfg['name'];g=cfg['g'];pats=cfg['patients'];n=len(pats)
 if n==0:continue
 inside=[];outside=[];detection=[];cellns=[]
 for p in pats:
  iz=g.index[(g.patient==p)&g.cluster.isin(cfg['target'])];oz=g.index[(g.patient==p)&g.cluster.isin(cfg['background'])]
  inside.append(counts[iz].sum(0));outside.append(counts[oz].sum(0))
  ni=int(g.loc[iz,'cells'].sum());no=int(g.loc[oz,'cells'].sum())
  detection.append([det[iz].sum(0)/ni,det[oz].sum(0)/no]);cellns.append([ni,no])
 a=np.asarray(inside);b=np.asarray(outside)
 ai=np.log2(a/a.sum(1)[:,None]*1e6+1);bo=np.log2(b/b.sum(1)[:,None]*1e6+1);eff=ai-bo
 observed=((a>=10)|(b>=10)).sum(0)>=min(n,max(3,math.ceil(.2*n)))
 test=observed&np.any(eff!=0,axis=0)&(n>=3)
 pv=paired_pvals(eff,test);med=np.median(eff,axis=0)
 tab=pd.DataFrame({'gene':genes,'n':n,'median_delta':med,'positive':(eff>0).sum(0),'negative':(eff<0).sum(0),'P':pv,'FDR':bh(pv),'test_eligible':test,'observed':observed})
 tab.to_csv(D/'core/contrasts'/f'{name}_genes.tsv.gz',sep='\t',index=False)
 np.savez_compressed(D/'core/contrasts'/f'{name}.npz',patients=np.asarray(pats),study=np.asarray([ptstudy[p] for p in pats]),genes=genes,effects=eff,target_log2CPM=ai,reference_log2CPM=bo,target_detection=np.asarray(detection)[:,0,:],reference_detection=np.asarray(detection)[:,1,:],cell_counts=np.asarray(cellns),tested=test,observed=observed)
 for _,row in tab[tab.gene.isin(focusgenes)].iterrows():markerfocus.append(dict(contrast=name,scope=cfg['scope'],group=cfg['label'],**row.to_dict()))
 # Study effects retained without treating study partitions as independent patients.
 for st in sorted(set(ptstudy[p] for p in pats)):
  ii=np.array([ptstudy[p]==st for p in pats])
  for gene in focusgenes:
   if gene in gidx:studyrows.append(dict(contrast=name,study=st,gene=gene,**summary(eff[ii,gidx[gene]])))
 reg.append(dict(contrast=name,scope=cfg['scope'],group=cfg['label'],target=';'.join(cfg['target']),background=';'.join(sorted(cfg['background'])),n=n,patient_list=';'.join(pats),genes_tested=int(test.sum())))
 # Full functional libraries, preserving original per-library/contrast/version family.
 for lib,sets in libs.items():
  for version,remove in [('complete',set()),('without_SPP1_TYMP',{'TYMP','SPP1'})]:
   universe=set(genes[observed])-remove;names=[];vals=[];rr=[];oo=[]
   up=set(tab.loc[(tab.P<.05)&(tab.median_delta>=.25),'gene'])&universe
   down=set(tab.loc[(tab.P<.05)&(tab.median_delta<=-.25),'gene'])&universe
   for term,full in sets.items():
    gs=full&universe
    if not 10<=len(gs)<=500:continue
    v=eff[:,[gidx[g] for g in sorted(gs)]].mean(1);names.append(term);vals.append(v)
    rr.append(dict(contrast=name,library=lib,version=version,term=term,mapped_genes=len(gs),**summary(v)))
    if name in ['all_U','all_C4','tumor_U','tumor_C4']:
     for direction,selected in [('up',up),('down',down)]:
      hits=selected&gs;p=float(hypergeom.sf(len(hits)-1,len(universe),len(gs),len(selected))) if n>=3 else np.nan
      oo.append(dict(contrast=name,library=lib,version=version,term=term,direction=direction,universe=len(universe),selected=len(selected),overlap=len(hits),P=p,genes=';'.join(sorted(hits))))
   for row,q in zip(rr,bh([v['P'] for v in rr])):row['FDR']=q
   for row,q in zip(oo,bh([v['P'] for v in oo])):row['FDR']=q
   prows.extend(rr);orows.extend(oo)
   np.savez_compressed(D/'functional/patient_effects'/f'{name}__{lib}__{version}.npz',patients=np.asarray(pats),terms=np.asarray(names),effects=np.asarray(vals).T)
 print('CONTRAST_DONE',name,'n',n,flush=True)
pd.DataFrame(reg).to_csv(D/'core/CONTRAST_REGISTRY.tsv',sep='\t',index=False)
pd.DataFrame(markerfocus).to_csv(D/'core/KEY_MARKER_COMPARISONS.tsv',sep='\t',index=False)
pd.DataFrame(studyrows).to_csv(D/'core/STUDY_KEY_MARKERS.tsv',sep='\t',index=False)
fun=pd.DataFrame(prows);fun.to_csv(D/'functional/ALL_PROGRAM_RESULTS.tsv.gz',sep='\t',index=False)
fun[fun.term.str.contains('CHOLESTEROL|LYSOSOM|PHAGOS|ANTIGEN|GLYCOL|HYPOXIA|E2F|G2|FATTY|ANGIOGEN|INFLAM|INTERFERON',case=False)].to_csv(D/'functional/FOCUS_PROGRAM_RESULTS.tsv',sep='\t',index=False)
pd.DataFrame(orows).to_csv(D/'functional/ORA_ALL_RESULTS.tsv.gz',sep='\t',index=False)
# Check original C4 marker estimates before proceeding.
orig=pd.read_csv(R/'part2_rebuild_v1/04_markers_annotation/run_v4_macrophage_only/full_marker_tables/cluster_4_all_genes.tsv.gz',sep='\t')
new=pd.read_csv(D/'core/contrasts/all_C4_genes.tsv.gz',sep='\t')
assert np.array_equal(orig.gene,new.gene) and np.allclose(orig.median_patient_log2CPM_difference,new.median_delta,atol=1e-6)
assert np.allclose(orig.patient_wilcoxon_p,new.P,equal_nan=True,atol=1e-10)
(D/'qc/CORE_FUNCTIONAL_STATUS.json').write_text(json.dumps({'status':'PASS','contrasts':len(reg),'functional_rows':len(fun),'original_C4_gene_effects_P_reproduced':True,'no_plots':True},indent=2))
print('CORE_FUNCTIONAL_COMPLETE',flush=True)
