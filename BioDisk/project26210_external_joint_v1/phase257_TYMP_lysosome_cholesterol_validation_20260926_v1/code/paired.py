from common import *
from scipy import stats
from itertools import permutations,product
from math import factorial
EXCLUDE={'TYMP','SPP1','GPNMB','LGALS3','FABP5','CD9','LPL'}
FUNCTIONS={'Lysosome':'lys','Cholesterol metabolism':'chol'}
def permutation_test(detail,studies,seed):
 a=detail['a'];b=detail['b'];Q=detail['Q'];studies=np.asarray(studies)[detail['keep']]
 groups=[np.flatnonzero(studies==s) for s in sorted(set(studies))]
 total=1
 for ix in groups:total*=factorial(len(ix))
 n=len(a);exact=total<=20000
 if exact:
  choices=[list(permutations(ix.tolist())) for ix in groups]
  perm=np.tile(np.arange(n),(total,1))
  for i,picks in enumerate(product(*choices)):
   for ix,pick in zip(groups,picks):perm[i,ix]=pick
 else:
  rng=np.random.default_rng(seed);perm=np.tile(np.arange(n),(19999,1))
  for ix in groups:perm[:,ix]=ix[np.argsort(rng.random((19999,len(ix))),axis=1)]
 A=a[perm];A=A-(A@Q)@Q.T
 denom=np.sqrt(np.sum(A*A,axis=1)*np.dot(b,b));null=(A@b)/denom
 valid=np.isfinite(null);assert valid.all()
 observed=float(np.dot(a,b)/np.sqrt(np.dot(a,a)*np.dot(b,b)))
 hits=int((np.abs(null)>=abs(observed)-1e-12).sum())
 pv=hits/len(null) if exact else (hits+1)/(len(null)+1)
 # Independently recompute first fixed 50 permutations, row by row.
 err=0.
 for i in range(min(50,len(perm))):
  v=a[perm[i]];v=v-Q@(Q.T@v);ref=float(stats.pearsonr(v,b).statistic);err=max(err,abs(ref-null[i]))
 assert err<1e-10
 return dict(P=pv,permutations=len(null),exact_permutation=exact,possible_permutations=str(total),movable_patients=sum(len(i) for i in groups if len(i)>1),permutation_max_check_error=err,seed=seed)
def main():
 gate('paired_start')
 source=ROOT/'phase154_same_library_TYMP_C4_test_20260915_v1/SAME_LIBRARY_PATIENT_METRICS.tsv'
 groupfile=ROOT/'phase157_ordered_bridge_reassessment_20260915_v1/02_inputs/GROUP_INDEX.tsv'
 arrayfile=ROOT/'phase157_ordered_bridge_reassessment_20260915_v1/02_inputs/GROUP_PSEUDOBULK.npz'
 gmt=ROOT/'phase256_spatial_region_function_20260926_v1/resources/KEGG_2021_Human.gmt'
 paths=[source,groupfile,arrayfile,gmt];before=[fingerprint(f) for f in paths]
 b=pd.read_csv(source,sep='\t').reset_index(drop=True);G=pd.read_csv(groupfile,sep='\t')
 assert len(b)==16 and b.patient.is_unique and b[['patient','sample_id']].duplicated().sum()==0
 Z=np.load(arrayfile,allow_pickle=True);genes=np.asarray(Z['genes'],dtype=str);counts=Z['counts']
 assert len(set(genes))==len(genes) and np.isfinite(counts).all() and (counts>=0).all() and np.allclose(counts,np.round(counts))
 gi={g:i for i,g in enumerate(genes)};arrays={};audit=[]
 for compartment,roles in [('tumor',['Malignant']),('macro',['C4','OtherMacrophage'])]:
  values=[]
  for r in b.itertuples():
   rows=G[(G.patient==r.patient)&(G.sample_id==r.sample_id)&G.role.isin(roles)]
   assert len(rows)>0 and set(rows.study)=={r.study} and set(rows.source_tissue)=={'Tumor'}
   c=counts[rows.group_index.to_numpy(int)].sum(0);total=float(c.sum());cells=int(rows.cells.sum())
   original_total=r.malignant_total_counts if compartment=='tumor' else r.macrophage_total_counts
   original_cells=r.malignant_cells if compartment=='tumor' else r.macrophage_cells
   assert cells==original_cells and total==original_total,(r.patient,compartment,total,original_total)
   if compartment=='tumor':assert c[gi['TYMP']]==r.malignant_TYMP_counts
   else:assert c[gi['SPP1']]==r.macrophage_SPP1_counts
   values.append(c);audit.append(dict(compartment=compartment,patient=r.patient,sample_id=r.sample_id,cells=cells,total_counts=total,TYMP_counts=float(c[gi['TYMP']]),SPP1_counts=float(c[gi['SPP1']])))
  C=np.asarray(values);CPM=C/C.sum(1)[:,None]*1e6
  arrays[compartment]=dict(counts=C,CPM=CPM,log=np.log2(CPM+1))
  for g in ['TYMP','SPP1']:b[compartment+'_'+g]=arrays[compartment]['log'][:,gi[g]]
 np.testing.assert_allclose(b.tumor_TYMP,b.TYMP_pseudobulk,atol=1e-10)
 np.testing.assert_allclose(b.macro_SPP1,b.SPP1_pseudobulk,atol=1e-10)
 sets={f[0]:set(f[2:])-EXCLUDE for f in [l.split('\t') for l in gmt.read_text().splitlines()] if f[0] in FUNCTIONS}
 shared=sets['Lysosome']&sets['Cholesterol metabolism'];coverage=[]
 for comp,arr in arrays.items():
  for term,short in FUNCTIONS.items():
   original=sets[term]
   for version in ['primary','no_shared']:
    allowed=original-shared if version=='no_shared' else original
    used=sorted(g for g in allowed if g in gi and np.mean(arr['CPM'][:,gi[g]]>=1)>=.5 and arr['log'][:,gi[g]].std()>1e-10)
    ok=len(used)>=10 and len(used)>=len(allowed)*.2
    coverage.append(dict(compartment=comp,function=term,version=version,listed=len(allowed),used=len(used),coverage=len(used)/len(allowed),eligible=ok,members=';'.join(used),removed_shared=';'.join(sorted(shared))))
    key=comp+'_'+short+('' if version=='primary' else '__no_shared')
    if ok:
     E=arr['log'][:,[gi[g] for g in used]];s=(E-E.mean(0))/E.std(0);b[key]=s.mean(1)
     # Independent reconstruction from per-gene series.
     check=np.column_stack([(arr['log'][:,gi[g]]-arr['log'][:,gi[g]].mean())/arr['log'][:,gi[g]].std() for g in used]).mean(1)
     np.testing.assert_allclose(b[key],check,atol=1e-12)
    else:b[key]=np.nan
 out('tables/PAIRED_COUNTS_AUDIT.tsv',pd.DataFrame(audit));out('tables/PAIRED_PROGRAM_COVERAGE.tsv',pd.DataFrame(coverage));out('scores/PAIRED_PATIENT_VALUES.tsv',b)
 endpoints=[]
 for f in ['lys','chol']:
  endpoints.extend([('tumor_intracellular','tumor_TYMP','tumor_'+f),('tumor_to_macro','tumor_TYMP','macro_'+f),('macro_SPP1_function','macro_SPP1','macro_'+f),('macro_TYMP_function','macro_TYMP','macro_'+f)])
 endpoints.append(('tumor_TYMP_macro_SPP1','tumor_TYMP','macro_SPP1'))
 out('inputs/PAIRED_FROZEN_ENDPOINTS.tsv',pd.DataFrame(endpoints,columns=['family','x','y']))
 study=pd.get_dummies(b.study,drop_first=True,dtype=float);rows=[];loo=[]
 for fam,x,y in endpoints:
  gate('paired_'+fam+'_'+y)
  for adj,version,role in [('study','primary','primary'),('raw','primary','sensitivity'),('study_covariate','primary','sensitivity'),('study','no_shared','sensitivity')]:
   if version=='no_shared' and fam=='tumor_TYMP_macro_SPP1':continue
   yy=y+'__no_shared' if version=='no_shared' else y
   cv={} if adj=='raw' else {str(c):study[c].to_numpy() for c in study}
   if adj=='study_covariate':
    cv['log1p_tumor_counts' if fam=='tumor_intracellular' else 'C4_fraction']=np.log1p(b.malignant_total_counts).to_numpy() if fam=='tumor_intracellular' else b.C4_fraction.to_numpy()
   info,detail=correlation(b[x],b[yy],cv,min_n=10)
   r=dict(family=fam,x=x,y=y,adjustment=adj,version=version,role=role,patient_n=len(b),**info)
   if detail is None:r['P']=np.nan
   elif adj=='raw':
    r['P']=float(stats.spearmanr(b[x],b[yy]).pvalue);r['permutations']=0;r['P_method']='patient_Spearman'
   else:
    seed=int(hashlib.sha256((fam+x+y+adj+version).encode()).hexdigest()[:8],16)^262257
    r.update(permutation_test(detail,b.study.to_numpy(),seed));r['P_method']='within_study_residual_permutation'
   rows.append(r)
   if adj=='study' and version=='primary':
    for omitted in sorted(b.study.unique()):
     sub=b[b.study!=omitted];D=pd.get_dummies(sub.study,drop_first=True,dtype=float)
     rv,_=correlation(sub[x],sub[y],{str(c):D[c].to_numpy() for c in D},min_n=5)
     loo.append(dict(family=fam,x=x,y=y,omitted_study=omitted,patient_n=len(sub),**rv))
 result=correction(pd.DataFrame(rows),['family','role']);out('tables/PAIRED_ASSOCIATIONS.tsv',result);out('tables/PAIRED_LEAVE_STUDY.tsv',pd.DataFrame(loo))
 assert before==[fingerprint(f) for f in paths]
 jwrite('qc/PAIRED_SOURCE_MANIFEST.json',before)
 jwrite('PAIRED_COMPLETE.json',dict(status='PASS',time=stamp(),patients=len(b),studies=b.study.value_counts().to_dict(),endpoints=len(endpoints),associations=len(result),QR_checks=len(CHECKS),max_QR_error=max(CHECKS),programs_eligible=int(pd.DataFrame(coverage).eligible.sum())))
 print(result[result.role=='primary'][['family','x','y','patient_n','rho','P','FDR']].to_string(index=False),flush=True)
if __name__=='__main__':main()
