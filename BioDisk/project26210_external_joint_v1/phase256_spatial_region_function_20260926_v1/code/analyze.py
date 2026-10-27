from pathlib import Path
import sys,json,datetime,hashlib,subprocess,re
import numpy as np,pandas as pd
from scipy import sparse,stats
P=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(P/'code'))
from extract import gate,stamp,fingerprint
MODELS=['absolute','adjusted','macro_share'];EXCLUDE={'TYMP','SPP1','GPNMB','LGALS3','FABP5','CD9','LPL'}
LIBRARIES={'GO_BP':'GO_BP_local.gmt','Hallmark':'MSigDB_Hallmark_2020.gmt','KEGG':'KEGG_2021_Human.gmt'}
def family(model,version):return 'primary' if model=='absolute' and version=='exclude_definition_genes' else 'sensitivity'
def es(rank,hit):
 hit=np.sort(np.asarray(hit,int));n=len(rank);k=len(hit);w=np.abs(rank[hit]);sw=w.sum()
 if not k or k>=n or sw==0:return np.nan
 top=np.cumsum(w)/sw-(hit-np.arange(k))/(n-k)
 bottom=np.r_[0,np.cumsum(w)[:-1]]/sw-(hit-np.arange(k))/(n-k)
 mx=top.max();mn=bottom.min();return float(mx if mx>=-mn else mn)
def main():
 gate('aggregate_start')
 assert json.loads((P/'EXTRACTION_COMPLETE.json').read_text())['status']=='PASS'
 reg=pd.read_csv(P/'tables/SECTION_IDENTITY.tsv',sep='\t').fillna('')
 eligibility=pd.read_csv(P/'tables/REGION_ELIGIBILITY.tsv',sep='\t')
 libraries={}
 for lib,name in LIBRARIES.items():
  lines=(P/'resources'/name).read_text().splitlines();terms={}
  for line in lines:
   f=line.rstrip('\r\n').split('\t');genes=set(g for g in f[2:] if g)
   if f[0] in terms:raise AssertionError(('duplicate term',lib,f[0]))
   terms[f[0]]=genes
  libraries[lib]=terms
 (P/'resources/GENESET_MANIFEST.json').write_text(json.dumps({lib:dict(sets=len(s),**fingerprint(P/'resources'/LIBRARIES[lib])) for lib,s in libraries.items()},indent=2))
 sample_arrays={s:dict(np.load(P/'effects'/f'{s}.npz',allow_pickle=False)) for s in reg['sample']}
 summary=[];coverage=[];members=[];jobs=[];cohortqc=[];rankqc=[]
 for dataset,df in reg.groupby('dataset',sort=True):
  for mi,model in enumerate(MODELS):
   gate('aggregate_'+dataset+'_'+model)
   ss=eligibility[(eligibility.dataset==dataset)&(eligibility.model==model)&eligibility.eligible]['sample'].tolist()
   if not ss:continue
   common=set.intersection(*[set(sample_arrays[s]['genes']) for s in ss]);genes=np.asarray(sorted(common))
   eff=[];det=[];pb=[]
   for s in ss:
    arr=sample_arrays[s];ix=pd.Index(arr['genes']).get_indexer(genes);assert (ix>=0).all()
    eff.append(arr['effects'][mi,ix]);det.append(arr['detected_spots'][mi,ix]);pb.append(arr['pseudobulk_log2CPM_difference'][mi,ix])
   eff=np.asarray(eff);det=np.asarray(det);pb=np.asarray(pb)
   keep=np.mean(det>=10,axis=0)>=.5
   genes=genes[keep];eff=eff[:,keep];pb=pb[:,keep]
   assert np.isfinite(eff).all()
   sid=pd.DataFrame({'sample':ss}).merge(df[['sample','patient']],on='sample',validate='one_to_one')
   patients=sorted(sid.patient.unique());E=[];PB=[]
   for patient in patients:
    ix=np.flatnonzero(sid.patient.eq(patient));E.append(eff[ix].mean(0));PB.append(pb[ix].mean(0))
    members.append(dict(dataset=dataset,model=model,patient=patient,sections=';'.join(sid.iloc[ix]['sample']),n_sections=len(ix)))
   E=np.asarray(E);PB=np.asarray(PB);n=len(patients);mu=E.mean(0);sd=E.std(0,ddof=1) if n>1 else np.full(len(genes),np.nan)
   rank=mu.copy() if n<10 else np.divide(mu,sd/np.sqrt(n),out=np.zeros(len(mu)),where=sd>1e-12)
   bad=(sd<=1e-12)&(np.abs(mu)>1e-12) if n>=10 else np.zeros(len(mu),bool)
   rank[bad]=np.nan
   table=pd.DataFrame({'gene':genes,'patient_n':n,'mean_Z_effect':mu,'median_Z_effect':np.median(E,0),'positive_patients':(E>1e-12).sum(0),'negative_patients':(E< -1e-12).sum(0),'rank':rank,'excluded_definition_gene':[g in EXCLUDE for g in genes],'median_pseudobulk_log2CPM_difference':np.median(PB,0)})
   table.to_csv(P/'tables'/f'GENE_RANK_{dataset}_{model}.tsv.gz',sep='\t',index=False)
   np.savez_compressed(P/'effects'/f'PATIENT_{dataset}_{model}.npz',genes=genes,patients=np.asarray(patients),effects=E,pseudobulk=PB)
   cohortqc.append(dict(dataset=dataset,model=model,patients=n,sections=len(ss),common_measured_genes=len(common),expression_eligible_genes=len(genes),non_estimable_constant_rank=int(bad.sum())))
   for version in ['exclude_definition_genes','complete']:
    gate('pathways_'+dataset+'_'+model+'_'+version)
    good=np.isfinite(rank)&np.array([g not in EXCLUDE or version=='complete' for g in genes])
    gg=genes[good];EE=E[:,good];rr=rank[good];lookup={g:i for i,g in enumerate(gg)}
    order=np.lexsort((gg,-rr));rg=gg[order];rv=rr[order];ranklookup={g:i for i,g in enumerate(rg)}
    tag=f'{dataset}__{model}__{version}'
    pd.DataFrame({'gene':rg,'rank':rv}).to_csv(P/'gsea'/f'{tag}.rnk',sep='\t',index=False)
    rankqc.append(dict(tag=tag,genes=len(rg),zero_ranks=int((rv==0).sum()),nonzero_ties=int(len(rv[rv!=0])-len(np.unique(rv[rv!=0]))),ranking='paired_patient_t_metric' if n>=10 else 'mean_patient_Z_effect_descriptive'))
    for lib,sets in libraries.items():
     terms=[];gene_lists=[];rows=[];cols=[];vals=[];gmt=[]
     for term,original in sorted(sets.items()):
      allowed=original-EXCLUDE if version=='exclude_definition_genes' else original
      used=sorted(allowed&set(lookup));num=len(used);frac=num/len(allowed) if allowed else 0
      ok=10<=num<=500 and frac>=.2
      coverage.append(dict(dataset=dataset,model=model,version=version,library=lib,term=term,original_size=len(original),after_exclusion_size=len(allowed),measured_size=num,coverage=frac,eligible=ok))
      if not ok:continue
      j=len(terms);terms.append(term);gene_lists.append(used)
      ix=[lookup[g] for g in used];rows.extend(ix);cols.extend([j]*len(ix));vals.extend([1/len(ix)]*len(ix))
      gmt.append('\t'.join([term,'frozen_measured_members']+used))
     if not terms:continue
     W=sparse.csc_matrix((vals,(rows,cols)),shape=(len(gg),len(terms)));D=np.asarray(EE@W)
     assert D.shape==(n,len(terms)) and np.isfinite(D).all()
     pos=(D>1e-12).sum(0);neg=(D< -1e-12).sum(0);nz=pos+neg
     pv=np.minimum(1,2*stats.binom.cdf(np.minimum(pos,neg),nz,.5)) if n>=10 else np.full(len(terms),np.nan)
     if n>=10:pv[nz==0]=1
     med=np.median(D,0);mean=D.mean(0)
     for j,term in enumerate(terms):
      summary.append(dict(dataset=dataset,model=model,version=version,family=family(model,version),library=lib,term=term,patient_n=n,genes=len(gene_lists[j]),mean_Z_difference=mean[j],median_Z_difference=med[j],positive=int(pos[j]),negative=int(neg[j]),zero=int(n-nz[j]),P_patient_sign=pv[j],GSEA_ES_descriptive=es(rv,[ranklookup[g] for g in gene_lists[j]]) if n<10 else np.nan))
     np.savez_compressed(P/'effects'/f'PROGRAM_{tag}__{lib}.npz',terms=np.asarray(terms),patients=np.asarray(patients),effects=D)
     (P/'gsea'/f'{tag}__{lib}.gmt').write_text('\n'.join(gmt)+'\n')
     jobs.append(dict(dataset=dataset,model=model,version=version,library=lib,family=family(model,version),patient_n=n,rank_path=str(P/'gsea'/f'{tag}.rnk'),gmt_path=str(P/'gsea'/f'{tag}__{lib}.gmt'),output=str(P/'gsea'/f'RESULT_{tag}__{lib}.tsv')))
 su=pd.DataFrame(summary);su['FDR_patient_sign']=np.nan
 for key,ids in su.groupby(['dataset','library','family']).groups.items():
  ids=np.asarray(list(ids));p=su.loc[ids,'P_patient_sign'];ok=p.notna()
  su.loc[ids,'registered_family_size']=len(ids)
  if ok.any():su.loc[ids[ok],'FDR_patient_sign']=stats.false_discovery_control(p[ok].to_numpy())
 su['patient_evidence']=np.where(su.P_patient_sign.isna(),'DESCRIPTIVE_N_LT_10',np.where(su.FDR_patient_sign<.05,'FDR_SUPPORTED',np.where(su.P_patient_sign<.05,'NOMINAL_ONLY','NOT_SUPPORTED')))
 su.to_csv(P/'tables/PATIENT_FUNCTION_SUMMARY.tsv.gz',sep='\t',index=False)
 pd.DataFrame(coverage).to_csv(P/'tables/GENESET_COVERAGE.tsv.gz',sep='\t',index=False)
 pd.DataFrame(members).to_csv(P/'tables/PATIENT_MEMBERSHIP.tsv',sep='\t',index=False)
 pd.DataFrame(cohortqc).to_csv(P/'qc/COHORT_GENE_UNIVERSE.tsv',sep='\t',index=False)
 pd.DataFrame(rankqc).to_csv(P/'qc/RANK_AUDIT.tsv',sep='\t',index=False)
 pd.DataFrame(jobs).to_csv(P/'gsea/JOBS.tsv',sep='\t',index=False)
 (P/'AGGREGATION_COMPLETE.json').write_text(json.dumps(dict(status='PASS',time=stamp(),program_rows=len(su),gsea_jobs=sum(j['patient_n']>=10 for j in jobs),descriptive_jobs=sum(j['patient_n']<10 for j in jobs)),indent=2))
 print('AGGREGATION_COMPLETE',len(su),flush=True)
if __name__=='__main__':main()
