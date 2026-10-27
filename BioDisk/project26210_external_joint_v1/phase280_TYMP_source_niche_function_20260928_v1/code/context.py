from source import *
def main():
 src=R/'phase256_spatial_region_function_20260926_v1/resources/MSigDB_Hallmark_2020.gmt';m=read(M/'inputs/MACRO_METADATA.tsv.gz',index_col=0,dtype={'candidate':str});X=sparse.load_npz(M/'inputs/MACRO_COUNTS.npz').tocsr();genes=read(M/'inputs/GENES.tsv').gene.tolist();gi={g:i for i,g in enumerate(genes)}
 targets=json.loads((P/'resources/PROGRAMS.json').read_text())
 excluded=set(['TYMP','SPP1','LGALS3','CTSB','CTSD']+sum(list(targets.values()),[]))
 names=['Interferon Alpha Response','Interferon Gamma Response','Hypoxia'];sets={l.split('\t')[0]:[g for g in l.split('\t')[2:] if g not in excluded and g in gi] for l in src.read_text().splitlines() if l.split('\t')[0] in names}
 dump(P/'resources/CONTEXT_PROGRAMS.json',dict(genes=sets,excluded=sorted(excluded),source=str(src),sha256=sha(src)))
 union=sorted(set(sum(list(sets.values()),[])));raw=X[:,[gi[g] for g in union]].toarray();total=np.asarray(X.sum(1)).ravel();log=np.log1p(raw/total[:,None]*1e4);index={g:i for i,g in enumerate(union)}
 cache=read(P/'inputs/CELL_CACHE.tsv.gz',dtype={'candidate':str});ef=read(P/'tables/WITHIN_MACRO_PATIENT_EFFECTS.tsv');rows=[];ca=[]
 for pat,z in cache.groupby('patient'):
  ids=m.index.get_indexer(z.cell_id);assert (ids>=0).all();scores=[];eligible=np.ones(len(z),bool)
  rec=z[['cell_id','patient']].copy()
  for name in names:
   jj=[index[g] for g in sets[name]];v=log[ids][:,jj];sd=v.std(0);sc=((v-v.mean(0))/np.where(sd>0,sd,1)).mean(1);scores.append(sc);rec[name]=sc
   # Each context must be measurable across at least 20 cells with >=2 genes.
   eligible &= np.repeat(int(((raw[ids][:,jj]>0).sum(1)>=2).sum())>=20,len(z))
  ca.append(rec)
  for r in ef[(ef.patient==pat)&(ef.model=='adjusted')].itertuples():
   mask=z.candidate.eq('4').to_numpy() if r.group=='C4' else np.ones(len(z),bool)
   cov=np.column_stack([z.log_depth,z.n_genes]+scores)
   if r.group=='WholeMacro':cov=np.column_stack([cov,pd.get_dummies(z.candidate,dtype=float,drop_first=True).to_numpy()])
   if np.isfinite(r.effect) and eligible.all():eff,n,reason=fit(z.TYMP,z[r.endpoint],cov,mask,50 if r.group=='C4' else 100)
   else:eff,n,reason=np.nan,int(mask.sum()),'ORIGINAL_OR_CONTEXT_INELIGIBLE'
   rows.append(dict(patient=pat,study=z.study.iloc[0],group=r.group,endpoint=r.endpoint,model='depth_complexity_cluster_IFNa_IFNg_hypoxia',effect=eff,cells=int(mask.sum()),reason=reason))
 e=pd.DataFrame(rows);tab(e,'CONTEXT_SENSITIVITY_PATIENT_EFFECTS.tsv');pd.concat(ca).to_csv(P/'inputs/CONTEXT_CELL_CACHE.tsv.gz',sep='\t',index=False);su=[]
 for key,z in e.groupby(['group','endpoint']):su.append(dict(zip(['group','endpoint'],key),**infer(z.effect,280600+len(su))))
 # Single frozen sensitivity family: 2 groups x all 6 original outcomes.
 su=pd.DataFrame(su);su['family']='context_sensitivity_12';tab(bh(su,['family']),'CONTEXT_SENSITIVITY_SUMMARY.tsv')
 from validate import independent_rho
 maxerr=0.;ncheck=0
 merged=cache.merge(pd.concat(ca),on=['cell_id','patient'],validate='one_to_one')
 for r in e[e.effect.notna()].itertuples():
  z=merged[merged.patient==r.patient];mask=z.candidate.eq('4').to_numpy() if r.group=='C4' else np.ones(len(z),bool)
  cov=np.column_stack([z.log_depth,z.n_genes]+[z[name] for name in names])
  if r.group=='WholeMacro':cov=np.column_stack([cov,pd.get_dummies(z.candidate,dtype=float,drop_first=True).to_numpy()])
  rho=independent_rho(z.TYMP,z[r.endpoint],cov,mask,50 if r.group=='C4' else 100);assert abs(rho-r.effect)<1e-10
  maxerr=max(maxerr,abs(rho-r.effect));ncheck+=1
 dump(P/'qc/CONTEXT_VALIDATION.json',dict(status='PASS',independent_effects=ncheck,max_error=maxerr,source_sha256=sha(src)))
 print('CONTEXT_COMPLETE',flush=True)
if __name__=='__main__':main()
