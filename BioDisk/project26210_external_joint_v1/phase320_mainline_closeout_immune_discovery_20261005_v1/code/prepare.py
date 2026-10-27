from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='2'
import json,subprocess,hashlib,gc
import numpy as np,pandas as pd,anndata as ad
P=Path(__file__).resolve().parents[1];R=P.parent;S=R/'phase305_larger_immune_sequential_20261002_v1';pr=json.loads((P/'resources/ANALYSIS_PROTOCOL.json').read_text());panels=pr['programs'];assert len(panels)==50 and all(v['collection']=='Hallmark2020' for v in panels.values())
old=pd.read_csv(S/'tables/PATIENT_MEASURES.tsv',sep='\t');sources=old[old.measure.isin(['TYMP','SPP1'])].copy();sources=sources[(sources.lineage=='Macrophage')|((sources.lineage=='Tumour_like')&(sources.measure=='TYMP'))].copy();sources['eligible_strict']=sources.eligible;sources['eligible_primary']=(sources.cells>=30)&(sources.genes_measured==1)&(sources.detected_cells_or_genes>=3);sources['source']=sources.lineage+':'+sources.measure;sources.to_csv(P/'tables/SOURCE_PATIENT_MEASURES.tsv',sep='\t',index=False)
out=[];cover=[];members=[];cachemanifest=[];identities=[]
def sha(f):
 h=hashlib.sha256()
 with open(f,'rb') as z:
  for b in iter(lambda:z.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
for ds in ['OBRADOVIC','HCRN']:
 q=subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True);(P/'logs'/f'memory_gate_{ds}.log').write_text(q.stdout+q.stderr);assert q.returncode==0
 hpath=S/f'inputs/{ds}_counts_qc.h5ad';mpath=S/f'tables/{ds}_QUALIFIED_CELL_METADATA.tsv.gz';identities.extend([dict(path=str(f),bytes=f.stat().st_size,sha256=sha(f)) for f in [hpath,mpath]])
 a=ad.read_h5ad(hpath);m=pd.read_csv(mpath,sep='\t',index_col=0);assert a.obs_names.equals(m.index);x=a.X.tocsr();present={}
 for lib in m.library.unique():
  b=ad.read_h5ad(S/f'inputs/{ds}_{lib}_qc.h5ad',backed='r');present[lib]=set(b.var_names);b.file.close()
 for unit,g in m.groupby('analysis_unit',sort=True):
  ui=m.analysis_unit.eq(unit).to_numpy();sel=ui&m.lineage_CD8.to_numpy(bool);msel=ui&m.lineage_Macrophage.to_numpy(bool);assert not (sel&msel).any();n=int(sel.sum());nm=int(msel.sum());idx=np.flatnonzero(sel);libs=m.loc[sel,'library'].unique();available=set.intersection(*[present[l] for l in libs]) if len(libs) else set();summ=np.asarray(x[sel].sum(axis=0)).ravel().astype(float);det=np.asarray((x[sel]>0).sum(axis=0)).ravel();pb=np.log1p(1e4*summ/max(1,summ.sum()));oldn=old[(old.dataset==ds)&(old.analysis_unit==unit)&(old.lineage=='CD8')].cells.unique();assert len(oldn)==1 and oldn[0]==n
  b=sources[(sources.dataset==ds)&(sources.analysis_unit==unit)].iloc[0];base={k:b[k] for k in ['dataset','analysis_unit','patient','stratum','technical_batch','pT_category']};base['CD8_cells']=n;base['macrophage_cells']=nm;base['receiver_median_UMI']=float(m.loc[sel,'total_counts'].median()) if n else np.nan
  for name,pp in panels.items():
   gs=pp['genes'];keep=[z for z in gs if z in available];jj=a.var_names.get_indexer(keep);assert (jj>=0).all();measured=len(jj);seen=int((det[jj]>=3).sum());strictseen=int((det[jj]>=max(5,int(np.ceil(.05*n)))).sum());ok=n>=20 and measured>=.8*len(gs) and seen>=.2*len(gs);oks=ok and n>=100 and strictseen>=.5*len(gs)
   reason=[]
   if n<20:reason.append('CD8_LT20')
   if measured<.8*len(gs):reason.append('MEASURED_LT80PCT')
   if seen<.2*len(gs):reason.append('DETECTED_LT20PCT')
   score=float(pb[jj].mean()) if measured else np.nan
   out.append(base|dict(endpoint=name,family='CD8_Hallmark',value=score,eligible_primary=bool(ok),eligible_strict=bool(oks),genes_expected=len(gs),genes_measured=measured,genes_detected=seen,genes_detected_strict=strictseen,reason=';'.join(reason)))
   for gene,j in zip(keep,jj):cover.append(dict(dataset=ds,analysis_unit=unit,endpoint=name,gene=gene,sum_count=float(summ[j]),detected_cells=int(det[j]),logCP10K=float(pb[j])))
  out.append(base|dict(endpoint='CD8_share_of_CD8_plus_macrophages',family='capture_composition',value=n/(n+nm) if n+nm else np.nan,eligible_primary=bool(n>=20 and nm>=30),eligible_strict=bool(n>=100 and nm>=100),genes_expected=0,genes_measured=0,genes_detected=0,genes_detected_strict=0,reason='' if n>=20 and nm>=30 else 'CAPTURE_COUNT_GATE'))
  for i in idx:members.append(dict(dataset=ds,analysis_unit=unit,patient=str(m.iloc[i].patient),cell_id=str(m.index[i]),original_row=int(i)))
  fn=unit.replace(':','_').replace('/','_')+'.npz';np.savez_compressed(P/'inputs'/fn,genes=np.array(a.var_names),raw_sum=summ,detected_cells=det,available=np.array(sorted(available)),raw_rows=idx)
  cachemanifest.append(dict(dataset=ds,analysis_unit=unit,path=str(P/'inputs'/fn),cells=n))
 print('PREPARED',ds,flush=True);del a,m,x;gc.collect()
pd.DataFrame(out).to_csv(P/'tables/CD8_PATIENT_MEASURES.tsv',sep='\t',index=False);pd.DataFrame(cover).to_csv(P/'tables/CD8_GENE_COVERAGE.tsv.gz',sep='\t',index=False);pd.DataFrame(members).to_csv(P/'tables/CD8_CELL_MEMBERSHIP.tsv.gz',sep='\t',index=False);(P/'inputs/CACHE_MANIFEST.json').write_text(json.dumps(cachemanifest,indent=2));(P/'resources/INPUT_IDENTITIES.json').write_text(json.dumps(identities,indent=2))
(P/'qc/PREPARATION_COMPLETE.json').write_text(json.dumps({'status':'COMPLETE','receiver_rows':len(out),'source_rows':len(sources),'cells':len(members),'units':len(cachemanifest)},indent=2));print('PREPARATION_COMPLETE',len(out),len(sources),len(members),flush=True)
