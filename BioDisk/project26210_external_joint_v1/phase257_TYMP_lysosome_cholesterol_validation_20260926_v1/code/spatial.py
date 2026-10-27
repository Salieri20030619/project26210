from common import *
import h5py
from anndata.io import read_elem
A=ROOT/'phase234_part4_replan_execution_20260921_v1';PRE=ROOT/'phase256_spatial_region_function_20260926_v1'
FUNC=['Lysosome','Cholesterol metabolism']
SPP=['V1_SPP1_q05','V3_SPP1_macro_share','SPP1_RNA']
def main():
 gate('spatial_start')
 reg=pd.read_csv(PRE/'inputs/INPUT_SECTION_REGISTER.tsv',sep='\t').fillna('')
 assert len(reg)==38 and reg[['dataset','patient']].drop_duplicates().shape[0]==37
 out('inputs/SPATIAL_INPUT_REGISTER.tsv',reg)
 manifest=[];rows=[];scorechecks=[];coverage=[];leaveblock=[]
 endpoints=[('TYMP_function','TYMP_RNA',t) for t in FUNC]
 endpoints += [('function_SPP1',f,t) for f in FUNC for t in SPP]
 endpoints += [('TYMP_SPP1','TYMP_RNA',t) for t in SPP]
 out('inputs/SPATIAL_FROZEN_ENDPOINTS.tsv',pd.DataFrame(endpoints,columns=['family','x','y']))
 for rec in reg.itertuples():
  gate('spatial_'+rec.sample)
  path=A/'inputs'/f'{rec.sample}_INPUTS.tsv.gz'
  gmt=PRE/'gsea'/f'{rec.dataset}__absolute__exclude_definition_genes__KEGG.gmt'
  src=[path,Path(rec.input_path),gmt,PRE/'effects'/f'{rec.sample}.npz'];before=[fingerprint(f) for f in src];manifest+=before
  d=pd.read_csv(path,sep='\t');allsets={f[0]:f[2:] for f in [l.split('\t') for l in gmt.read_text().splitlines()]}
  sets={t:allsets[t] for t in FUNC};shared=set(sets[FUNC[0]])&set(sets[FUNC[1]])
  with h5py.File(rec.input_path) as h:
   obs=read_elem(h['obs']);var=read_elem(h['var']);assert np.array_equal(obs.index.astype(str),d.spot_id.astype(str))
   X=read_elem(h['layers/counts'] if 'counts' in h.get('layers',{}) else h['X']).tocsr().astype(np.float64)
  assert X.shape[0]==len(d) and np.isfinite(X.data).all() and (X.data>=0).all()
  totals=np.asarray(X.sum(1)).ravel();np.testing.assert_allclose(np.log1p(totals),d.log1p_UMI,atol=1e-10)
  need=sorted(set(sum(sets.values(),[]))|{'TYMP','SPP1'});ix=var.index.get_indexer(need);assert (ix>=0).all()
  log=np.log1p(X[:,ix].toarray()/totals[:,None]*10000);lookup={g:i for i,g in enumerate(need)}
  sd=log.std(0);Z=(log-log.mean(0))/np.where(sd>1e-10,sd,1)
  for g in ['TYMP','SPP1']:np.testing.assert_allclose(log[:,lookup[g]],d[g+'_RNA'],rtol=1e-9,atol=1e-10)
  previous=np.load(PRE/'effects'/f'{rec.sample}.npz',allow_pickle=False);pi={g:i for i,g in enumerate(previous['genes'])}
  for f in FUNC:
   genes=sets[f];d[f]=Z[:,[lookup[g] for g in genes]].mean(1)
   reduced=[g for g in genes if g not in shared]
   d[f+'__no_shared']=Z[:,[lookup[g] for g in reduced]].mean(1) if len(reduced)>=10 and len(reduced)>=len(genes)*.2 else np.nan
   threshold=d.V1_SPP1_q05.quantile(.8);high=d.V1_SPP1_q05>=threshold
   delta=float(d.loc[high,f].mean()-d.loc[~high,f].mean())
   expected=float(previous['effects'][0,[pi[g] for g in genes]].mean());assert abs(delta-expected)<1e-10
   scorechecks.append(dict(sample=rec.sample,function=f,expected=expected,actual=delta,error=abs(delta-expected)))
   coverage.append(dict(sample=rec.sample,dataset=rec.dataset,function=f,genes=len(genes),members=';'.join(genes),shared_genes=';'.join(sorted(shared)),no_shared_genes=len(reduced),no_shared_members=';'.join(reduced)))
  columns=['spot_id','block_4','TYMP_RNA']+SPP+FUNC+[f+'__no_shared' for f in FUNC]+['log1p_UMI','total_macro_fraction','other_macro_fraction','malignant_fraction','endothelial_fraction','pericyte_fraction']
  d[columns].to_csv(P/'scores'/f'{rec.sample}_SCORES.tsv.gz',sep='\t',index=False)
  blocks=d.groupby('block_4')[columns[2:]].median().reset_index()
  for scale,frame,adj,version,role in [('whole',d,'adjusted','primary','primary'),('whole',d,'raw','primary','sensitivity'),('block4',blocks,'raw','primary','sensitivity'),('block4',blocks,'adjusted','primary','sensitivity'),('whole',d,'adjusted','no_shared','sensitivity')]:
   gate('associations_'+rec.sample+'_'+scale+'_'+adj+'_'+version)
   for fam,x,y in endpoints:
    xx=x+'__no_shared' if version=='no_shared' and x in FUNC else x
    yy=y+'__no_shared' if version=='no_shared' and y in FUNC else y
    # No duplicated sensitivity contrast for the direct TYMP-SPP1 endpoint.
    if version=='no_shared' and fam=='TYMP_SPP1':continue
    cs=['log1p_UMI','other_macro_fraction' if y=='V1_SPP1_q05' else 'total_macro_fraction','malignant_fraction','endothelial_fraction','pericyte_fraction'] if adj=='adjusted' else []
    info,_=correlation(frame[xx],frame[yy],{c:frame[c].to_numpy() for c in cs},min_n=50 if scale=='whole' else 8)
    rows.append(dict(sample=rec.sample,dataset=rec.dataset,patient=rec.patient,family=fam,x=x,y=y,scale=scale,adjustment=adj,version=version,role=role,**info))
    if scale=='whole' and adj=='adjusted' and version=='primary' and fam=='TYMP_function':
     br=[]
     for block in sorted(d.block_4.unique()):
      sub=d[d.block_4!=block];r,_=correlation(sub[x],sub[y],{c:sub[c].to_numpy() for c in cs},min_n=50);br.append(r['rho'])
     br=np.asarray(br);valid=np.isfinite(br)
     leaveblock.append(dict(sample=rec.sample,dataset=rec.dataset,patient=rec.patient,function=y,full_rho=info['rho'],valid_leave_blocks=int(valid.sum()),same_sign=int((np.sign(br[valid])==np.sign(info['rho'])).sum()),min_rho=float(np.nanmin(br)) if valid.any() else np.nan,max_rho=float(np.nanmax(br)) if valid.any() else np.nan))
  assert before==[fingerprint(f) for f in src]
  print('SPATIAL_DONE',rec.sample,flush=True)
 sec=pd.DataFrame(rows);out('tables/SPATIAL_SECTION_EFFECTS.tsv',sec);out('qc/SPATIAL_SCORE_RECONSTRUCTION.tsv',pd.DataFrame(scorechecks));out('tables/SPATIAL_PROGRAM_COVERAGE.tsv',pd.DataFrame(coverage));out('tables/SPATIAL_LEAVE_BLOCK.tsv',pd.DataFrame(leaveblock))
 keys=['dataset','family','x','y','scale','adjustment','version','role'];patient=[]
 for key,g in sec.groupby(keys+['patient'],dropna=False):
  vals=g.rho.dropna().to_numpy();rho=float(np.tanh(np.mean(np.arctanh(np.clip(vals,-.999999999,.999999999))))) if len(vals) else np.nan
  patient.append(dict(zip(keys+['patient'],key),rho=rho,sections=len(vals),section_ids=';'.join(g.loc[g.rho.notna(),'sample'])))
 pat=pd.DataFrame(patient);out('tables/SPATIAL_PATIENT_EFFECTS.tsv',pat)
 result=[]
 for key,g in pat.groupby(keys,dropna=False):result.append(dict(zip(keys,key),**cohort_summary(g.rho,'|'.join(key))))
 sm=correction(pd.DataFrame(result),['dataset','family','role']);out('tables/SPATIAL_COHORT_SUMMARY.tsv',sm)
 jwrite('qc/SPATIAL_SOURCE_MANIFEST.json',manifest)
 jwrite('SPATIAL_COMPLETE.json',dict(status='PASS',time=stamp(),sections=38,registered_patients=37,comparisons=len(sec),cohort_results=len(sm),QR_checks=len(CHECKS),max_QR_error=max(CHECKS)))
if __name__=='__main__':main()
