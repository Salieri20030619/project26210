from common import *
from scipy import sparse
import h5py
from anndata.io import read_elem
def main():
 gate('state_alignment')
 old=R/'phase287_stromal_refinement_image_morphology_20260929_v1'
 m=read(old/'tables/STROMAL_ANNOTATION_V1.tsv.gz');raw=sparse.load_npz(old/'inputs/STROMAL_COUNTS.npz');genes=np.load(old/'inputs/GENES.npy').tolist()
 panels=json.loads((P/'PANELS_FROZEN.json').read_text());union=sorted(set(sum(panels.values(),[])));present=[g for g in union if g in genes]
 jj=[genes.index(g) for g in present];x=raw[:,jj].toarray().astype(float);ix={g:i for i,g in enumerate(present)}
 norm=np.log1p(x/m.total_counts_raw.to_numpy()[:,None]*1e4)
 cov=[]
 for name,gs in panels.items():
  kept=[g for g in gs if g in ix];ok=len(kept)>=.8*len(gs)
  cov.append(dict(program=name,present=len(kept),requested=len(gs),missing=';'.join(g for g in gs if g not in ix),eligible=ok))
  m['program_'+name]=norm[:,[ix[g] for g in kept]].mean(1) if ok else np.nan
  m['det_'+name]=(x[:,[ix[g] for g in kept]]>0).sum(1) if ok else np.nan
 bad=m.competing_immune|m.competing_endo|m.competing_epi
 malignant=m.malignant_conservative_v1.astype(str).str.lower().isin(['true','1'])
 pool=m.cell_type_broad_v1.isin(['Fibroblast','Pericyte/SMC'])&~bad&~malignant
 m['stromal_pool']=pool
 m['activated_candidate']=pool&m.det_CAF_activation.ge(2)&(m.fibro_candidate|m.cluster_primary.isin([9,15]))
 m['transition_candidate']=pool&m.det_transition_RCC.ge(2)&m.det_perivascular.ge(1)
 m['overlap_candidate']=m.activated_candidate&m.transition_candidate
 m.to_csv(P/'tables/STROMAL_STATE_MEMBERSHIP.tsv.gz',sep='\t',index=False)
 np.savez_compressed(P/'inputs/STROMAL_PANEL_COUNTS.npz',cell_id=m.cell_id.to_numpy(str),genes=np.array(present),counts=x,total_counts=m.total_counts_raw.to_numpy())
 tab(cov,'STATE_PANEL_COVERAGE.tsv')
 rows=[];pb=[];pvec=[];sizes=[]
 for (st,pat,lib,tissue),g in m.groupby(['study','curated_patient_id','sample_id','source_tissue']):
  rec=dict(study=st,patient=pat,library=lib,tissue=tissue,selected=len(g),stromal_pool=int(g.stromal_pool.sum()),activated_candidate=int(g.activated_candidate.sum()),transition_candidate=int(g.transition_candidate.sum()),overlap_candidate=int(g.overlap_candidate.sum()))
  rows.append(rec)
  for label,mask in [('stromal_pool',g.stromal_pool),('activated',g.activated_candidate),('transition',g.transition_candidate),('conservative',g.conservative_fibro_reference)]:
   ii=g.index[mask].to_numpy()
   if not len(ii):continue
   counts=x[ii].sum(0);total=float(m.total_counts_raw.iloc[ii].sum());log=np.log1p(counts/total*1e4)
   rec2=dict(study=st,patient=pat,library=lib,tissue=tissue,population=label,cells=len(ii),total_UMI=total)
   for name,gs in panels.items():
    inds=[ix[k] for k in gs if k in ix];rec2[name]=float(log[inds].mean()) if len(inds)>=.8*len(gs) else np.nan
   for gene in ['CD44','ITGAV','ITGB1','TGFBR1','TGFBR2','TYMP','SPP1']:
    rec2[gene+'_positive_cells']=int((x[ii,ix[gene]]>0).sum());rec2[gene+'_logCP10K']=float(log[ix[gene]])
   pb.append(rec2);pvec.append(counts);sizes.append(total)
 tab(rows,'STATE_LIBRARY_CAPACITY.tsv');tab(pb,'STROMAL_STATE_PSEUDOBULK.tsv')
 np.savez_compressed(P/'inputs/STROMAL_STATE_PSEUDOBULK_COUNTS.npz',genes=np.array(present),counts=np.array(pvec),total_counts=np.array(sizes))
 desc=[]
 for (cl,st,tissue),g in m.groupby(['cluster_primary','study','source_tissue']):
  for name in panels:
   desc.append(dict(cluster=int(cl),study=st,tissue=tissue,program=name,cells=len(g),patients=g.curated_patient_id.nunique(),mean_score=float(g['program_'+name].mean()),median_detected=float(g['det_'+name].median())))
 tab(desc,'CLUSTER_STATE_DESCRIPTIVE.tsv')
 # Compartment source evidence in exact frozen mature macrophages and original H5 cell types.
 sourceH=H
 with h5py.File(sourceH,'r') as h:
  fields=['cell_type_broad_v1','study','curated_patient_id','sample_id','source_tissue','analysis_keep','malignant_conservative_v1']
  orig=pd.DataFrame({k:np.asarray(read_elem(h['obs'][k])) for k in fields});orig['cell_id']=h['obs/cell_id'].asstr()[:]
  gg=h['var/gene'].asstr()[:].tolist();gi=[gg.index(k) for k in ['TYMP','SPP1']]
  xx=h['layers/counts'];ptr=xx['indptr'][:];out=np.zeros((len(orig),2));tot=np.zeros(len(orig))
  for start in range(0,len(orig),4096):
   end=min(start+4096,len(orig));lo,hi=ptr[start],ptr[end]
   block=sparse.csr_matrix((xx['data'][lo:hi],xx['indices'][lo:hi],ptr[start:end+1]-lo),shape=(end-start,len(gg)))
   out[start:end]=block[:,gi].toarray();tot[start:end]=np.asarray(block.sum(1)).ravel()
 mm=read(R/'phase200_union_C4_C10_C15_reassessment_20260918_v1/inputs/MACRO_METADATA.tsv.gz',dtype={'candidate':str})
 # Inspect schema explicitly, match original cell ids, never infer by order.
 idcol='cell_id' if 'cell_id' in mm else 'cell'
 assert idcol in mm.columns
 matched=pd.Index(orig.cell_id).get_indexer(mm[idcol]);assert (matched>=0).all()
 orig['source_group']=orig.cell_type_broad_v1
 orig.loc[matched,'source_group']='Frozen_other_mature_macrophage'
 orig.loc[matched[mm.candidate.eq('4')],'source_group']='Frozen_C4'
 pos=orig.malignant_conservative_v1.astype(str).str.lower().isin(['true','1']);orig.loc[pos,'source_group']='Conservative_malignant'
 stm=pd.Index(orig.cell_id).get_indexer(m.loc[m.stromal_pool,'cell_id']);assert (stm>=0).all();orig.loc[stm,'source_group']='Refined_stromal_pool'
 # These labels are only source-audit columns; original labels/C4 file are untouched.
 sr=[]
 for (st,pat,lib,tissue,label),g in orig.groupby(['study','curated_patient_id','sample_id','source_tissue','source_group']):
  ii=g.index.to_numpy()
  for k,gene in enumerate(['TYMP','SPP1']):
   sr.append(dict(study=st,patient=pat,library=lib,tissue=tissue,source_group=label,cells=len(ii),gene=gene,positive_cells=int((out[ii,k]>0).sum()),raw_counts=float(out[ii,k].sum()),total_UMI=float(tot[ii].sum()),pseudobulk_logCP10K=float(np.log1p(out[ii,k].sum()/tot[ii].sum()*1e4))))
 tab(sr,'SOURCE_GENE_PATIENT_LIBRARY.tsv')
 dump(P/'qc/STATE_ALIGNMENT_COMPLETE.json',dict(status='COMPUTED',cells=len(m),stromal_pool=int(pool.sum()),activated=int(m.activated_candidate.sum()),transition=int(m.transition_candidate.sum()),overlap=int(m.overlap_candidate.sum()),original_C4_members_unchanged=True,new_figures=0,source_cells=len(orig),frozen_macro_cells=len(mm)))
 print('STATE_ALIGNMENT_COMPLETE',int(pool.sum()),int(m.activated_candidate.sum()),int(m.transition_candidate.sum()),flush=True)
if __name__=='__main__':main()
