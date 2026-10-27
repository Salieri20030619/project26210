from pathlib import Path
import json,hashlib,datetime,traceback,importlib.util
import numpy as np,pandas as pd,h5py
from scipy import sparse
from anndata.io import read_elem
P=Path(__file__).resolve().parents[1];R=P.parent
spec=importlib.util.spec_from_file_location('source',P/'code/source_analysis.py');b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
def main():
 b.gate('receiver_extract');cfg=json.loads((P/'PROGRAMS_FROZEN.json').read_text());outputs=cfg['outputs']
 inter=pd.read_csv(P/'resources/CELLCHAT_OUTPUT_INTERACTIONS.tsv',sep='\t').fillna('');comp=pd.read_csv(P/'resources/CELLCHAT_COMPLEXES.tsv',sep='\t').fillna('')
 subcols=[x for x in comp if x.startswith('subunit')];complexes={r.complex_name:[str(r[c]).strip() for c in subcols if str(r[c]).strip()] for _,r in comp.iterrows()}
 candidates=[]
 for r in inter.itertuples():
  subs=complexes.get(r.receptor,[r.receptor]);assert all('_' not in g for g in subs)
  candidates.append(dict(interaction=r.interaction_name,ligand=r.ligand,receptor=r.receptor,subunits=';'.join(subs),annotation=r.annotation,evidence=r.evidence))
 b.tab(candidates,'FROZEN_RECEPTOR_CANDIDATES.tsv');genes=sorted(set(g for c in candidates for g in c['subunits'].split(';')))
 b.tab([dict(gene=g,database_interactions=int((inter.ligand==g).sum()),status='PRESENT' if (inter.ligand==g).any() else 'NO_LOCKED_DB_ENTRY_NOT_BIOLOGICAL_ABSENCE') for g in outputs],'DATABASE_COVERAGE.tsv')
 fmeta=R/'phase211_immune_receiver_refreeze_20260919_v1/freeze_v2/DISCOVERY_RECEIVER_CELL_FREEZE.tsv.gz'
 cd=pd.read_csv(fmeta,sep='\t');cd=cd[(cd.receiver_label_v1=='CD8_T')&cd.source_tissue.str.lower().eq('tumor')]
 cd=cd[['cell_id','study','curated_patient_id','sample_id']].rename(columns={'curated_patient_id':'patient','sample_id':'library'});cd['receiver']='CD8_T'
 ff=R/'phase219_spatial_reference_identifiability_20260920_v1/reference/ALL_CELL_REFERENCE_DISPOSITION.tsv.gz'
 fib=pd.read_csv(ff,sep='\t');fib=fib[(fib.label=='Fibroblast')&fib.tissue.str.lower().eq('tumor')][['cell_id','study','patient']];fib['receiver']='Fibroblast';fib['library']=''
 selected=pd.concat([cd,fib],ignore_index=True).set_index('cell_id');assert selected.index.is_unique
 src=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/03_malignancy_cnv/cnv_annotated_h5ad_v1/part1_scvi_integrated_cell_annotated_cnv_v1.h5ad')
 st=src.stat();result=np.zeros((len(selected),len(genes)));total=np.zeros(len(selected));seen=np.zeros(len(selected),bool)
 with h5py.File(src,'r',rdcc_nbytes=128*1024*1024) as f:
  ids=f['obs/cell_id'].asstr()[:];var=f['var/gene'].asstr()[:];libraries=np.asarray(read_elem(f['obs/sample_id'])).astype(str)
  positions=pd.Index(ids).get_indexer(selected.index);assert (positions>=0).all()
  selected.loc[selected.receiver=='Fibroblast','library']=libraries[positions[selected.receiver.to_numpy()=='Fibroblast']]
  gidx=pd.Index(var).get_indexer(genes);present=gidx>=0
  x=f['layers/counts'];ptr=x['indptr'][:];lookup={int(v):i for i,v in enumerate(positions)}
  for start in range(0,len(ids),4096):
   end=min(start+4096,len(ids));loc=[j for j in range(start,end) if j in lookup]
   if not loc:continue
   lo,hi=int(ptr[start]),int(ptr[end]);mat=sparse.csr_matrix((x['data'][lo:hi],x['indices'][lo:hi],ptr[start:end+1]-lo),shape=(end-start,len(var)))
   a=mat[np.array(loc)-start];rows=np.array([lookup[j] for j in loc]);total[rows]=np.asarray(a.sum(1,dtype=np.float64)).ravel()
   result[np.ix_(rows,np.flatnonzero(present))]=a[:,gidx[present]].toarray();seen[rows]=True
 assert seen.all() and (total>0).all() and np.isfinite(result).all() and (result>=0).all()
 selected.reset_index().to_csv(P/'inputs/RECEIVER_CELL_METADATA.tsv.gz',sep='\t',index=False)
 np.savez_compressed(P/'inputs/RECEIVER_COUNTS.npz',cell_id=selected.index.to_numpy(str),genes=np.array(genes),gene_present=present,counts=result,total=total)
 log=np.log1p(result/total[:,None]*10000);gmap={g:i for i,g in enumerate(genes)};rows=[];coverage=[]
 sc=pd.read_csv(P/'inputs/CELL_SCORE_CACHE.tsv.gz',sep='\t').set_index('cell_id');mm=pd.read_csv(R/'phase200_union_C4_C10_C15_reassessment_20260918_v1/inputs/MACRO_METADATA.tsv.gz',sep='\t',index_col=0,dtype={'candidate':str})
 sc['library']=mm.loc[sc.index,'sample_id'];sc=sc[sc.candidate.astype(str)=='4']
 for key,grp in selected.groupby(['receiver','study','patient','library'],sort=True):
  ii=selected.index.get_indexer(grp.index);n=len(ii);receiver,study,patient,library=key;sender=sc[(sc.patient==patient)&(sc.library==library)]
  coverage.append(dict(receiver=receiver,study=study,patient=patient,library=library,receiver_cells=n,C4_cells=len(sender)))
  for c in candidates:
   jj=[gmap[g] for g in c['subunits'].split(';')];all_present=bool(present[jj].all());joint=(result[np.ix_(ii,jj)]>0).all(1) if all_present else np.zeros(n,bool);positive=int(joint.sum());ok=n>=20 and all_present and positive>=5 and positive/n>=.1
   sender_positive=int((sender[c['ligand']]>0).sum())
   rows.append(dict(**c,receiver=receiver,study=study,patient=patient,library=library,receiver_cells=n,complete_genes_present=all_present,joint_positive_cells=positive,joint_fraction=positive/n,receptor_eligible=ok,min_subunit_mean_log=float(log[np.ix_(ii,jj)].mean(0).min()) if all_present else np.nan,C4_sender_cells=len(sender),same_library_C4_ge20=len(sender)>=20,C4_ligand_detected_cells=sender_positive,C4_ligand_detected_fraction=sender_positive/len(sender) if len(sender) else np.nan))
 rr=pd.DataFrame(rows);b.tab(rr,'RECEIVER_COMPLEX_LIBRARY_QUALIFICATION.tsv');b.tab(coverage,'RECEIVER_CELL_LIBRARY_COUNTS.tsv')
 summaries=[]
 for key,g in rr.groupby(['ligand','receptor','receiver']):
  good=g[g.receptor_eligible];matched=good[good.same_library_C4_ge20];detected=matched[matched.C4_ligand_detected_cells>0]
  summaries.append(dict(zip(['ligand','receptor','receiver'],key),total_patients=g.patient.nunique(),receiver_eligible_patients=good.patient.nunique(),same_library_C4_patients=matched.patient.nunique(),same_library_C4_ligand_detected_patients=detected.patient.nunique(),eligible_libraries=len(good),inference='DESCRIPTIVE_ONLY'))
 b.tab(summaries,'RECEIVER_PATIENT_COVERAGE.tsv')
 assert src.stat().st_size==st.st_size and src.stat().st_mtime_ns==st.st_mtime_ns
 b.dump(P/'qc/RECEIVER_SOURCE_MANIFEST.json',dict(path=str(src),bytes=st.st_size,mtime_ns=st.st_mtime_ns,metadata_sources=[dict(path=str(p),sha256=b.sha(p)) for p in [fmeta,ff]],resource_sha256=b.sha(R/'phase35_part1_c01_c04_cell_communication_v1/resource_lock_v1/CellChatDB.human_v1.rds')))
 b.dump(P/'qc/RECEIVER_COMPUTATION_COMPLETED.json',dict(at=b.now(),cells=len(selected),genes=len(genes),interactions=len(candidates),qualification_rows=len(rr),library_groups=len(coverage)))
 b.status('COMPUTED_AWAITING_VALIDATION');print('RECEIVER_COMPLETE',len(selected),len(rr),flush=True)
if __name__=='__main__':
 try:main()
 except Exception:b.status('FAILED',traceback=traceback.format_exc());raise
