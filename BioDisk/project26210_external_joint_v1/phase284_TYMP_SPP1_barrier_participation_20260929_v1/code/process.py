from common import *
NEW=['RAB37','MMP9','MMP14','TGFB1','PDGFB','VEGFA']
LINKS=[('TYMP',g) for g in NEW]+[('RAB37','SPP1'),('CTSB','SPP1')]
CTX=['Interferon Alpha Response','Interferon Gamma Response','Hypoxia']
def main():
 gate('process_start')
 assert sha(P/'PROTOCOL_BEFORE_TESTING.md')==json.loads((P/'FREEZE.json').read_text())['protocol_sha256']
 paths=[M/'inputs/MACRO_METADATA.tsv.gz',M/'inputs/MACRO_COUNTS.npz',M/'inputs/GENES.tsv',PRE/'resources/CONTEXT_PROGRAMS.json']
 manifest=[dict(path=str(f),sha256=sha(f)) for f in paths]
 m=read(paths[0],index_col=0,dtype={'candidate':str});X=sparse.load_npz(paths[1]).tocsr();genes=read(paths[2]).gene.tolist();gi={g:i for i,g in enumerate(genes)}
 total=np.asarray(X.sum(1)).ravel();assert np.allclose(total,m.library_counts) and len(m)==22444 and m.candidate.eq('4').sum()==2674
 ctx=json.loads(paths[3].read_text())['genes'];excluded=set(NEW+['TYMP','SPP1','CTSB'])
 ctx={k:[g for g in v if g not in excluded and g in gi] for k,v in ctx.items()}
 dump(P/'resources/PROCESS_CONTEXT_PROGRAMS.json',dict(programs=ctx,excluded=sorted(excluded),single_genes=NEW,ordered_links=LINKS))
 union=sorted(set(sum(list(ctx.values()),[])+NEW+['TYMP','SPP1','CTSB','LGALS3']))
 assert all(g in gi for g in union);raw=X[:,[gi[g] for g in union]].toarray().astype(float);log=np.log1p(raw/total[:,None]*1e4);ix={g:i for i,g in enumerate(union)}
 rows=[];cache=[];qualification=[]
 for patient,g in m[m.source_tissue.str.lower().eq('tumor')].groupby('curated_patient_id'):
  ids=m.index.get_indexer(g.index);rv=raw[ids];lv=log[ids];vals={gene:lv[:,ix[gene]] for gene in NEW+['TYMP','SPP1','CTSB','LGALS3']}
  for name,panel in ctx.items():
   v=lv[:,[ix[x] for x in panel]];sd=v.std(0);vals[name]=((v-v.mean(0))/np.where(sd>0,sd,1)).mean(1)
  dep=np.log1p(total[ids]);ng=g.n_genes_raw.to_numpy()
  cache.append(pd.DataFrame(dict(cell_id=g.index,patient=patient,study=g.study.iloc[0],candidate=g.candidate.to_numpy(),log_depth=dep,n_genes=ng,**vals,**{'det_'+gene:rv[:,ix[gene]]>0 for gene in NEW+['TYMP','SPP1','CTSB']})))
  for scope,mask,minn in [('C4',g.candidate.eq('4').to_numpy(),50),('WholeMacro',np.ones(len(g),bool),100)]:
   cov=np.column_stack([dep,ng])
   if scope=='WholeMacro':cov=np.column_stack([cov,pd.get_dummies(g.candidate,dtype=float,drop_first=True).to_numpy()])
   for source,target in LINKS:
    nx=int(((rv[:,ix[source]]>0)&mask).sum());ny=int(((rv[:,ix[target]]>0)&mask).sum());ok=mask.sum()>=minn and nx>=10 and ny>=10
    qualification.append(dict(patient=patient,study=g.study.iloc[0],scope=scope,source=source,target=target,cells=int(mask.sum()),source_positive=nx,target_positive=ny,eligible=ok))
    for model in ['core','context']:
     cc=cov if model=='core' else np.column_stack([cov]+[vals[n] for n in CTX])
     if ok:e,n,reason=fit(vals[source],vals[target],cc,mask,min_n=minn)
     else:e,n,reason=np.nan,int(mask.sum()),'CELL_OR_DETECTION_GATE'
     rows.append(dict(patient=patient,study=g.study.iloc[0],scope=scope,source=source,target=target,model=model,effect=e,cells=n,reason=reason))
  print('PROCESS_PATIENT',patient,flush=True)
 e=pd.DataFrame(rows);tab(e,'PROCESS_PATIENT_EFFECTS.tsv');tab(qualification,'PROCESS_QUALIFICATION.tsv')
 pd.concat(cache,ignore_index=True).to_csv(P/'inputs/PROCESS_CELL_CACHE.tsv.gz',sep='\t',index=False)
 summary(e,['scope','source','target','model'],['scope'],'PROCESS_SUMMARY.tsv')
 describe_study(e,['scope','source','target','model'],'PROCESS')
 # Inherited molecular links stay in original family and are not re-labelled as independent validation.
 old=read(PRE/'tables/WITHIN_MACRO_SUMMARY.tsv');old=old[old.endpoint.isin(['SPP1','CTSB','LGALS3'])].copy();old['inheritance']='Phase280 original effects and family';tab(old,'INHERITED_TYMP_PROCESS_LINKS.tsv')
 # Exact same-library malignant-to-macrophage correspondence.
 src=R/'phase157_ordered_bridge_reassessment_20260915_v1/02_inputs';groupfile=src/'GROUP_INDEX.tsv';arrayfile=src/'GROUP_PSEUDOBULK.npz';basefile=R/'phase154_same_library_TYMP_C4_test_20260915_v1/SAME_LIBRARY_PATIENT_METRICS.tsv'
 G=read(groupfile);zz=np.load(arrayfile,allow_pickle=True);gg=list(zz['genes']);lookup={g:i for i,g in enumerate(gg)};C=zz['counts'];base=read(basefile);out=[]
 assert len(base)==16 and base.patient.is_unique
 for r in base.itertuples():
  o=dict(patient=r.patient,study=r.study,library=r.sample_id)
  for compartment,roles in [('malignant',['Malignant']),('macro',['C4','OtherMacrophage'])]:
   gr=G[(G.patient==r.patient)&(G.sample_id==r.sample_id)&G.role.isin(roles)]
   assert gr.source_tissue.eq('Tumor').all() and set(gr.study)=={r.study}
   v=C[gr.group_index.to_numpy(int)].sum(0);n=int(gr.cells.sum());tot=float(v.sum());lp=np.log2(v/tot*1e6+1)
   expected_n=r.malignant_cells if compartment=='malignant' else r.macrophage_cells;expected_tot=r.malignant_total_counts if compartment=='malignant' else r.macrophage_total_counts
   assert n==expected_n and tot==expected_tot
   o[compartment+'_TYMP']=float(lp[lookup['TYMP']]);o[compartment+'_cells']=n;o[compartment+'_total_counts']=tot
   if compartment=='macro':
    for gene in NEW:o[gene]=float(lp[lookup[gene]])
   else:assert abs(o['malignant_TYMP']-r.TYMP_pseudobulk)<1e-10
  out.append(o)
 d=pd.DataFrame(out);tab(d,'MALIGNANT_TO_MACRO_PATIENT_INPUTS.tsv');results=[];leave=[]
 for target in NEW:
  for model in ['study','study_macroTYMP']:
   extra=d.macro_TYMP.to_numpy() if model=='study_macroTYMP' else None
   result=cross(d.malignant_TYMP,d[target],d.study,extra=extra,seed=262128400+len(results))
   results.append(dict(target=target,model=model,**result,family='malignant_to_macro_12'))
   for st in sorted(d.study.unique()):
    s=d[d.study!=st];q=cross(s.malignant_TYMP,s[target],s.study,extra=s.macro_TYMP.to_numpy() if extra is not None else None,B=0)
    leave.append(dict(target=target,model=model,omitted_study=st,n=q['n'],rho=q['rho']))
 tab(bh(pd.DataFrame(results),['family']),'MALIGNANT_TO_MACRO_SUMMARY.tsv');tab(leave,'MALIGNANT_TO_MACRO_LEAVE_STUDY.tsv')
 for f in [groupfile,arrayfile,basefile]:manifest.append(dict(path=str(f),sha256=sha(f)))
 for rec in manifest:assert sha(rec['path'])==rec['sha256']
 dump(P/'qc/PROCESS_INPUTS.json',manifest);dump(P/'qc/PROCESS_COMPLETE.json',dict(status='COMPUTED',patient_effects=len(e),same_library_patients=len(d),new_cross_endpoints=len(results)))
 print('PROCESS_COMPLETE',flush=True)
if __name__=='__main__':main()
