from common import *
def main():
 gate('regulation_start')
 net=R/'part2_rebuild_v1/08_regulons/network_v2_20260915/REGULON_TARGET_EDGES.tsv'
 ed=read(net);focus=ed[ed.target.isin(['TYMP','SPP1'])];names=sorted(focus.regulon.unique());assert len(names)==14
 ctx=json.loads((PRE/'resources/CONTEXT_PROGRAMS.json').read_text())['genes']
 exclusion=set(['TYMP','SPP1','CTSB','CTSD','LGALS3']+sum(list(ctx.values()),[]))
 meta=read(M/'inputs/MACRO_METADATA.tsv.gz',index_col=0,dtype={'candidate':str});genes=read(M/'inputs/GENES.tsv').gene.tolist();gi={g:i for i,g in enumerate(genes)}
 sets={};qual=[]
 for name in names:
  orig=set(ed.loc[ed.regulon==name,'target']);remain=sorted(orig-exclusion);measured=[g for g in remain if g in gi]
  ok=len(measured)>=10 and len(measured)>=.8*len(remain);sets[name]=measured
  qual.append(dict(regulon=name,TF=ed.loc[ed.regulon==name,'TF'].iloc[0],network_anchor=';'.join(sorted(focus.loc[focus.regulon==name,'target'])),original_targets=len(orig),removed=len(orig&exclusion),remaining=len(remain),measured=len(measured),eligible=ok,genes=';'.join(measured)))
 assert all(not set(v)&exclusion for v in sets.values())
 tab(qual,'REGULON_GENE_QUALIFICATION.tsv');tab(focus,'LOCKED_NETWORK_TYMP_SPP1_EDGES.tsv')
 dump(P/'resources/REGULON_PROGRAMS_FROZEN.json',dict(source=str(net),source_sha256=sha(net),programs=sets,excluded=sorted(exclusion),score='per-patient mean of standardized log1pCP10K target expression; not AUCell/ATAC activity'))
 z=data();X=sparse.load_npz(M/'inputs/MACRO_COUNTS.npz').tocsr();assert len(meta)==22444
 total=np.asarray(X.sum(1)).ravel();assert np.allclose(total,meta.library_counts)
 union=sorted(set(sum(list(sets.values()),[])));cols=[gi[g] for g in union];look={g:i for i,g in enumerate(union)};rr=[];condition=[];cache=[];detcache=[];elig=[]
 globally={r['regulon']:r['eligible'] for r in qual}
 for patient,h in z.groupby('patient'):
  ids=meta.index.get_indexer(h.cell_id);assert (ids>=0).all();raw=X[ids][:,cols].toarray().astype(float);log=np.log1p(raw/total[ids,None]*1e4)
  values={};dets={}
  for name,gs in sets.items():
   jj=[look[g] for g in gs];v=log[:,jj];sd=v.std(0)
   values[name]=((v-v.mean(0))/np.where(sd>0,sd,1)).mean(1) if globally[name] else np.full(len(h),np.nan)
   dets[name]=(raw[:,jj]>0).sum(1)>=2
  cache.append(pd.DataFrame(dict(cell_id=h.cell_id.to_numpy(),patient=patient,**values)))
  detcache.append(pd.DataFrame(dict(cell_id=h.cell_id.to_numpy(),patient=patient,**dets)))
  for scope,mask in [('C4',h.candidate.eq('4').to_numpy()),('OtherMacro',h.candidate.ne('4').to_numpy())]:
   hh=h.loc[mask];n=len(hh)
   for name in names:
    sc=values[name][mask];nd=int(dets[name][mask].sum());progok=globally[name] and nd>=20 and n>=50
    elig.append(dict(patient=patient,study=h.study.iloc[0],scope=scope,regulon=name,cells=n,program_detected=nd,eligible=progok))
    for gene in ['TYMP','SPP1']:
     ok=progok and int(hh['det_'+gene].sum())>=10
     for model in ['core','context']:
      if ok:e,nn,why=correlation(sc,hh[gene],covariates(hh,model,scope=='OtherMacro'))
      else:e,nn,why=np.nan,n,'PROGRAM_OR_CELL_DETECTION_GATE'
      rr.append(dict(patient=patient,study=h.study.iloc[0],scope=scope,regulon=name,gene=gene,model=model,effect=e,cells=nn,reason=why))
    ok=progok and hh.det_TYMP.sum()>=10 and hh.det_SPP1.sum()>=10
    if ok:
     base,nn,_=correlation(hh.TYMP,hh.SPP1,covariates(hh,'context',scope=='OtherMacro'))
     effect,nn,why=correlation(hh.TYMP,hh.SPP1,covariates(hh,'context',scope=='OtherMacro',sc))
    else:base,effect,nn,why=np.nan,np.nan,n,'PROGRAM_OR_CELL_DETECTION_GATE'
    condition.append(dict(patient=patient,study=h.study.iloc[0],scope=scope,regulon=name,effect=effect,baseline_rho=base,rho_change=effect-base,cells=nn,reason=why))
  print('REGULON_PATIENT',patient,flush=True)
 pd.concat(cache,ignore_index=True).to_csv(P/'inputs/REGULON_CELL_SCORES.tsv.gz',sep='\t',index=False)
 pd.concat(detcache,ignore_index=True).to_csv(P/'inputs/REGULON_CELL_DETECTION.tsv.gz',sep='\t',index=False)
 e=pd.DataFrame(rr);tab(e,'REGULON_PATIENT_EFFECTS.tsv');tab(elig,'REGULON_PATIENT_QUALIFICATION.tsv')
 summarize(e,['scope','regulon','gene','model'],['scope'],'REGULON_ASSOCIATION_SUMMARY.tsv')
 study_descriptions(e,['scope','regulon','gene','model'],'REGULON_ASSOCIATIONS')
 c=pd.DataFrame(condition);c['family']='conditional_remaining_28';tab(c,'CONDITIONAL_REMAINING_PATIENT_EFFECTS.tsv')
 ss=summarize(c,['family','scope','regulon'],['family'],'CONDITIONAL_REMAINING_SUMMARY.tsv')
 delta=c.groupby(['scope','regulon']).agg(n_complete=('rho_change',lambda v:int(v.notna().sum())),median_baseline=('baseline_rho','median'),median_remaining=('effect','median'),median_within_patient_change=('rho_change','median'),lower_remaining=('rho_change',lambda v:int((v<0).sum())),higher_remaining=('rho_change',lambda v:int((v>0).sum()))).reset_index()
 tab(delta,'CONDITIONAL_CHANGE_DESCRIPTIVE.tsv')
 manifest=[dict(path=str(f),sha256=sha(f)) for f in [net,M/'inputs/MACRO_COUNTS.npz',M/'inputs/MACRO_METADATA.tsv.gz',M/'inputs/GENES.tsv',PRE/'resources/CONTEXT_PROGRAMS.json']]
 dump(P/'qc/REGULATION_INPUTS.json',manifest)
 dump(P/'qc/REGULATION_COMPLETE.json',dict(status='COMPUTED',regulons=len(names),qualified_regulons=sum(globally.values()),patient_effects=len(e),conditional_effects=len(c),figures=0))
 print('REGULATION_COMPLETE',flush=True)
if __name__=='__main__':main()
