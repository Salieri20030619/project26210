from common import *
def main():
 gate('specificity_start')
 assert sha(P/'PROTOCOL_BEFORE_TESTING.md')==json.loads((P/'FREEZE.json').read_text())['protocol_sha256']
 z=data();rows=[];grid=[];parts=[];counts=[]
 for patient,g in z.groupby('patient'):
  study=g.study.iloc[0]
  for scope,h in [('C4',g[g.candidate=='4']),('OtherMacro',g[g.candidate!='4'])]:
   for gene in GENES:
    nx=int(h.det_TYMP.sum());ny=int(h['det_'+gene].sum());eligible=len(h)>=50 and nx>=10 and ny>=10
    for model in ['core','context']:
     if eligible:e,n,why=correlation(h.TYMP,h[gene],covariates(h,model,scope=='OtherMacro'))
     else:e,n,why=np.nan,len(h),'CELL_OR_DETECTION_GATE'
     rows.append(dict(patient=patient,study=study,scope=scope,gene=gene,model=model,effect=e,cells=len(h),TYMP_detected=nx,target_detected=ny,reason=why))
   # Decomposition only TYMP-SPP1, exactly fixed gates.
   x=h.det_TYMP.to_numpy(bool);y=h.det_SPP1.to_numpy(bool);n11=int((x&y).sum());n10=int((x&~y).sum());n01=int((~x&y).sum());n00=int((~x&~y).sum())
   counts.append(dict(patient=patient,study=study,scope=scope,cells=len(h),n11=n11,n10=n10,n01=n01,n00=n00,raw_odds_ratio=n11*n00/(n10*n01) if n10*n01 else np.nan))
   base=len(h)>=50 and x.sum()>=10 and y.sum()>=10
   for mode in ['co_detection','double_positive_intensity']:
    hh=h if mode=='co_detection' else h[x&y]
    ok=base and ((len(h)-x.sum()>=10 and len(h)-y.sum()>=10) if mode=='co_detection' else len(hh)>=20)
    for model in ['core','context']:
     if ok:
      xx=hh.det_TYMP if mode=='co_detection' else hh.TYMP;yy=hh.det_SPP1 if mode=='co_detection' else hh.SPP1
      e,n,why=correlation(xx,yy,covariates(hh,model,scope=='OtherMacro'),min_n=50 if mode=='co_detection' else 20,binary=mode=='co_detection')
     else:e,n,why=np.nan,len(hh),'PART_SPECIFIC_CELL_OR_DETECTION_GATE'
     parts.append(dict(patient=patient,study=study,scope=scope,part=mode,model=model,effect=e,cells=len(hh),n_double_positive=n11,reason=why))
  for cluster,h in g.groupby('candidate'):
   eligible=len(h)>=50 and h.det_TYMP.sum()>=10 and h.det_SPP1.sum()>=10
   for model in ['core','context']:
    if eligible:e,n,why=correlation(h.TYMP,h.SPP1,covariates(h,model))
    else:e,n,why=np.nan,len(h),'CELL_OR_DETECTION_GATE'
    grid.append(dict(patient=patient,study=study,cluster=cluster,model=model,effect=e,cells=len(h),reason=why))
  print('SPECIFICITY_PATIENT',patient,flush=True)
 e=pd.DataFrame(rows);tab(e,'SCOPE_PATIENT_EFFECTS.tsv')
 ss=summarize(e,['scope','gene','model'],['scope'],'SCOPE_SUMMARY.tsv')
 pairs=e[e.scope=='C4'].merge(e[e.scope=='OtherMacro'],on=['patient','study','gene','model'],suffixes=('_C4','_Other'),validate='one_to_one')
 pairs['effect']=pairs.effect_C4-pairs.effect_Other;pairs['fisher_z_difference']=np.arctanh(pairs.effect_C4.clip(-.99999999,.99999999))-np.arctanh(pairs.effect_Other.clip(-.99999999,.99999999));pairs['family']='paired_specificity_8'
 tab(pairs,'C4_OTHER_PAIRED_PATIENT_EFFECTS.tsv')
 summarize(pairs,['family','gene','model'],['family'],'C4_OTHER_PAIRED_SUMMARY.tsv')
 study_descriptions(pairs,['gene','model'],'PAIRED_SPECIFICITY')
 gr=pd.DataFrame(grid)
 # Full 12-cluster x36-patient x2-model grid, absent clusters retained.
 idx=pd.MultiIndex.from_product([sorted(z.patient.unique()),sorted(z.candidate.unique()),['core','context']],names=['patient','cluster','model'])
 gr=gr.set_index(['patient','cluster','model']).reindex(idx).reset_index();gr['study']=gr.patient.map(z.groupby('patient').study.first());gr['reason']=gr.reason.fillna('NO_CELLS');gr['family']='state_grid_24'
 tab(gr,'STATE_GRID_PATIENT_EFFECTS.tsv');summarize(gr,['family','cluster','model'],['family'],'STATE_GRID_SUMMARY.tsv')
 pa=pd.DataFrame(parts);pa['family']='decomposition_8';tab(pa,'DECOMPOSITION_PATIENT_EFFECTS.tsv');tab(counts,'DETECTION_2X2_COUNTS.tsv')
 summarize(pa,['family','scope','part','model'],['family'],'DECOMPOSITION_SUMMARY.tsv');study_descriptions(pa,['scope','part','model'],'DECOMPOSITION')
 manifest=[dict(path=str(f),sha256=sha(f)) for f in [PRE/'inputs/CELL_CACHE.tsv.gz',PRE/'inputs/CONTEXT_CELL_CACHE.tsv.gz',M/'inputs/MACRO_METADATA.tsv.gz']]
 dump(P/'qc/SPECIFICITY_INPUTS.json',manifest);dump(P/'qc/SPECIFICITY_COMPLETE.json',dict(status='COMPUTED',scope_rows=len(e),paired_rows=len(pairs),grid_rows=len(gr),decomposition_rows=len(pa),figures=0))
 print('SPECIFICITY_COMPLETE',flush=True)
if __name__=='__main__':main()
