from common import *
def main():
 z=data();rows=[]
 for patient,g in z.groupby('patient'):
  for scope,h in [('C4',g[g.candidate=='4']),('OtherMacro',g[g.candidate!='4'])]:
   ok=len(h)>=50 and h.det_TYMP.sum()>=10 and h.det_SPP1.sum()>=10 and (~h.det_SPP1).sum()>=10 and h.TYMP.nunique()>=3
   for model in ['core','context']:
    if ok:e,n,why=correlation(stats.rankdata(h.TYMP),h.det_SPP1,covariates(h,model,scope=='OtherMacro'),binary=True)
    else:e,n,why=np.nan,len(h),'CELL_OR_SPP1_DETECTION_GATE'
    rows.append(dict(patient=patient,study=g.study.iloc[0],scope=scope,model=model,effect=e,cells=len(h),TYMP_detection_fraction=float(h.det_TYMP.mean()) if len(h) else np.nan,SPP1_detection_fraction=float(h.det_SPP1.mean()) if len(h) else np.nan,reason=why,family='continuous_TYMP_to_SPP1_detection_4'))
 d=pd.DataFrame(rows);tab(d,'TYMP_INTENSITY_SPP1_DETECTION_PATIENT_EFFECTS.tsv')
 summarize(d,['family','scope','model'],['family'],'TYMP_INTENSITY_SPP1_DETECTION_SUMMARY.tsv')
 study_descriptions(d,['scope','model'],'TYMP_INTENSITY_SPP1_DETECTION')
 # Program overlap is descriptive; no new hypothesis test.
 cfg=json.loads((P/'resources/REGULON_PROGRAMS_FROZEN.json').read_text())['programs'];out=[]
 for a in cfg:
  for bb in cfg:
   aa=set(cfg[a]);bset=set(cfg[bb]);inter=len(aa&bset)
   out.append(dict(program_a=a,program_b=bb,a_genes=len(aa),b_genes=len(bset),shared=inter,Jaccard=inter/len(aa|bset),fraction_a_shared=inter/len(aa)))
 tab(out,'REGULON_PROGRAM_OVERLAP.tsv')
 dump(P/'qc/DETECTION_EXTENSION_COMPLETE.json',dict(status='COMPUTED',patient_rows=len(d),new_hypotheses=4))
 print('DETECTION_EXTENSION_COMPLETE',flush=True)
if __name__=='__main__':main()
