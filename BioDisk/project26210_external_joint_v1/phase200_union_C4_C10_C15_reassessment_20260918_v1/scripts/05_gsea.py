from common import *
import json
models=['all_U','all_C4','all_U_common','all_C4_common','all_U_without_C10_dominant_patient']
reg=[]
for model in models:
 z=np.load(D/'core/contrasts'/f'{model}.npz',allow_pickle=True)
 e=z['effects'];sd=e.std(0,ddof=1);ok=z['tested']&(sd>0)&np.isfinite(sd)
 score=e[:,ok].mean(0)/(sd[ok]/np.sqrt(len(e)))
 a=pd.DataFrame({'gene':z['genes'][ok],'score':score}).sort_values(['score','gene'],ascending=[False,True],kind='stable')
 a.to_csv(D/'gsea'/f'{model}_rank.tsv',sep='\t',index=False)
 reg.append(dict(model=model,n=len(e),genes=len(a),ties=int(a.score.duplicated().sum())))
pd.DataFrame(reg).to_csv(D/'gsea/RANK_REGISTRY.tsv',sep='\t',index=False)
print('RANKS_READY')
