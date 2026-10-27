from pathlib import Path
import pandas as pd,json,subprocess,h5py,numpy as np
from anndata.io import read_elem
P=Path(__file__).resolve().parents[1];R=P.parent
subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],check=True)
features={};src=[]
for sheet,cl,label in [('MSC','FAP+ fibroblasts','FAP_CAF25'),('Myeloid Cells','MARCO+ Macrophage','SPP1_MAC25')]:
 d=pd.read_excel(P/'resources/author_subtype_markers.xlsx',sheet_name=sheet);d=d[d.cluster.eq(cl)].sort_values(['avg_logFC','gene'],ascending=[False,True],kind='stable');d['rank_avg_logFC']=np.arange(1,len(d)+1)
 d.to_csv(P/'resources'/f'{label}_FULL_ORIGINAL_MARKERS.tsv',sep='\t',index=False);top=d.head(25);features[label]=top.gene.tolist();print(label,top[['gene','avg_logFC','p_val','p_val_adj']].to_json(orient='records'))
 print('border',d.iloc[24:27][['gene','avg_logFC']].to_json(orient='records'))
features['FAP_CAF25_disjoint']=sorted(set(features['FAP_CAF25'])-set(features['SPP1_MAC25']));features['SPP1_MAC25_disjoint']=sorted(set(features['SPP1_MAC25'])-set(features['FAP_CAF25']));print('SHARED',set(features['FAP_CAF25'])&set(features['SPP1_MAC25']))
go=pd.read_csv(R/'phase315_CRC_SPP1_barrier_paper_reproduction_20261004_v1/tables/SOURCE_Figure_6l.tsv',sep='\t')
for label,term in [('ECM_org','GO:0030198'),('Collagen_org','GO:0030199'),('TGFb_response','GO:0071559')]:
 genes=sorted(set('/'.join(go.loc[go.ID.eq(term),'geneID']).split('/')));features[label+'_original']=genes;features[label+'_heldout']=sorted(set(genes)-set(features['FAP_CAF25'])-set(features['SPP1_MAC25']));print(label,'original',len(genes),'heldout',len(features[label+'_heldout']))
(P/'resources/FEATURES_FROM_PAPER.json').write_text(json.dumps(features,indent=2))
reg=pd.read_csv(R/'phase234_part4_replan_execution_20260921_v1/tables/INPUT_SECTION_REGISTER.tsv',sep='\t');out=[]
for a in reg.itertuples():
 with h5py.File(a.input_path,'r') as h:
  genes=np.asarray(read_elem(h['var/gene_symbol_original'])) if 'gene_symbol_original' in h['var'] else np.asarray(read_elem(h['var'][h['var'].attrs['_index']]))
  layers=list(h.get('layers',{}));ids=np.asarray(read_elem(h['obs'][h['obs'].attrs['_index']])).astype(str)
  for nm,gs in features.items():out.append({'sample':a.sample,'dataset':a.dataset,'patient':a.patient,'feature':nm,'present':sum(g in set(genes) for g in gs),'total':len(gs),'missing':';'.join(g for g in gs if g not in set(genes))})
  for g in ['CD3D','CD8A']:out.append({'sample':a.sample,'dataset':a.dataset,'patient':a.patient,'feature':g,'present':int(g in genes),'total':1,'missing':'' if g in genes else g})
  if a.sample in ['GLMF1','GSM5924030','H46_T']:print('H5_SCHEMA',a.sample,'layers',layers,'X',dict(h['X'].attrs),'obsfields',list(h['obs'])[:12])
pd.DataFrame(out).to_csv(P/'qc/FEATURE_COVERAGE_PREFLIGHT.tsv',sep='\t',index=False);print('COVERAGE',pd.DataFrame(out).groupby('feature').agg(min_present=('present','min'),max_present=('present','max'),total=('total','first')).to_json())
print('PATIENT_PREP',reg[['sample','dataset','patient','preparation','treatment']].to_json(orient='records'))
