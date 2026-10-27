from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='4'
import json,subprocess
import numpy as np,pandas as pd,anndata as ad
P=Path(__file__).resolve().parents[1]
proto=json.loads((P/'resources/QUALIFICATION_PROTOCOL.json').read_text());pan=proto['identity_panels']
maps={'GSE210038':{'Macrophage':[17],'Fibro_candidate':[2],'Mural':[1,4],'RCC_like':[16,21],'Tcell':[3,6,7,8,9,10,11,12,18]},'GSE224630':{'Macrophage':[16],'Fibro_candidate':[18],'Mural':[5,6,7],'RCC_like':[0,1,2,3,8,10,11,13,14,15,19],'Tcell':[12]}}
freeze=P/'resources/LINEAGE_RULES_FROZEN.json'
rules=dict(cluster_maps=maps,author='GSE222703 Celltype_Harmony: Mac_C1QA, Fibro, ccRCC, CD8T',fibro='Reviewed stromal cluster or author Fibro; >=2 core & >=1 support; <2 immune/endo/epithelial; exclude >=2 pericyte OR (MYH11/CNN1 and ACTA2/TAGLN). Overlap retained as separate transition pool.',CD8='T cluster and >=1 CD3D/CD3E/TRAC and >=1 CD8A/CD8B; author CD8T for GSE222703',gates={'Macrophage':100,'Fibro_strict':20,'Fibro_transition':20,'Mural':50,'RCC_like':50,'CD8':20},associations_inspected=False)
if not freeze.exists():freeze.write_text(json.dumps(rules,indent=2))
else:assert json.loads(freeze.read_text())==rules
counts=[];gene_rows=[]
for ds in ['GSE210038','GSE224630','GSE222703']:
 subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],check=True,stdout=subprocess.DEVNULL)
 a=ad.read_h5ad(P/f'inputs/{ds}_counts_qc.h5ad');x=a.X.tocsr()
 def ng(gs):return np.asarray((x[:,a.var_names.isin(gs)]>0).sum(axis=1)).ravel()
 masks={}
 if ds in maps:
  for k,cl in maps[ds].items():masks[k]=a.obs.cluster.astype(str).isin([str(c) for c in cl]).to_numpy()
  masks['CD8']=masks['Tcell']&(ng(['CD3D','CD3E','TRAC'])>=1)&(ng(['CD8A','CD8B'])>=1)
 else:
  for k,lab in [('Macrophage','Mac_C1QA'),('Fibro_candidate','Fibro'),('RCC_like','ccRCC'),('CD8','CD8T')]:masks[k]=a.obs.Celltype_Harmony.eq(lab).to_numpy()
  masks['Mural']=np.zeros(a.n_obs,dtype=bool)
 fibro_ok=(ng(pan['fibro_core'])>=2)&(ng(pan['fibro_support'])>=1)
 other=(ng(pan['immune'])>=2)|(ng(pan['endothelial'])>=2)|(ng(pan['epithelial'])>=2)
 mural=(ng(pan['pericyte'])>=2)|((ng(['MYH11','CNN1'])>=1)&(ng(['ACTA2','TAGLN'])>=1))
 # Include reviewed mural clusters in the identity audit: a true fibro-like subset can pass the same gate.
 stromal=masks['Fibro_candidate']|masks['Mural']
 masks['Fibro_strict']=stromal&fibro_ok&~other&~mural
 masks['Fibro_transition']=stromal&fibro_ok&~other&mural
 for k,m in masks.items():a.obs['lineage_'+k]=m
 a.obs.to_csv(P/f'tables/{ds}_QUALIFIED_CELL_METADATA.tsv.gz',sep='\t')
 for sample,g in a.obs.groupby('sample',observed=True):
  ii=a.obs['sample'].eq(sample).to_numpy();rr=dict(dataset=ds,sample=sample,patient=str(g.patient.iloc[0]),main=bool(g.main.iloc[0]),cells=int(ii.sum()))
  for k,m in masks.items():rr[k]=int((ii&m).sum())
  counts.append(rr)
  for k in ['Fibro_candidate','Fibro_strict','Fibro_transition','Mural','Macrophage','RCC_like','CD8']:
   sel=ii&masks[k]
   for gene in sorted(set(sum(pan.values(),[]))):
    if gene not in a.var_names or not sel.any():continue
    j=a.var_names.get_loc(gene);v=x[sel,j].toarray().ravel()
    gene_rows.append(dict(dataset=ds,sample=sample,lineage=k,n=int(sel.sum()),gene=gene,detected_fraction=float((v>0).mean())))
 print(ds,'DONE',flush=True)
pd.DataFrame(counts).fillna(0).to_csv(P/'tables/LINEAGE_COUNTS_BY_PATIENT.tsv',sep='\t',index=False)
pd.DataFrame(gene_rows).to_csv(P/'tables/LINEAGE_IDENTITY_GENE_COVERAGE.tsv',sep='\t',index=False)
print(pd.DataFrame(counts).fillna(0).to_string(index=False),flush=True)
(P/'qc/LINEAGE_QUALIFICATION_COMPLETE.json').write_text(json.dumps(dict(status='COMPLETE',association_tests_run=False),indent=2))
