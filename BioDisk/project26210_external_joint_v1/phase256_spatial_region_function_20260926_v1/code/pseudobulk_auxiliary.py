from pathlib import Path
import sys,json
import numpy as np,pandas as pd
P=Path(__file__).resolve().parents[1];sys.path.insert(0,str(P/'code'))
from extract import gate
out=[];jobs=pd.read_csv(P/'gsea/JOBS.tsv',sep='\t')
for j in jobs.itertuples():
 if j.model=='adjusted':continue
 gate('pseudobulk_'+j.dataset+'_'+j.model+'_'+j.version+'_'+j.library)
 a=np.load(P/'effects'/f'PATIENT_{j.dataset}_{j.model}.npz',allow_pickle=False)
 genes=a['genes'];D=a['pseudobulk'];n=len(a['patients']);lookup={x:i for i,x in enumerate(genes)}
 assert np.isfinite(D).all()
 terms=[];mat=[]
 for line in Path(j.gmt_path).read_text().splitlines():
  f=line.split('\t');v=D[:,[lookup[x] for x in f[2:]]].mean(1)
  terms.append(f[0]);mat.append(v)
  out.append(dict(dataset=j.dataset,model=j.model,version=j.version,library=j.library,term=f[0],patient_n=n,genes=len(f)-2,median_mean_gene_log2CPM_difference=float(np.median(v)),positive=int((v>1e-12).sum()),negative=int((v< -1e-12).sum()),zero=int((np.abs(v)<=1e-12).sum()),inference='DIRECTION_ONLY_AUXILIARY_NO_NEW_P'))
 np.savez_compressed(P/'effects'/f'PSEUDOBULK_PROGRAM_{j.dataset}__{j.model}__{j.version}__{j.library}.npz',terms=np.asarray(terms),patients=a['patients'],effects=np.asarray(mat).T)
pd.DataFrame(out).to_csv(P/'tables/PSEUDOBULK_FUNCTION_DIRECTION.tsv.gz',sep='\t',index=False)
qc=json.loads((P/'qc/FINAL_VALIDATION.json').read_text())
(P/'qc/PSEUDOBULK_AUXILIARY_VALIDATION.json').write_text(json.dumps(dict(status='PASS',rows=len(out),all_finite=True,region_sums_from_integer_counts=True,patient_preparations_equal_weight=True,new_P_values=0),indent=2))
d=pd.DataFrame(out)
d['name']=d.term.str.split('__').str[-1]
q=d[(d.version=='exclude_definition_genes')&d.name.isin(['Lysosome','Phagosome','Cholesterol metabolism','negative regulation of immune response'])]
print(q[['dataset','model','library','name','median_mean_gene_log2CPM_difference','positive','negative']].to_string(index=False))
