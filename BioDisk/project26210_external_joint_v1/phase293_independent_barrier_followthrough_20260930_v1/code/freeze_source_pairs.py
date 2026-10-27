from pathlib import Path
import pandas as pd,json,hashlib,datetime
P=Path(__file__).resolve().parents[1]
d=pd.read_csv(P/'tables/EXTERNAL_LIBRARY_INPUTS.tsv',sep='\t')
panels=pd.read_csv(P/'tables/FROZEN_FUNCTION_GENES.tsv',sep='\t').program.unique()
rows=[]
for src in ['tumor_like_TYMP','macro_TYMP','macro_SPP1']:
 source_count='tumor_like_cells' if src=='tumor_like_TYMP' else 'macro_cells';minn=50 if src=='tumor_like_TYMP' else 100
 ok=d[source_count].ge(minn)&((d.receiver.eq('strict_fibro')&d.receiver_cells.ge(20))|(d.receiver.eq('old_fibro_pool')&d.receiver_cells.ge(50)))
 for (ds,pat,layer,recv),g in d[ok].groupby(['dataset','patient','stratum','receiver']):
  for e in panels:
   z=g[g[e].notna()&g[src].notna()]
   if len(z):
    rows.append(dict(dataset=ds,patient=pat,stratum=layer,receiver=recv,source=src,endpoint=e,libraries=';'.join(z.library),n_libraries=len(z),x=z[src].mean(),y=z[e].mean(),receiver_depth=z.receiver_depth.mean(),all_libraries_triple_qualified=bool(z.qualified.all())))
out=pd.DataFrame(rows)
out.to_csv(P/'tables/EXTERNAL_SOURCEPAIR_FREEZE.tsv',sep='\t',index=False)
files=['EXTERNAL_IDENTITY_AND_STATISTICS_PROTOCOL.md','SOURCE_PAIR_ELIGIBILITY_ADDENDUM.md','tables/EXTERNAL_SOURCEPAIR_FREEZE.tsv','tables/EXTERNAL_PATIENT_FREEZE.tsv']
(P/'qc/SOURCEPAIR_FREEZE.json').write_text(json.dumps(dict(at=datetime.datetime.now().astimezone().isoformat(),outcomes_examined=False,files={n:hashlib.sha256((P/n).read_bytes()).hexdigest() for n in files}),indent=2))
print(out.groupby(['dataset','stratum','receiver','source']).patient.nunique().to_string())
