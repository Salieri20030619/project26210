from pathlib import Path
import pickle,pandas as pd,json,hashlib
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/part2_rebuild_v1');O=R/'08_regulons/network_v2_20260915';assert (O/'RUN_COMPLETE').exists()
with (O/'regulons.dat').open('rb') as f:reg=pickle.load(f)
rows=[];edges=[]
for r in reg:
 rows.append({'regulon':r.name,'TF':r.name.split('(')[0],'n_targets':len(r.gene2weight),'targets':';'.join(sorted(r.gene2weight))})
 for g,w in r.gene2weight.items():edges.append({'regulon':r.name,'TF':r.name.split('(')[0],'target':g,'weight':float(w)})
pd.DataFrame(rows).to_csv(O/'REGULON_REGISTRY.tsv',sep='\t',index=False);pd.DataFrame(edges).to_csv(O/'REGULON_TARGET_EDGES.tsv',sep='\t',index=False)
refs=Path('/mnt/BioDisk/ludingyu/project26210/references/pyscenic/hg38');records=[]
for p in refs.iterdir():
 if p.is_file():
  h=hashlib.sha256()
  with p.open('rb') as f:
   for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
  records.append({'path':str(p),'bytes':p.stat().st_size,'sha256':h.hexdigest()})
(O/'REFERENCE_HASHES.json').write_text(json.dumps(records,indent=2));(O/'REGISTRY_STATUS.json').write_text(json.dumps({'status':'COMPLETE_NEW_REGULON_NETWORK','regulons':len(rows),'edges':len(edges),'figures_created':False},indent=2));print('REGULONS',len(rows),'EDGES',len(edges))
