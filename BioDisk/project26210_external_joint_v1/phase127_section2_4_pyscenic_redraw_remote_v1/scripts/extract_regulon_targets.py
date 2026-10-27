from pathlib import Path
import pickle,csv,json
w=Path(__file__).resolve().parents[1]
base=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/07_pyscenic_v1/intermediate')
rows=[];stats=[]
for method in ['GENIE3','GRNBoost2']:
 p=base/('regulons_'+method.lower()+'_v1.dat')
 with p.open('rb') as f:regs=pickle.load(f)
 for reg in regs:
  if reg.name in ['ILF2(+)','HOXD3(+)','HNF4A(+)']:
   stats.append({'method':method,'regulon':reg.name,'targets':len(reg.gene2weight)})
   for target,weight in reg.gene2weight.items():rows.append([method,reg.name[:-3],reg.name,target,float(weight)])
assert len(stats)==6
with (w/'data/CANDIDATE_REGULON_TARGETS.tsv').open('w') as f:
 wr=csv.writer(f,delimiter='\t');wr.writerow(['method','tf','regulon','target','weight']);wr.writerows(rows)
(w/'qc/REGULON_EXTRACTION.json').write_text(json.dumps(stats,indent=2))
print(json.dumps(stats))
