import pathlib,json,hashlib,importlib.util,sys
import numpy as np
import pandas as pd
import pysam
P=pathlib.Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('phase240_pipeline',P/'code/pipeline.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
bed=P/'reference/SNP_GRCh38_hg38_wChr.bed'
b=bed.read_bytes()
tree=json.loads((P/'reference/SOURCE_TREE.json').read_text())
expected=next(x['sha'] for x in tree['tree'] if x['path']=='SNP/SNP_GRCh38_hg38_wChr.bed')
assert hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()==expected
assert hashlib.sha256(b).hexdigest()==m.CFG['reference_sha256']
assert np.allclose(m.rank(np.array([1,2,2,4])),[1,2.5,2.5,4])
assert np.allclose(m.bh([.01,.04,.03]),[.03,.04,.04])
assert abs(m.spear(np.array([1,2,3]),np.array([3,2,1]))+1)<1e-12

m.CFG['statistics']['permutations']=199
m.CFG['statistics']['bootstrap']=100
r=m.evaluate(np.arange(12.),np.arange(12.),np.random.default_rng(123))
assert r['rho']>.999 and r['p']<=.01 and r['ci_low']>.999
qc=pd.read_csv(m.CONV/'tables/sample_QC_summary.tsv',sep='\t')
assert len(qc)==60 and qc.run.is_unique
for col in ['run','sample_alias','group','assigned_fragments','BAM']: assert col in qc
ref=list((m.CONV/'reference').glob('*.fa'))[0]
panel=[]
with pysam.FastaFile(str(ref)) as fa:
 for line in b.decode().splitlines():
  c,s,e,rs,a,alt=line.split()[:6];s=int(s);e=int(e)
  if c in ['chr'+str(i) for i in range(1,23)] and e==s+1 and a in 'ACGT' and alt in 'ACGT' and a!=alt:
   base=fa.fetch(c,s,e).upper()
   if base not in [a,alt]:
    comp=str.maketrans('ACGT','TGCA');a,alt=a.translate(comp),alt.translate(comp)
   if base not in [a,alt]:continue
   if base==alt:a,alt=alt,a
   panel.append((c,s,a,alt,rs))
assert len(panel)>10000
row=qc.iloc[0]
with pysam.AlignmentFile(row.BAM,'rb') as bam:
 assert bam.has_index()
 for c,s,a,alt,rs in panel[:20]:
  for col in bam.pileup(c,s,s+1,truncate=True,stepper='samtools',min_base_quality=20,min_mapping_quality=30,max_depth=2000,ignore_overlaps=True,ignore_orphans=True,compute_baq=False,flag_filter=3844):
   assert col.reference_pos==s
result={'status':'PASS','reference_git_blob_verified':True,'reference_matched_autosomal_rows':len(panel),'rank_ties_test':True,'BH_test':True,'Spearman_permutation_bootstrap_test':True,'BAM_index_pileup_smoke':True,'libraries':len(qc),'versions':{'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__,'pysam':pysam.__version__}}
m.dump(P/'PREFLIGHT_VALIDATION.json',result)
print(json.dumps(result),flush=True)
