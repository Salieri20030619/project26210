import pysam,collections,pathlib,json
p=pathlib.Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase240_bulk_immune_analysis_20260924_v1')
fa=pysam.FastaFile(str(p.parent/'phase236_public_bulk_conversion_20260923_v1/reference/GRCh38.primary_assembly.genome.fa'))
c=collections.Counter();examples=[]
for line in (p/'reference/SNP_GRCh38_hg38_wChr.bed').read_text().splitlines()[:300]:
 chrom,s,e,rs,a,b=line.split()[:6];s=int(s);e=int(e);c['rows']+=1;c['length_'+str(e-s)]+=1
 for delta in [-1,0,1]:
  base=fa.fetch(chrom,s+delta,s+delta+1).upper()
  c['delta'+str(delta)+'_ref']+=int(base==a)
  c['delta'+str(delta)+'_either']+=int(base in [a,b])
 if len(examples)<8:examples.append([line,fa.fetch(chrom,s-1,e+1)])
print(json.dumps({'counts':dict(c),'examples':examples}))
