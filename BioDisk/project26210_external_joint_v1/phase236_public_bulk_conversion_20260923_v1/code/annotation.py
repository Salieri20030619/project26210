from common import *
import re,collections,csv

patch=json.loads((P/'software/HISAT2_PLATFORM_PATCH.json').read_text())
source_file=P/'reference/SOURCE_DOWNLOADS.json'
sources=json.loads(source_file.read_text())
for source in sources:
 if source['path'].endswith('/hisat2-2.1.0.zip'):
  source['path']=patch['official_archive'];source['derived_archive']=patch['derived_archive'];source['derived_sha256']=patch['derived_sha256'];source['platform_patch']=str(P/'software/HISAT2_PLATFORM_PATCH.json')
save(source_file,sources)
ref=json.loads((P/'reference/REFERENCE_FROZEN.json').read_text())
contigs={line.split('\t')[0] for line in pathlib.Path(ref['fasta']+'.fai').read_text().splitlines()}
genes={};tx={};skipped=set()
def merged(intervals):
 out=[]
 for a,b in sorted(set(intervals)):
  if out and a<=out[-1][1]:out[-1][1]=max(out[-1][1],b)
  else:out.append([a,b])
 return out
with open(ref['gtf']) as f:
 for line in f:
  if line.startswith('#'):continue
  v=line.rstrip('\n').split('\t')
  if len(v)!=9 or v[2]!='exon':continue
  if v[0] not in contigs:skipped.add(v[0]);continue
  a=dict(re.findall(r'(\w+) "([^"]*)"',v[8]));g=a['gene_id'];t=a['transcript_id'];start=int(v[3])-1;end=int(v[4]);assert start<end
  rec=genes.setdefault(g,dict(gene_id=g,gene_symbol=a.get('gene_name',''),gene_type=a.get('gene_type',''),chromosome=v[0],strand=v[6],intervals=[]));assert rec['chromosome']==v[0] and rec['strand']==v[6];rec['intervals'].append((start,end))
  tr=tx.setdefault(t,dict(chromosome=v[0],strand=v[6],intervals=[]));tr['intervals'].append((start,end))
assert not skipped,'GTF contigs absent from genome: '+str(skipped)
with (P/'tables/gene_annotation.tsv').open('w') as f:
 w=csv.DictWriter(f,['gene_id','gene_symbol','gene_type','chromosome','strand','exon_union_bp'],delimiter='\t');w.writeheader()
 for g,rec in sorted(genes.items()):
  x=dict(rec);intervals=x.pop('intervals');x['exon_union_bp']=sum(b-a for a,b in merged(intervals));assert x['exon_union_bp']>0;w.writerow(x)
with (P/'reference/transcripts.bed12').open('w') as f:
 for name,rec in sorted(tx.items()):
  iv=merged(rec['intervals']);start=iv[0][0];end=iv[-1][1]
  f.write('\t'.join(map(str,[rec['chromosome'],start,end,name,0,rec['strand'],start,end,0,len(iv),','.join(str(b-a) for a,b in iv),','.join(str(a-start) for a,b in iv)]))+'\n')
save(P/'reference/ANNOTATION_QA.json',dict(genes=len(genes),transcripts=len(tx),all_contigs_in_reference=True,gene_annotation_sha256=sha(P/'tables/gene_annotation.tsv'),bed_sha256=sha(P/'reference/transcripts.bed12')))
