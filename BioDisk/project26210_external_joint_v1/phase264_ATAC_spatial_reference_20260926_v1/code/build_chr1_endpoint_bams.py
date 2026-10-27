from pathlib import Path
import json,collections,bisect,hashlib
import numpy as np,pandas as pd,pysam
P=Path(__file__).resolve().parents[1];D=P/'downloads/RCC100';O=P/'atac_chr1_pilot';O.mkdir(exist_ok=True)
assign=pd.read_csv(P/'tables/RCC100_ATAC_LINEAGE_ASSIGNMENTS.tsv',sep='\t')
qc=pd.read_csv(P/'tables/RCC100_LINEAGE_QC.tsv',sep='\t')
eligible=qc.loc[qc.n_eligible_for_pseudobulk,'lineage'].tolist()
assert len(eligible)>=2,'Fewer than two lineages meet the predefined cell-count gate'
lookup=assign[assign.lineage.isin(eligible)].set_index('barcode').lineage.to_dict()
frag=D/'fragments.tsv.gz';assert frag.is_file()
if not Path(str(frag)+'.tbi').exists():pysam.tabix_index(str(frag),preset='bed',force=False)
length=248956422
header={'HD':{'VN':'1.6','SO':'unsorted'},'SQ':[{'SN':'chr1','LN':length}],
'CO':['Coordinate-only compatibility records derived from 10x fragments; not original sequenced read alignments.',
'One-base CIGAR encodes required end positions only; SEQ/QUAL absent. Original mate identities and MAPQ are not reconstructed.',
'10x shifts reversed (+4/-5); HINT cut positions and original fragment lengths retained. Each unique fragment counted once; PCR support counts not replicated.']}
handles={};unsorted={}
for group in eligible:
 path=O/f'{group}_endpoints.unsorted.bam';unsorted[group]=path
 handles[group]=pysam.AlignmentFile(path,'wb',header=header)
counts=collections.Counter();bcseen=set();excluded=collections.Counter()
centers=np.linspace(100000,length-100000,100,dtype=int)
windows=np.column_stack([centers-250,centers+250])
direct={g:np.zeros((100,2),dtype=np.int64) for g in eligible}
source_entries=0
with pysam.TabixFile(str(frag)) as tb:
 assert 'chr1' in tb.contigs
 for serial,line in enumerate(tb.fetch('chr1')):
  fields=line.split('\t');start,end=int(fields[1]),int(fields[2]);barcode=fields[3]
  source_entries+=1
  group=lookup.get(barcode)
  if group is None:continue
  bcseen.add(barcode)
  left=start-4;right_end=end+5;insert=right_end-left
  if left<0 or right_end>length or insert<=0:excluded['invalid_end_position']+=1;continue
  name='fragment_'+str(serial)
  r1=pysam.AlignedSegment();r1.query_name=name;r1.flag=99;r1.reference_id=0;r1.reference_start=left;r1.mapping_quality=0;r1.cigar=[(0,1)];r1.next_reference_id=0;r1.next_reference_start=right_end-1;r1.template_length=insert;r1.set_tag('CB',barcode);r1.set_tag('ZT','fragment_endpoint')
  r2=pysam.AlignedSegment();r2.query_name=name;r2.flag=147;r2.reference_id=0;r2.reference_start=right_end-1;r2.mapping_quality=0;r2.cigar=[(0,1)];r2.next_reference_id=0;r2.next_reference_start=left;r2.template_length=-insert;r2.set_tag('CB',barcode);r2.set_tag('ZT','fragment_endpoint')
  assert r1.reference_start+4==start and r2.reference_end-5-1==end-1
  assert abs(r1.template_length)==end-start+9
  handles[group].write(r1);handles[group].write(r2);counts[group]+=1
  for strand,cut in [(0,start),(1,end-1)]:
   j=np.searchsorted(windows[:,0],cut,side='right')-1
   if j>=0 and cut<windows[j,1]:direct[group][j,strand]+=1
for h in handles.values():h.close()
results=[]
for group in eligible:
 target=O/f'{group}_chr1_fragment_endpoints.bam'
 pysam.sort('-@','2','-m','512M','-o',str(target),str(unsorted[group]))
 pysam.index(str(target))
 observed=np.zeros((100,2),dtype=np.int64)
 with pysam.AlignmentFile(target,'rb') as b:
  assert b.mapped==2*counts[group],(group,b.mapped,counts[group])
  for j,(left,right) in enumerate(windows):
   for r in b.fetch('chr1',max(0,left-10),min(length,right+10)):
    cut=r.reference_end-6 if r.is_reverse else r.reference_start+4
    if left<=cut<right:observed[j,int(r.is_reverse)]+=1
 assert np.array_equal(observed,direct[group]),group
 results.append(dict(lineage=group,n_cells=int((assign.lineage==group).sum()),unique_chr1_fragments=counts[group],coordinate_records=2*counts[group],window_check=200,window_max_abs_error=int(np.abs(observed-direct[group]).max()),eligible_for_footprinting=bool(counts[group]>=20000),path=str(target)))
pd.DataFrame(results).to_csv(P/'tables/CHR1_ENDPOINT_INPUT_QC.tsv',sep='\t',index=False)
(P/'qc/ENDPOINT_COMPATIBILITY_CHECK.json').write_text(json.dumps(dict(status='PASS',source_fragments=str(frag),source_chr1_entries=source_entries,accepted_fragments=dict(counts),excluded=dict(excluded),cell_barcodes_with_fragments=len(bcseen),checks=results,raw_sequenced_BAM=False,read_sequence_reconstructed=False,forward_cut='fragment start',reverse_cut='fragment end - 1',original_insert_length='fragment end - fragment start + 9',mapq='unknown represented by 0; published fragment filtering retained'),indent=2))
print('ENDPOINT_CHECK_PASS',json.dumps(results),flush=True)

