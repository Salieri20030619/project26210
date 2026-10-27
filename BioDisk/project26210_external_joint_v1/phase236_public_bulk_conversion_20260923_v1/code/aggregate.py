from common import *
import pandas as pd,numpy as np,csv,gzip
samples=json.loads((P/'SAMPLES.json').read_text());ann=pd.read_csv(P/'tables/gene_annotation.tsv',sep='\t',index_col='gene_id',keep_default_na=False);frames={};qcs=[]
for s in samples:
 d=P/'samples'/s['run'];q=json.loads((d/'QA.json').read_text());assert q['status']=='PASS';qcs.append(dict(s,**q))
 x=pd.read_csv(d/'gene_counts.tsv',sep='\t',header=None,names=['gene_id','count'],index_col=0);x=x.loc[~x.index.str.startswith('__'),'count'];assert x.index.is_unique;assert set(x.index)==set(ann.index)
 assert (x>=0).all() and np.equal(x,np.floor(x)).all();frames[s['run']]=x
counts=pd.DataFrame(frames).reindex(ann.index).astype('int64');length=ann['exon_union_bp'].astype(float)
rpk=counts.div(length/1000,axis=0);assert (rpk.sum()>0).all();tpm=rpk.div(rpk.sum(),axis=1)*1e6
den=pd.Series({q['run']:q['mapped_fragments'] for q in qcs});assert (den>0).all();fpkm=counts.div(length,axis=0).div(den,axis=1)*1e9
assert np.allclose(tpm.sum(),1e6,rtol=1e-10) and np.isfinite(tpm.to_numpy()).all() and np.isfinite(fpkm.to_numpy()).all()
outputs=[]
for group in ['all_libraries_technical_only']+sorted({s['group'] for s in samples}):
 selected=[s['run'] for s in samples if group=='all_libraries_technical_only' or s['group']==group]
 out=P/'matrices'/group;out.mkdir(exist_ok=True)
 for name,table in [('gene_counts',counts),('TPM',tpm),('FPKM',fpkm)]:
  path=out/(name+'.tsv.gz');table[selected].to_csv(path,sep='\t',compression='gzip');outputs.append(dict(path=str(path),genes=len(table),libraries=len(selected),sha256=sha(path)))
 pd.DataFrame([s for s in samples if s['run'] in selected]).to_csv(out/'samples.tsv',sep='\t',index=False)
pd.DataFrame(qcs).to_csv(P/'tables/sample_QC_summary.tsv',sep='\t',index=False)
qa=dict(status='PASS',sample_libraries=60,independent_donors='not inferred',genes=len(ann),counts_nonnegative_integers=True,TPM_sums_one_million=True,FPKM_definition='mapped paired fragments and union exon length',all_samples_counted_after_empirical_strand_inference=True,outputs=outputs,biological_comparisons_started=False,donor_fingerprint_analysis_started=False,finished_at=now())
save(P/'MATRIX_QA.json',qa)
