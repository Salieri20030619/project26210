from source import *
import tarfile,gzip,re,io
def main():
 d=P/'resources/GSE171542';f=d/'GSE171542_RAW_proxy.tar'
 assert sha(f)=='c42963b5bba65559bc7153a0c43b441823965ab5eaa41b0bb778fda517229849'
 soft=(d/'GSE171542_samples.soft').read_text();meta={}
 for block in soft.split('^SAMPLE = ')[1:]:
  gsm=block.splitlines()[0].strip();title=re.search(r'^!Sample_title = (.+)$',block,re.M).group(1)
  tr=re.search(r'^!Sample_characteristics_ch1 = treatment: (.+)$',block,re.M).group(1)
  meta[gsm]=dict(sample=gsm,title=title,treatment=tr,donor_pairing='NOT_DECLARED_DO_NOT_INFER_FROM_SUFFIX')
 assert len(meta)==9
 rows=[];qc=[]
 with tarfile.open(f) as t:
  members=[m for m in t.getmembers() if m.isfile()];assert len(members)==9
  for m in members:
   gsm=m.name.split('_')[0];assert gsm in meta
   with gzip.open(t.extractfile(m),'rt') as stream:df=pd.read_csv(stream,sep='\t')
   assert len(df)==227818 and df.transcript_id.is_unique and np.isfinite(df[['expected_count','TPM','FPKM']]).all().all() and (df[['expected_count','TPM','FPKM']]>=0).all().all()
   sums=df.groupby('gene_id')[['expected_count','TPM']].sum();sums['sample']=gsm;rows.append(sums.reset_index())
   qc.append(dict(sample=gsm,transcripts=len(df),genes=len(sums),sum_TPM=float(sums.TPM.sum()),sum_expected_counts=float(sums.expected_count.sum()),gzip_integrity='PASS',**{k:v for k,v in meta[gsm].items() if k!='sample'}))
 allv=pd.concat(rows,ignore_index=True);allv.to_csv(d/'ENSEMBL_SAMPLE_VALUES.tsv.gz',sep='\t',index=False)
 pd.DataFrame({'gene_id':sorted(allv.gene_id.unique())}).to_csv(d/'ENSEMBL_IDS.tsv',sep='\t',index=False);tab(qc,'EXTERNAL_SAMPLE_QC.tsv')
 run=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript','--vanilla',str(P/'code/map_external_genes.R'),str(P)],capture_output=True,text=True)
 (P/'logs/external_annotation.log').write_text(run.stdout+run.stderr)
 if run.returncode:raise RuntimeError('ANNOTATION_MAPPING_FAILED')
 mp=read(d/'ENSEMBL_SYMBOL_MAPPING.tsv').dropna().drop_duplicates()
 nmap=mp.groupby('ENSEMBL').SYMBOL.nunique();unique=mp[mp.ENSEMBL.isin(nmap[nmap==1].index)]
 nrev=unique.groupby('SYMBOL').ENSEMBL.nunique();unique=unique[unique.SYMBOL.isin(nrev[nrev==1].index)]
 tab(mp,'EXTERNAL_FULL_MAPPING.tsv');tab(unique,'EXTERNAL_ONE_TO_ONE_MAPPING.tsv')
 out=allv.merge(unique,left_on='gene_id',right_on='ENSEMBL',validate='many_to_one');out=out.merge(pd.DataFrame(meta.values()),on='sample',validate='many_to_one')
 out['group']=np.where(out.title.str.startswith('MCSF'),'MCSF',np.where(out.title.str.startswith('RANKL'),'RANKL','TYMP'))
 piv=out.pivot(index='SYMBOL',columns='sample',values='TPM')
 md=pd.DataFrame(meta.values());md['group']=np.where(md.title.str.startswith('MCSF'),'MCSF',np.where(md.title.str.startswith('RANKL'),'RANKL','TYMP'))
 groups={g:md.loc[md.group==g,'sample'].tolist() for g in ['MCSF','RANKL','TYMP']};assert all(len(v)==3 for v in groups.values())
 desc=pd.DataFrame(index=piv.index)
 for g,ids in groups.items():
  desc[g+'_mean_TPM']=piv[ids].mean(1);desc[g+'_min_TPM']=piv[ids].min(1);desc[g+'_max_TPM']=piv[ids].max(1)
 for g in ['TYMP','RANKL']:
  desc[g+'_vs_MCSF_log2_meanTPMplus1']=np.log2(desc[g+'_mean_TPM']+1)-np.log2(desc.MCSF_mean_TPM+1)
  desc[g+'_samples_above_control_mean']=(piv[groups[g]].gt(desc.MCSF_mean_TPM,axis=0)).sum(1)
 desc['TYMP_vs_RANKL_log2_meanTPMplus1']=np.log2(desc.TYMP_mean_TPM+1)-np.log2(desc.RANKL_mean_TPM+1)
 desc.reset_index().rename(columns={'SYMBOL':'gene'}).to_csv(P/'tables/EXTERNAL_ALL_GENE_DESCRIPTIVE.tsv.gz',sep='\t',index=False)
 focus=['TYMP','SPP1','LGALS3','CTSB','CTSD','APOE','ITGB1','FYN','IRF8','NFATC1','ACP5','CTSK']
 tab(desc.reindex(focus).reset_index().rename(columns={'SYMBOL':'gene'}),'EXTERNAL_FIXED_GENE_RESULTS.tsv')
 tab(piv.reindex(focus).reset_index().rename(columns={'SYMBOL':'gene'}),'EXTERNAL_FIXED_GENE_SAMPLE_TPM.tsv')
 cfg=json.loads((P/'resources/PROGRAMS.json').read_text());hp=R/'phase277_metabolic_macrophage_barrier_cooperation_20260928_v1/harmonized_platform_v1/ADDENDUM_BEFORE_HARMONIZED_MODELS.json'
 cfg['C4_platform17']=json.loads(hp.read_text())['C4_platform17']
 cov=[];programrows=[]
 for name,gs in cfg.items():
  present=[g for g in gs if g in piv.index];cov.append(dict(program=name,requested=len(gs),present=len(present),missing=';'.join(g for g in gs if g not in piv.index)))
  if len(present)<.8*len(gs):continue
  v=np.log2(piv.loc[present]+1).mean(0)
  for g,ids in groups.items():
   for sid in ids:programrows.append(dict(program=name,group=g,sample=sid,score=float(v[sid])))
 prog=pd.DataFrame(programrows);tab(prog,'EXTERNAL_PROGRAM_SAMPLE_SCORES.tsv');tab(cov,'EXTERNAL_PROGRAM_COVERAGE.tsv')
 pr=[]
 for name,g in prog.groupby('program'):
  control=g.loc[g.group=='MCSF','score']
  for treat in ['TYMP','RANKL']:
   v=g.loc[g.group==treat,'score'];pr.append(dict(program=name,treatment=treat,control_mean=control.mean(),treated_mean=v.mean(),mean_difference=v.mean()-control.mean(),treated_samples_above_control_mean=int((v>control.mean()).sum()),samples_each=3,inference='DESCRIPTIVE_ONLY_DONOR_PAIRING_UNCONFIRMED'))
 tab(pr,'EXTERNAL_PROGRAM_DESCRIPTIVE.tsv')
 # Independent aggregation of complete transcript rows for prespecified genes.
 mapping=unique.set_index('SYMBOL').ENSEMBL.to_dict();errors=[]
 with tarfile.open(f) as t:
  for m in t.getmembers():
   if not m.isfile():continue
   gsm=m.name.split('_')[0]
   with gzip.open(t.extractfile(m),'rt') as stream:df=pd.read_csv(stream,sep='\t')
   for gene in focus:
    if gene not in mapping:continue
    value=sum(df.loc[df.gene_id==mapping[gene],'TPM'].tolist());errors.append(abs(value-piv.at[gene,gsm]))
 assert max(errors)<1e-8
 dump(P/'qc/EXTERNAL_VALIDATION.json',dict(status='PASS',archive_sha256=sha(f),metadata_sha256=sha(d/'GSE171542_samples.soft'),samples=9,transcripts_each=227818,genes_one_to_one=len(unique),independent_focus_aggregates=len(errors),max_error=max(errors),inference='NO_P_VALUES_SMALL_N_AND_PAIRING_UNCONFIRMED',annotation_source='org.Hs.eg.db; full mapping and R session retained'))
 print('EXTERNAL_COMPLETE',flush=True)
if __name__=='__main__':main()
