from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase186_MTA_fibroblast_panel_20260916_v1')
df=pd.read_csv(r/'sources/counts_complete.txt.gz',sep='\t',comment='#').set_index('Geneid')
cols=['C1_1','C2_1','C3_1','C4_2','M1_1','M2_1','M3_1','M4_1'];x=df[cols].to_numpy(dtype=float)
assert df.index.is_unique and np.isfinite(x).all() and (x>=0).all() and (x==np.floor(x)).all()
positive=(x>0).all(axis=1)
gm=np.exp(np.mean(np.log(x[positive]),axis=1))
sf=np.median(x[positive]/gm[:,None],axis=0);sf/=np.exp(np.mean(np.log(sf)))
uq=np.array([np.quantile(x[x[:,i]>0,i],.75) for i in range(8)]);uq/=np.exp(np.mean(np.log(uq)))
maps=json.loads((r/'GENE_MAPPING_HGNC.json').read_text());assert len(maps)==11
values=[];results=[];qc=[]
for j,c in enumerate(cols):
 qc.append({'library':c,'group':'Control' if j<4 else 'MTA','sum_counts':x[:,j].sum(),'detected_genes':(x[:,j]>0).sum(),'median_ratio_factor':sf[j],'upper_quartile_factor':uq[j],'donor_id':'unknown','paired_id':'not_assigned'})
for a in maps:
 g=a['symbol'];eid=a['ensembl_id'];raw=df.loc[eid,cols].to_numpy(dtype=float)
 for j,c in enumerate(cols):values.append({'gene':g,'ensembl_id':eid,'library':c,'group':'Control' if j<4 else 'MTA','raw_count':raw[j],'median_ratio_normalized':raw[j]/sf[j],'upper_quartile_normalized':raw[j]/uq[j]})
 for meth,factors in [('median_ratio',sf),('upper_quartile',uq)]:
  n=raw/factors;l=np.log2(n+1);c=l[:4];t=l[4:]
  results.append({'gene':g,'method':meth,'control_median_normalized':np.median(n[:4]),'MTA_median_normalized':np.median(n[4:]),'delta_log2_median':np.median(t)-np.median(c),'control_log2_min':min(c),'control_log2_max':max(c),'MTA_log2_min':min(t),'MTA_log2_max':max(t),'complete_separation':'all_MTA_above' if min(t)>max(c) else 'all_MTA_below' if max(t)<min(c) else 'overlap','minimum_raw':min(raw)})
pd.DataFrame(values).to_csv(r/'PANEL_SAMPLE_VALUES.tsv',sep='\t',index=False)
pd.DataFrame(results).to_csv(r/'PANEL_EFFECTS.tsv',sep='\t',index=False)
pd.DataFrame(qc).to_csv(r/'LIBRARY_QC.tsv',sep='\t',index=False)
corr=pd.DataFrame(np.corrcoef(np.log2(x[positive]/sf+1).T),index=cols,columns=cols)
corr.to_csv(r/'LIBRARY_CORRELATION.tsv',sep='\t')
val={'genes':len(df),'libraries':8,'all_positive_reference_genes':int(positive.sum()),'no_negative_missing_fractional_counts':True,'unique_geneids':True,'panel_coverage':11,'no_inferential_tests':True,'min_library_correlation':float(corr.to_numpy()[np.triu_indices(8,1)].min()),'input_sha256':hashlib.sha256((r/'sources/counts_complete.txt.gz').read_bytes()).hexdigest()}
(r/'ANALYSIS_VALIDATION.json').write_text(json.dumps(val,indent=2))
print(json.dumps(val));print(pd.DataFrame(results).query("method=='median_ratio'")[['gene','control_median_normalized','MTA_median_normalized','delta_log2_median','complete_separation','minimum_raw']].to_string(index=False))
