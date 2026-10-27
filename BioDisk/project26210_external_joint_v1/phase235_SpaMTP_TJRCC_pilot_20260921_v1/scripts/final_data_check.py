from pathlib import Path
import json,hashlib,datetime
import numpy as np,pandas as pd
P=Path(__file__).resolve().parents[1];R=P.parent
def save(n,o):(P/n).write_text(json.dumps(o,ensure_ascii=False,indent=2))
def tab(n,d):d.to_csv(P/'tables'/n,sep='\t',index=False)
def main():
 q=pd.read_csv(P/'tables/SPECTRAL_TECHNICAL_QC.tsv',sep='\t')
 cand=pd.read_csv(P/'tables/TARGETED_MASS_CANDIDATES_ALL.tsv',sep='\t')
 anno=pd.read_excel(P/'downloads/spatial_metabolites_annotation.xlsx')
 summaries=[];regionrows=[];anchecks=[]
 for row in q.itertuples():
  label=row.sample.replace('_','-')+'-'+row.polarity;z=np.load(P/'results/matrices'/(label+'.npz'));X=z['intensity'];mz=z['mz'];xy=z['xy'];tic=z['tic'];mask=tic>0
  if row.polarity=='neg':
   np.testing.assert_allclose(mz,anno.mz.to_numpy(),atol=5.01e-6,rtol=0)
   anchecks.append(dict(sample=row.sample,polarity='neg',author_rows=len(anno),axis_rows=len(mz),max_mass_rounding_error=float(np.max(abs(mz-anno.mz.to_numpy()))),matched=True))
  selected=cand[(cand['sample']==row.sample)&(cand.polarity==row.polarity)]
  for molecule,g in selected.groupby('molecule'):
   hits=g[g.within_5ppm]
   summaries.append(dict(sample=row.sample,polarity=row.polarity,molecule=molecule,eligible_mass_peak_count=len(hits),status='mass-compatible candidate;not confirmed identity' if len(hits) else 'no candidate in deposited feature list under frozen adducts/5ppm',minimum_abs_ppm=float(g.ppm.abs().min()),best_observed_mz=float(g.iloc[g.ppm.abs().argmin()].observed_mz)))
  matches=selected[selected.within_5ppm]
  reg=pd.read_csv(P/'tables'/(label+'_SPAMTP_GRAPH_REGIONS.tsv'),sep='\t').set_index(['x','y'])
  for c in matches.itertuples():
   j=c.feature_index;intensity=X[:,j];norm=np.divide(intensity,tic,out=np.zeros_like(tic),where=mask)*len(mz)
   bins=np.floor((xy-1)/4).astype(int)
   counts={}
   for region in [1,2,3,4]:
    chosen=mask & np.array([int(reg.loc[(a,b),'region'])==region if (a,b) in reg.index else False for a,b in bins])
    vals=norm[chosen]
    regionrows.append(dict(sample=row.sample,polarity=row.polarity,molecule=c.molecule,mz=c.observed_mz,ppm=c.ppm,technical_region=region,pixels=int(chosen.sum()),positive_pixels=int((intensity[chosen]>0).sum()),positive_fraction=float((intensity[chosen]>0).mean()) if chosen.any() else np.nan,mean_TIC_normalized=float(vals.mean()) if len(vals) else np.nan,median_TIC_normalized=float(np.median(vals)) if len(vals) else np.nan,unit='descriptive pixels, not biological replicates',circularity='regions derived from SM features;no differential P or independent enrichment claim'))
 tab('TARGET_MOLECULE_COVERAGE.tsv',pd.DataFrame(summaries));tab('CANDIDATE_REGION_DESCRIPTIVE.tsv',pd.DataFrame(regionrows));tab('AUTHOR_NEGATIVE_ANNOTATION_AXIS_CHECK.tsv',pd.DataFrame(anchecks))
 # Read-only preservation of pre-existing RDS and spatial inputs.
 old=json.loads((R/'phase234_part4_replan_execution_20260921_v1/audit/RDS_BEFORE.json').read_text())
 for x in old:
  s=Path(x['path']).stat();assert(s.st_size,s.st_mtime_ns,s.st_ino)==(x['bytes'],x['mtime_ns'],x['inode'])
 for x in json.loads((P/'audit/READONLY_ST_SOURCES.json').read_text()):
  f=Path(x['path']);s=f.stat();assert s.st_size==x['bytes'] and s.st_mtime_ns==x['mtime_ns'];assert hashlib.sha256(f.read_bytes()).hexdigest()==x['sha256']
 # Check downloaded raw file identity and source code lock once at completion.
 for x in json.loads((P/'audit/EXTRACTED_FILE_INTEGRITY.json').read_text()):
  assert hashlib.sha256((P/'extracted'/x['file']).read_bytes()).hexdigest()==x['sha256']
 for x in json.loads((P/'audit/OFFICIAL_SOURCE_HASHES.json').read_text()):
  assert hashlib.sha256((P/'metadata'/x['file']).read_bytes()).hexdigest()==x['sha256']
 save('PILOT_DATA_VALIDATION.json',dict(status='PASS',time=datetime.datetime.now().astimezone().isoformat(),raw_files=9,imzml_ibd_pairs=4,pixels=int(q.pixels.sum()),nonzero_pixels=int(q.nonzero_pixels.sum()),official_source_files=6,author_annotation_axis_pairs=2,old_RDS_unchanged=len(old),original_ST_unchanged=True,biological_inference=False,new_figures=0))
 print('PILOT_DATA_VALIDATION PASS',flush=True)
if __name__=='__main__':main()
