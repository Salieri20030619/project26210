from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='4'
import json,datetime,subprocess,hashlib
import numpy as np,pandas as pd,anndata as ad
from sklearn.neighbors import BallTree
P=Path(__file__).resolve().parents[1];checks=[]
def ck(name,ok):
 checks.append({'check':name,'passed':bool(ok)})
 if not ok:raise AssertionError(name)
def eq(a,b):return np.allclose(a,b,rtol=1e-8,atol=1e-9,equal_nan=True)
ck('analysis_exit',(P/'logs/spatial_analysis.exit').read_text().strip()=='0')
r=pd.read_csv(P/'tables/SPATIAL_PERTURBATION_EFFECTS.tsv',sep='\t');qual=pd.read_csv(P/'tables/SECTION_QUALIFICATION.tsv',sep='\t');cov=pd.read_csv(P/'tables/RECIPIENT_GENE_COUNTS.tsv',sep='\t');proto=json.loads((P/'resources/SPATIAL_PERTURBATION_PROTOCOL.json').read_text())
for sample in qual['sample'].unique():
 subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],check=True,stdout=subprocess.DEVNULL)
 z=np.load(P/f'inputs/{sample}_audit.npz');xy=z['coordinates'];genes=z['genes'].tolist();gg=z['guide_names'].tolist();a=ad.read_h5ad(P/f'downloads/{sample}.h5');g=ad.read_h5ad(P/f'downloads/{sample}.guide.h5',backed='r');x=a.X.tocsr();tot=np.asarray(x.sum(1)).ravel();gene_mat=x[:,a.var_names.get_indexer(genes)].toarray();ck(sample+'_original_counts',eq(gene_mat,z['gene_counts']));ck(sample+'_original_RNA_total',eq(tot,z['RNA_total']));ck(sample+'_RNA_coordinates',np.array_equal(xy,a.obsm['spatial']));ck(sample+'_RNA_barcodes',np.array_equal(a.obs_names.to_numpy(dtype=str),z['barcodes']))
 coord_map={tuple(v):i for i,v in enumerate(g.obsm['spatial'])};indices=np.array([coord_map[tuple(v)] for v in xy]);ck(sample+'_coordinate_join',np.array_equal(indices,z['guide_row_index']));order=np.argsort(indices);inv=np.argsort(order);guide=np.asarray(g.X[indices[order],:])[inv];ck(sample+'_guide_counts',np.array_equal(guide,z['guide_counts']));g.file.close()
 ann=json.loads((P/f'downloads/{sample}.annotation.json').read_text());labels=np.array([ann[b] for b in a.obs_names]);ck(sample+'_annotation_join',np.array_equal(labels,z['annotation']));support=lambda gs:(gene_mat[:,[genes.index(s) for s in gs]]>0).sum(1)>0;m=(labels=='Macrophages')&support(['C1qa','C1qb','C1qc','Csf1r','Adgre1'])&(tot>=20);mp=m&(gene_mat[:,genes.index('Spp1')]>0);t=(labels=='T')&support(['Cd3d','Cd3e','Trac'])&(tot>=20);ck(sample+'_identity',np.array_equal(m,z['macro']) and np.array_equal(mp,z['Spp1_detected_macro']) and np.array_equal(t,z['T']))
 for rule in ['single_target','dominant80_min2']:
  groups={}
  for name,col in [('sgCd44','sgCd44'),('NTC','sgnon-targeting')]:
   j=gg.index(col);mask=t&((guide[:,j]>0)&((guide>0).sum(1)==1) if rule=='single_target' else (guide[:,j]>=2)&(guide[:,j]>=.8*guide.sum(1)));groups[name]=mask;ck(sample+rule+name+'_assignment',np.array_equal(mask,z['group_'+rule]==name))
  qu=qual[qual['sample'].eq(sample)&qual.assignment.eq(rule)].iloc[0];ok=min(v.sum() for v in groups.values())>=20 and mp.sum()>=20;ck(sample+rule+'_qualification',qu.eligible==ok and qu.Cd44_bins==groups['sgCd44'].sum() and qu.NTC_bins==groups['NTC'].sum());rows=r[r['sample'].eq(sample)&r.assignment.eq(rule)].set_index('endpoint');ratios={}
  for endpoint,source_mask in [('Spp1_detected_M',mp),('all_M',m),('Spp1_undetected_M',m&~mp),('Fibroblast',(labels=='Fibroblast')&(tot>=20))]:
   source_xy=xy[source_mask];tree=BallTree(source_xy) if len(source_xy) else None;med={}
   for name,mask in groups.items():
    dist=tree.query(xy[mask],k=1)[0].ravel() if tree is not None and mask.any() else np.full(mask.sum(),np.nan);ck(sample+rule+name+endpoint+'_independent_geometry',eq(dist,z['distance_'+endpoint][mask]));med[name]=np.median(dist) if len(dist) else np.nan
   rho=np.log2(med['sgCd44']/med['NTC']) if min(med.values())>0 else np.nan;ratios[endpoint]=rho;row=rows.loc['distance_'+endpoint];ck(sample+rule+endpoint+'_effect',eq(rho if ok else np.nan,row.effect));ck(sample+rule+endpoint+'_groupvalues',eq([med['sgCd44'],med['NTC']],[row.Cd44_value,row.NTC_value]))
  ck(sample+rule+'_specificity',eq(ratios['Spp1_detected_M']-ratios['all_M'] if ok else np.nan,rows.loc['Spp1_specific_excess_distance','effect']))
  for endpoint,gs in (proto['recipient_programs']|{'Cd44_target_RNA':['Cd44']}).items():
   values={}
   for name,mask in groups.items():
    sums=np.asarray(x[mask][:,a.var_names.get_indexer(gs)].sum(0)).ravel();values[name]=float(np.log1p(sums*1e4/max(tot[mask].sum(),1)).mean()) if mask.any() else np.nan
   row=rows.loc[endpoint];ck(sample+rule+endpoint+'_raw_pseudobulk',eq([values['sgCd44'],values['NTC']],[row.Cd44_value,row.NTC_value]));ck(sample+rule+endpoint+'_effect',eq(values['sgCd44']-values['NTC'] if ok else np.nan,row.effect))
  for row in cov[cov['sample'].eq(sample)&cov.assignment.eq(rule)].itertuples():
   v=gene_mat[groups[row.group],genes.index(row.gene)];ck(sample+rule+row.group+row.gene+'_coverage',v.sum()==row.total_gene_counts and (v>0).sum()==row.gene_detected_bins and eq(tot[groups[row.group]].sum(),row.total_RNA_counts))
 print('VALIDATED_SECTION',sample,flush=True)
s=pd.read_csv(P/'tables/DESCRIPTIVE_SECTION_SUMMARY.tsv',sep='\t')
for i,row in s.iterrows():
 v=r[r.assignment.eq(row.assignment)&r.endpoint.eq(row.endpoint)&r.eligible].effect.dropna();ck(f'summary{i}',len(v)==row.qualified_sections and (v>0).sum()==row.positive_sections and (v<0).sum()==row.negative_sections and eq(v.median(),row.median_effect))
ck('96_complete_effect_rows',len(r)==96);ck('no_bin_or_biological_P',r.P.isna().all() and r.FDR.isna().all() and s.P.isna().all());ck('24_officialMD5',len(json.loads((P/'qc/OFFICIAL_MD5_VALIDATION.json').read_text())['files'])==24)
(P/'qc/VALIDATION_CHECKS.json').write_text(json.dumps(checks,indent=2));(P/'qc/FINAL_VALIDATION.json').write_text(json.dumps({'status':'PASS','checks':len(checks),'time':datetime.datetime.now().astimezone().isoformat(),'method':'Original H5 counts and independent tuple coordinate join; BallTree versus cKDTree geometry; raw pseudobulk; complete descriptive effects and no inferential P'},indent=2));print('VALIDATION_PASS',len(checks),flush=True)
