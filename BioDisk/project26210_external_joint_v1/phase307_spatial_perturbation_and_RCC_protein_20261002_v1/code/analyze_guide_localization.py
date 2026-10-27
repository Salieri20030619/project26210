from pathlib import Path
import json,datetime,hashlib
import numpy as np,pandas as pd
from scipy.spatial import cKDTree
from sklearn.neighbors import BallTree
P=Path(__file__).resolve().parents[1]
proto={
'created':datetime.datetime.now().astimezone().isoformat(),
'relationship_to_v1':'V1 phenotype eligibility failed in all6sections, before any eligible effect. This is a distinct spatial-signal localization analysis; V1 gates/results retained. Assay known to deliver guides in OT1, but deposited mixed bins cannot establish identity of every signal.',
'question':'Do CD44-targeting sgRNA transcript locations differ from non-targeting transcript locations relative to Spp1-detected macrophage-associated bins?',
'unit':'sgRNA-bearing spatial bin, not a verified KO T cell. Section/timepoint effects only; donor independence unresolved. No cell/bin/section P/FDR.',
'groups':'Use all RNA-QC spatial coordinates, no T expression gate for localization. Main exactly1 positive target gene among34 and targetsgRNA>0. Preserved sensitivity requires targetguide>=2UMIs and targetshare>=80%. All overlap/ambiguous sites excluded by these fixed rules.',
'source':'Reuse V1 qualified macro and Spp1-detected source masks without modification. Main>=20locations in eachguide group and >=20Spp1-detected source bins. All six sections retained including late zero capture.',
'primary':'Difference of median log1p(nearest Spp1-detected macrophage distance / gridstep): Cd44 minus NTC. Gridstep=median positive unique-x-coordinate spacing computed from stored RNA geometry, recorded. Positive means farther.',
'specificity':'Subtract analogous all-macrophage distance effect from primary. Also retain raw medians, Spp1-undetected source and fibroblast distances, same-bin Spp1-macro colocalization fractions. No choosing best reference.',
'limitations':'Most guide locations have1UMI and lack purified recipient transcriptomes. Ambient/cross-bin/mixed-cell signal possible; no new functional scores, no causal receptor-response signature. Ratio if denominatorzero not manufactured. No lowered V1 threshold.'}
f=P/'resources/GUIDE_LOCALIZATION_PROTOCOL.json';assert not f.exists();f.write_text(json.dumps(proto,indent=2))
effects=[];counts=[];checks=[];source_effects=[]
def ck(n,ok):
 checks.append({'check':n,'passed':bool(ok)})
 if not ok:raise AssertionError(n)
for file in sorted((P/'inputs').glob('*_audit.npz')):
 sample=file.stem.replace('_audit','');z=np.load(file);xy=z['coordinates'];g=z['guide_counts'];names=z['guide_names'].tolist();total=g.sum(1);ng=(g>0).sum(1);step=float(np.median(np.diff(np.unique(xy[:,0]))));ck(sample+'_positive_gridstep',step>0);ci=names.index('sgCd44');ni=names.index('sgnon-targeting');m=z['macro'];mp=z['Spp1_detected_macro'];groups={};dist={}
 for endpoint,mask in [('Spp1_detected_M',mp),('all_M',m),('Spp1_undetected_M',m&~mp),('Fibroblast',(z['annotation']=='Fibroblast')&(z['RNA_total']>=20))]:
  dd=z['distance_'+endpoint];dist[endpoint]=dd
  selected=((g[:,ci]>0)|(g[:,ni]>0));new=BallTree(xy[mask]).query(xy[selected],k=1)[0].ravel();ck(sample+endpoint+'_independent_geometry',np.allclose(new,dd[selected],rtol=1e-8,atol=1e-9))
 for rule in ['single_target','dominant80_min2']:
  groups={}
  for name,j in [('sgCd44',ci),('NTC',ni)]:groups[name]=((g[:,j]>0)&(ng==1)) if rule=='single_target' else ((g[:,j]>=2)&(g[:,j]>=.8*total))
  n={k:int(v.sum()) for k,v in groups.items()};ok=min(n.values())>=20 and mp.sum()>=20
  counts.append({'sample':sample,'timepoint':sample.split('_')[0],'assignment':rule,'Cd44_sites':n['sgCd44'],'NTC_sites':n['NTC'],'Cd44_guideUMIs':int(g[groups['sgCd44'],ci].sum()),'NTC_guideUMIs':int(g[groups['NTC'],ni].sum()),'Cd44_singleUMI_fraction':float((g[groups['sgCd44'],ci]==1).mean()) if n['sgCd44'] else np.nan,'Spp1_detected_macro_bins':int(mp.sum()),'gridstep_native':step,'eligible':bool(ok),'reason':'' if ok else 'INSUFFICIENT_LOCATIONS_OR_SOURCE'})
  ee={}
  for endpoint,dd in dist.items():
   vals={k:float(np.median(np.log1p(dd[v]/step))) if n[k] else np.nan for k,v in groups.items()};raw={k:float(np.median(dd[v])) if n[k] else np.nan for k,v in groups.items()};delta=vals['sgCd44']-vals['NTC'];ee[endpoint]=delta
   effects.append({'sample':sample,'timepoint':sample.split('_')[0],'assignment':rule,'endpoint':'distance_'+endpoint,'Cd44_value':vals['sgCd44'],'NTC_value':vals['NTC'],'Cd44_raw_median_distance':raw['sgCd44'],'NTC_raw_median_distance':raw['NTC'],'raw_distance_ratio':raw['sgCd44']/raw['NTC'] if raw['NTC']>0 else np.nan,'effect':delta if ok else np.nan,'scale':'median_log1p_gridnormalized_distance_difference','eligible':bool(ok),'P':np.nan,'FDR':np.nan})
   for k,v in groups.items():
    # Independently reconstruct the transformed median from sorted source distances.
    values=sorted(dd[v].tolist());mid=len(values)//2
    calc=(np.log1p(values[mid]/step) if len(values)%2 else (np.log1p(values[mid-1]/step)+np.log1p(values[mid]/step))/2) if len(values) else np.nan;ck(sample+rule+endpoint+k+'_median',np.allclose(calc,vals[k],equal_nan=True))
  effects.append({'sample':sample,'timepoint':sample.split('_')[0],'assignment':rule,'endpoint':'Spp1_specific_excess_distance','Cd44_value':np.nan,'NTC_value':np.nan,'effect':ee['Spp1_detected_M']-ee['all_M'] if ok else np.nan,'scale':'difference_of_distance_effects','eligible':bool(ok),'P':np.nan,'FDR':np.nan})
  vv={k:float(mp[v].mean()) if n[k] else np.nan for k,v in groups.items()};effects.append({'sample':sample,'timepoint':sample.split('_')[0],'assignment':rule,'endpoint':'Spp1_macro_same_bin_fraction','Cd44_value':vv['sgCd44'],'NTC_value':vv['NTC'],'effect':vv['sgCd44']-vv['NTC'] if ok else np.nan,'scale':'fraction_difference','eligible':bool(ok),'P':np.nan,'FDR':np.nan})
 print('LOCALIZATION_DONE',sample,flush=True)
r=pd.DataFrame(effects);c=pd.DataFrame(counts);r.to_csv(P/'tables/GUIDE_LOCALIZATION_EFFECTS.tsv',sep='\t',index=False);c.to_csv(P/'tables/GUIDE_LOCALIZATION_QUALIFICATION.tsv',sep='\t',index=False);summary=[]
for (rule,ep),df in r.groupby(['assignment','endpoint']):
 q=df[df.eligible&df.effect.notna()];summary.append({'assignment':rule,'endpoint':ep,'qualified_sections':len(q),'positive_sections':int((q.effect>0).sum()),'negative_sections':int((q.effect<0).sum()),'median_effect':q.effect.median(),'P':np.nan,'FDR':np.nan,'independent_mouse_n':'unresolved'})
s=pd.DataFrame(summary);s.to_csv(P/'tables/GUIDE_LOCALIZATION_SUMMARY.tsv',sep='\t',index=False);ck('all72_endpoints',len(r)==72);ck('no_inference',r.P.isna().all() and r.FDR.isna().all());ck('no_relabeling',not (P/'inputs/relabeled.h5ad').exists());(P/'qc/LOCALIZATION_VALIDATION_CHECKS.json').write_text(json.dumps(checks,indent=2));(P/'qc/LOCALIZATION_COMPLETE.json').write_text(json.dumps({'status':'PASS','checks':len(checks),'time':datetime.datetime.now().astimezone().isoformat(),'protocol_sha256':hashlib.sha256(f.read_bytes()).hexdigest()},indent=2));print('LOCALIZATION_RESULTS',r[r.assignment.eq('single_target')&r.endpoint.isin(['distance_Spp1_detected_M','Spp1_specific_excess_distance'])].to_string(index=False),flush=True)
