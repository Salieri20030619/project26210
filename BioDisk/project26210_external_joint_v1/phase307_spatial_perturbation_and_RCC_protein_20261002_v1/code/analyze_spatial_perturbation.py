from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='4'
import json,datetime,hashlib,re,subprocess,gc,traceback
import numpy as np,pandas as pd,anndata as ad,openpyxl
from scipy.spatial import cKDTree
P=Path(__file__).resolve().parents[1]
for d in ['inputs','tables','qc','resources','logs']:(P/d).mkdir(exist_ok=True)
proto={
'created':datetime.datetime.now().astimezone().isoformat(),
'scope':'Independent descriptive reproduction of external mouse OT1 spatial CRISPR perturbation. All six public Day4/7/10 rep1/2 sections; no renal or TYMP causal claim.',
'biological_unit':'A section is an acquisition, not automatically an independent mouse. Underlying three BioSample IDs and two sections/day have unresolved donor mapping. Report each section and timepoint; no pooled cell/bin P, no new biological P/FDR, no six-mouse claim.',
'spatial_unit':'Author fixed bins annotated as cell type; official checksum workbook calls annotation10um/mapping10um-to50um while CSV names20um/100um. Native stored coordinate units are retained, distance ratios dimensionless. Do not call bins segmented single cells or infer micron scale.',
'alignment':'RNA QC row order and barcode suffix differ from guide matrix. Exact unique (x,y) within each named section, every RNA bin must map to guide coordinate. Annotation by original RNA barcode, never suffix-stripping or row-order join.',
'source':'Author Macrophages annotation with >=1 raw transcript among C1qa,C1qb,C1qc,Csf1r,Adgre1 and RNA total>=20. Spp1-detected if rawSpp1>0, undetected otherwise; detection is RNA not protein secretion.',
'recipient':'Author T annotation, >=1 raw transcript among Cd3d,Cd3e,Trac and RNA total>=20. sgCd44 versus sgnon-targeting. Main assignment exactly one gene-target column>0 across34; ambiguous perturbation bins excluded. Sensitivity dominant target >=2guideUMIs and >=80%of all guidecounts. Both rules fixed before effect calculation.',
'eligibility':'Each rule requires >=20 T bins/guide group and >=20 Spp1-detected macrophage bins for section endpoint. Retain unqualified and all opposite effects.',
'primary_endpoint':'log2 ratio of median nearest Spp1-detected macrophage distance (sgCd44 / non-targeting), expectedpositive.',
'density_control':'Also distance to all macrophages and to Spp1-undetected macrophages. Spp1-specific excess=log2(Cd44/NTC median distance to detected M) minus log2(Cd44/NTC median distance to all M); expectedpositive. AllM and undetected endpoints retained.',
'recipient_programs':{'effector':['Ifng','Gzmb','Prf1'],'dysfunction':['Pdcd1','Havcr2','Lag3','Tigit','Tox','Entpd1']},
'program_method':'Per section and perturbation aggregate rawRNA in qualified T bins, mean log1pCP10K over all fixed genes. Effect KO minuscontrol; no significance-based gene selection. Show gene coverage/detection and zero scores. Local mixed-bin RNA is not direct killing function.',
'additional_descriptive':'Cd44 RNA expression as target readout; fibroblast distance reported for context, not barrier structure; Spp1 detected M count and composition retained. No tuned radius or clustering.',
'statistics':'No inferential P or FDR at section/bin level. Median effects and section direction counts, with biological independence unresolved. Previous RCC negative/reverse results remain.'}
f=P/'resources/SPATIAL_PERTURBATION_PROTOCOL.json'
if f.exists():assert {k:v for k,v in json.loads(f.read_text()).items() if k!='created'}=={k:v for k,v in proto.items() if k!='created'}
else:f.write_text(json.dumps(proto,indent=2))
done=json.loads((P/'qc/OT1_DOWNLOAD_COMPLETE.json').read_text());assert done['files']==24;assert (P/'logs/download_ot1.exit').read_text().strip()=='0'
w=openpyxl.load_workbook(P/'metadata/OT1_md5.xlsx',read_only=True,data_only=True);official={}
for row in w.active.values:
 for val in row:
  if isinstance(val,str):
   match=re.fullmatch(r'MD5 \((.+)\) = ([a-f0-9]{32})',val.strip())
   if match:official[match[1]]=match[2]
ver=[]
for rr in done['records']:
 h=hashlib.md5()
 with (P/'downloads'/rr['file']).open('rb') as inf:
  for block in iter(lambda:inf.read(8*1024*1024),b''):h.update(block)
 assert h.hexdigest()==official[rr['file']],rr['file'];ver.append({'file':rr['file'],'official_MD5':official[rr['file']],'verified':True})
(P/'qc/OFFICIAL_MD5_VALIDATION.json').write_text(json.dumps({'status':'PASS','files':ver},indent=2));print('OFFICIAL_MD5_PASS',len(ver),flush=True)
genes=sorted(set(['Spp1','Cd44','Cd3d','Cd3e','Trac','C1qa','C1qb','C1qc','Csf1r','Adgre1']+sum(proto['recipient_programs'].values(),[])))
all_results=[];eligibility=[];coverage=[];identities=[]
for sample in ['Day4_rep1','Day4_rep2','Day7_rep1','Day7_rep2','Day10_rep1','Day10_rep2']:
 subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],check=True,stdout=subprocess.DEVNULL)
 a=ad.read_h5ad(P/f'downloads/{sample}.h5');g=ad.read_h5ad(P/f'downloads/{sample}.guide.h5');xy=np.asarray(a.obsm['spatial']);gi=pd.MultiIndex.from_arrays(g.obsm['spatial'].T);ai=pd.MultiIndex.from_arrays(xy.T);assert gi.is_unique and ai.is_unique;idx=gi.get_indexer(ai);assert (idx>=0).all();guide=np.asarray(g.X[idx,:]);assert np.isfinite(guide).all() and (guide>=0).all() and (guide==np.floor(guide)).all();names=g.var_names.tolist();assert names==sorted(names)
 annotation=json.loads((P/f'downloads/{sample}.annotation.json').read_text());labels=np.array([annotation.get(b,'UNANNOTATED') for b in a.obs_names]);assert (labels!='UNANNOTATED').all();x=a.X.tocsr();assert np.isfinite(x.data).all() and (x.data>=0).all() and (x.data==np.floor(x.data)).all();tot=np.asarray(x.sum(1)).ravel();jj=a.var_names.get_indexer(genes);assert (jj>=0).all();z=x[:,jj].toarray();gene=lambda name:z[:,genes.index(name)];support=lambda gs:np.any(z[:,[genes.index(s) for s in gs]]>0,axis=1)
 macro=(labels=='Macrophages')&support(['C1qa','C1qb','C1qc','Csf1r','Adgre1'])&(tot>=20);mp=macro&(gene('Spp1')>0);mn=macro&~mp;t=(labels=='T')&support(['Cd3d','Cd3e','Trac'])&(tot>=20);fib=(labels=='Fibroblast')&(tot>=20)
 distances={name:(cKDTree(xy[mask]).query(xy,workers=4)[0] if mask.any() else np.full(len(xy),np.nan)) for name,mask in [('Spp1_detected_M',mp),('all_M',macro),('Spp1_undetected_M',mn),('Fibroblast',fib)]}
 guide_n=(guide>0).sum(1);guide_total=guide.sum(1);ci=names.index('sgCd44');ni=names.index('sgnon-targeting');assignment={}
 for rule in ['single_target','dominant80_min2']:
  code=np.full(len(xy),'unassigned',dtype='U10')
  for col,label in [(ci,'sgCd44'),(ni,'NTC')]:
   mask=t&((guide[:,col]>0)&(guide_n==1) if rule=='single_target' else (guide[:,col]>=2)&(guide[:,col]/np.maximum(guide_total,1)>=.8));code[mask]=label
  assignment[rule]=code
  cnt={label:int((code==label).sum()) for label in ['sgCd44','NTC']};ok=min(cnt.values())>=20 and mp.sum()>=20;reason='' if ok else 'INSUFFICIENT_TARGET_OR_SOURCE_BINS';eligibility.append({'sample':sample,'timepoint':sample.split('_')[0],'assignment':rule,'RNA_bins':len(xy),'guide_bins':g.n_obs,'guide_aligned':len(idx),'T_qualified_bins':int(t.sum()),'macro_qualified_bins':int(macro.sum()),'Spp1_detected_macro_bins':int(mp.sum()),'Spp1_undetected_macro_bins':int(mn.sum()),'Cd44_bins':cnt['sgCd44'],'NTC_bins':cnt['NTC'],'eligible':bool(ok),'reason':reason})
  for group in ['sgCd44','NTC']:
   mask=code==group
   for j,gene_name in enumerate(genes):coverage.append({'sample':sample,'assignment':rule,'group':group,'gene':gene_name,'bins':int(mask.sum()),'total_gene_counts':int(z[mask,j].sum()),'gene_detected_bins':int((z[mask,j]>0).sum()),'total_RNA_counts':float(tot[mask].sum())})
  values={}
  for end,dd in distances.items():
   med={group:float(np.median(dd[code==group])) if cnt[group] else np.nan for group in cnt};val=float(np.log2(med['sgCd44']/med['NTC'])) if min(med.values())>0 else np.nan;values[end]=val
   all_results.append({'sample':sample,'timepoint':sample.split('_')[0],'assignment':rule,'endpoint':'distance_'+end,'Cd44_value':med['sgCd44'],'NTC_value':med['NTC'],'effect':val if ok else np.nan,'effect_scale':'log2_median_distance_ratio','eligible':bool(ok),'reason':reason,'P':np.nan,'FDR':np.nan})
  all_results.append({'sample':sample,'timepoint':sample.split('_')[0],'assignment':rule,'endpoint':'Spp1_specific_excess_distance','Cd44_value':np.nan,'NTC_value':np.nan,'effect':values['Spp1_detected_M']-values['all_M'] if ok else np.nan,'effect_scale':'difference_of_log2_ratios','eligible':bool(ok),'reason':reason,'P':np.nan,'FDR':np.nan})
  for end,gs in (proto['recipient_programs']|{'Cd44_target_RNA':['Cd44']}).items():
   val={group:float(np.mean(np.log1p(z[code==group][:,[genes.index(s) for s in gs]].sum(0)*1e4/max(tot[code==group].sum(),1)))) if cnt[group] else np.nan for group in cnt}
   all_results.append({'sample':sample,'timepoint':sample.split('_')[0],'assignment':rule,'endpoint':end,'Cd44_value':val['sgCd44'],'NTC_value':val['NTC'],'effect':val['sgCd44']-val['NTC'] if ok else np.nan,'effect_scale':'mean_logCP10K_difference','eligible':bool(ok),'reason':reason,'P':np.nan,'FDR':np.nan})
 np.savez_compressed(P/f'inputs/{sample}_audit.npz',coordinates=xy,barcodes=a.obs_names.to_numpy(dtype=str),annotation=labels.astype(str),RNA_total=tot,genes=np.array(genes),gene_counts=z,guide_names=np.array(names),guide_counts=guide,guide_row_index=idx,macro=macro,Spp1_detected_macro=mp,T=t,**{'group_'+k:v for k,v in assignment.items()},**{'distance_'+k:v for k,v in distances.items()})
 identities.append({'sample':sample,'RNA_marker':a.obs.marker.value_counts().to_dict(),'guide_marker':g.obs.marker.value_counts().to_dict(),'annotation_classes':pd.Series(labels).value_counts().to_dict(),'exact_coordinate_matches':len(idx),'RNA_bins':a.n_obs,'guide_bins':g.n_obs,'biological_mouse_n':'unresolved','unit':'author_annotated_spatial_bins'})
 print('SECTION_COMPLETE',sample,[(r['assignment'],r['Cd44_bins'],r['NTC_bins']) for r in eligibility if r['sample']==sample],flush=True);del a,g,x,z,guide,distances;gc.collect()
pd.DataFrame(all_results).to_csv(P/'tables/SPATIAL_PERTURBATION_EFFECTS.tsv',sep='\t',index=False);pd.DataFrame(eligibility).to_csv(P/'tables/SECTION_QUALIFICATION.tsv',sep='\t',index=False);pd.DataFrame(coverage).to_csv(P/'tables/RECIPIENT_GENE_COUNTS.tsv',sep='\t',index=False);(P/'qc/INPUT_IDENTITY.json').write_text(json.dumps(identities,indent=2))
r=pd.DataFrame(all_results);summ=[]
for (rule,ep),df in r.groupby(['assignment','endpoint']):
 q=df[df.eligible&df.effect.notna()];summ.append({'assignment':rule,'endpoint':ep,'qualified_sections':len(q),'positive_sections':int((q.effect>0).sum()),'negative_sections':int((q.effect<0).sum()),'median_effect':float(q.effect.median()),'range_min':float(q.effect.min()),'range_max':float(q.effect.max()),'independent_mouse_n':'unresolved','P':np.nan,'FDR':np.nan})
pd.DataFrame(summ).to_csv(P/'tables/DESCRIPTIVE_SECTION_SUMMARY.tsv',sep='\t',index=False);(P/'qc/ANALYSIS_COMPLETE.json').write_text(json.dumps({'status':'COMPLETE','time':datetime.datetime.now().astimezone().isoformat(),'effect_rows':len(r),'protocol_sha256':hashlib.sha256(f.read_bytes()).hexdigest()},indent=2));print('SPATIAL_ANALYSIS_COMPLETE',len(r),flush=True)
