from pathlib import Path
import pandas as pd,numpy as np,json,subprocess
from scipy.spatial import cKDTree
P=Path(__file__).resolve().parents[1];B=P.parent/'phase315_CRC_SPP1_barrier_paper_reproduction_20261004_v1'
subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],check=True)
def read(p):return pd.read_csv(p,sep='\t',float_precision='round_trip')
co=read(P/'tables/CANDIDATE_SCALED_COORDINATES.tsv');m=read(P/'tables/EXAMPLE_METADATA.tsv');qc=read(P/'tables/CANDIDATE_COUNT_QC.tsv');counts=read(P/'tables/CANDIDATE_COUNTS_TARGETS.tsv');data=read(P/'tables/CANDIDATE_DATA_TARGETS.tsv')
orig=read(B/'tables/SOURCE_Figure_6e.tsv');checks=[]
def ck(n,ok,detail=None):checks.append({'check':n,'pass':bool(ok),'detail':detail})
for name,a in [('coordinates',co),('metadata',m),('QC',qc),('counts',counts),('data',data)]:ck(name+'_unique_barcodes',a.barcode.is_unique and len(a)==4248,len(a))
dist,ix=cKDTree(orig[['imagerow','imagecol']]).query(co[['source_imagerow','source_imagecol']])
ck('all_P7_coordinates_match',np.max(dist)<1e-4,{'max_distance':float(max(dist)),'matched':int((dist<1e-4).sum())});ck('one_to_one_P7_coordinates',len(np.unique(ix))==4248,len(np.unique(ix)))
x=co.merge(m[['barcode','DefineTypes']],on='barcode',validate='one_to_one');ck('all_original_labels_identical',np.array_equal(x.DefineTypes.to_numpy(),orig.iloc[ix].DefineTypes.to_numpy()),int((x.DefineTypes.to_numpy()==orig.iloc[ix].DefineTypes.to_numpy()).sum()))
sourceqcs=list(B.glob('tables/SOURCE_Supplementary_Figure_8*.tsv'));qc_candidates=[]
original_p7_qc=pd.read_excel(B/'downloads/source.xlsx',sheet_name='Figure 8d-f').rename(columns={'imageRow':'imagerow','imageCol':'imagecol'})
original_p7_qc.to_csv(P/'tables/ORIGINAL_P7_QC_Figure_8d-f.tsv',sep='\t',index=False)
sourceqcs.append(P/'tables/ORIGINAL_P7_QC_Figure_8d-f.tsv')
for f in sourceqcs:
 z=read(f)
 if 'Patient' in z and z.Patient.eq('P7').all():qc_candidates.append((f,z))
ck('one_original_P7_QC_sheet',len(qc_candidates)==1,[str(a[0].name) for a in qc_candidates])
if len(qc_candidates)==1:
 f,z=qc_candidates[0];dd,jj=cKDTree(z[['imagerow','imagecol']]).query(co[['source_imagerow','source_imagecol']]);ck('QC_coordinates_match',np.max(dd)<1e-4,float(max(dd)))
 qq=co[['barcode']].merge(qc,on='barcode',validate='one_to_one')
 for field,reference in [('meta_library_size','nCount_Spatial'),('meta_detected_features','nFeature_Spatial')]:
  delta=np.abs(qq[field].to_numpy()-z.iloc[jj][reference].to_numpy());ck('original_'+reference+'_identity',np.max(delta)==0,{'max_abs':float(max(delta)),'exact':int((delta==0).sum())})
for field,ref in [('counts_library_size','meta_library_size'),('counts_detected_features','meta_detected_features')]:ck(field+'_equals_meta',np.array_equal(qc[field],qc[ref]),float(max(abs(qc[field]-qc[ref]))))
v=counts.drop(columns='barcode').to_numpy();ck('target_counts_nonnegative_integer',np.isfinite(v).all() and (v>=0).all() and (v==np.floor(v)).all())
vc=data.drop(columns='barcode').to_numpy();ck('target_data_nonnegative_finite',np.isfinite(vc).all() and (vc>=0).all())
cc=counts.set_index('barcode');dd=data.set_index('barcode');tot=qc.set_index('barcode').counts_library_size.loc[cc.index].to_numpy();recon=np.log1p(cc.to_numpy()/tot[:,None]*1e4);delta=np.max(abs(recon-dd.loc[cc.index].to_numpy()));ck('stored_data_log1pCP10K',delta<1e-6,float(delta))
pd.DataFrame(checks).to_csv(P/'qc/P7_IDENTITY_CHECKS.tsv',sep='\t',index=False)
result={'status':'PASS' if all(a['pass'] for a in checks) else 'FAIL','checks':checks,'max_coordinate_error':float(max(dist)),'normalization_max_error':float(delta),'primary_effect_not_yet_computed':True}
(P/'qc/P7_IDENTITY_VALIDATION.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
if result['status']=='PASS':
 x['source_row']=ix;x.to_csv(P/'tables/P7_VERIFIED_COORDINATE_MAPPING.tsv',sep='\t',index=False)
