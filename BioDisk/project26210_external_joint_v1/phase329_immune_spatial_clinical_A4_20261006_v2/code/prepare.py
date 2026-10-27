from pathlib import Path
import pandas as pd,numpy as np,h5py,json,hashlib,shutil
P=Path(__file__).resolve().parents[1];R=P.parent
sources=[]
def copy(f,name):
 sources.append({'path':str(f),'sha256':hashlib.sha256(f.read_bytes()).hexdigest()});shutil.copy2(f,P/'data'/name)
for f in (R/'phase328_immune_clinical_A4_20261006_v1/data').glob('*.tsv'):copy(f,f.name)
for phase,files,prefix in [('phase317_RCC_external_paper_relation_validation_20261005_v1',['PATIENT_EFFECTS_ALL.tsv','COHORT_SUMMARY_ALL.tsv'],'SPATIAL_'),('phase275_immune_interface_tests_20260928_v1',['PATIENT_EFFECTS.tsv','COHORT_SUMMARY.tsv'],'BOUNDARY_'),('phase319_SPP1_direct_barrier_paper_reproduction_20261005_v1',['BIOLOGICAL_OR_DECLARED_UNIT_SUMMARY.tsv','QUALIFIED_MODEL_CHECKS.tsv','ALL_GROUP_SUMMARIES.tsv'],'KO_')]:
 for f in files:copy(R/phase/'tables'/f,prefix+f)
S=R/'phase317_RCC_external_paper_relation_validation_20261005_v1';reg=pd.read_csv(S/'resources/EXTERNAL_SECTION_REGISTER.tsv',sep='\t');r=reg[reg['sample']=='GSM5924033'].iloc[0]
for f in ['METADATA.tsv','MODULE_SCORES.tsv']:copy(S/'inputs/GSM5924033'/f,'GSM5924033_'+f)
z=pd.read_csv(P/'data/GSM5924033_MODULE_SCORES.tsv',sep='\t')
with h5py.File(r.input_path) as h:
 ids=h['obs/_index'].asstr()[:];xy=h['obsm/spatial'][:]
assert list(ids)==z.spot_id.tolist();z['x_um']=xy[:,0]*r.pixel_to_um;z['y_um']=xy[:,1]*r.pixel_to_um
z.to_csv(P/'data/SPATIAL_MAP_VALUES.tsv',sep='\t',index=False)
sources.append({'path':r.input_path,'datasets':['obs/_index','obsm/spatial'],'coordinate_sha256':hashlib.sha256(xy.tobytes()).hexdigest(),'sample_selection':'Same GSM5924033 section shown in user-specified reference Figure 3; selected before viewing score maps, not by effect size.','pixel_to_um':r.pixel_to_um})
(P/'qc/SOURCE_MANIFEST.json').write_text(json.dumps(sources,indent=2));print('PREPARED',len(z),'spatial spots')
