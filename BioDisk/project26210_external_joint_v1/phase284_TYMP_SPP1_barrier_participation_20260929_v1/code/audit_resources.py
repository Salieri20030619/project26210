from common import *
import h5py,re
from anndata.io import read_elem
def main():
 f=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/03_malignancy_cnv/cnv_annotated_h5ad_v1/part1_scvi_integrated_cell_annotated_cnv_v1.h5ad')
 with h5py.File(f) as h:
  d=pd.DataFrame({k:np.asarray(read_elem(h['obs'][k])) for k in ['cell_type_broad_v1','study','curated_patient_id','sample_id','source_tissue','annotation_confidence_v1']})
 z=d[d.cell_type_broad_v1.eq('Fibroblast')&d.source_tissue.str.lower().eq('tumor')]
 assert len(z)==39 and z.curated_patient_id.nunique()==1
 tab(z.groupby(['study','curated_patient_id','sample_id','annotation_confidence_v1'],observed=True).size().reset_index(name='cells'),'FIBROBLAST_RAW_ANNOTATION_ASSET.tsv')
 st=f.stat()
 dump(P/'qc/FIBROBLAST_ASSET_LIMIT.json',dict(source=str(f),bytes=st.st_size,mtime_ns=st.st_mtime_ns,tumor_cells=len(z),independent_patients=z.curated_patient_id.nunique(),interpretation='no replicated CAF-response inference'))
 sf=P/'resources/GSE311076_samples.soft';s=sf.read_text();rows=[]
 for block in s.split('^SAMPLE = ')[1:]:
  gsm=block.splitlines()[0].strip()
  title=re.search(r'^!Sample_title = (.*)$',block,re.M).group(1);geno=re.search(r'^!Sample_characteristics_ch1 = genotype: (.*)$',block,re.M).group(1)
  rows.append(dict(GSM=gsm,title=title,genotype=geno,tissue_context='tumor' if 'tumor' in title.lower() else 'normal',identifiable_independent_replicates='not provided in this record'))
 assert len(rows)==4
 tab(rows,'GSE311076_SAMPLE_QUALIFICATION.tsv')
 q=pd.DataFrame(rows).groupby(['genotype','tissue_context']).size()
 assert len(q)==4 and (q==1).all()
 dump(P/'qc/FUNCTIONAL_RESOURCE_AUDIT.json',dict(status='PASS',GSE311076_metadata_sha256=sha(sf),submitted_libraries=4,one_library_per_genotype_tissue=True,biological_replicate_count='not assumed from library number',new_differential_tests=False))
 print('RESOURCE_AUDIT_PASS',flush=True)
if __name__=='__main__':main()
