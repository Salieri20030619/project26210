from pathlib import Path
import pandas as pd,h5py,json
from anndata.io import read_elem
O=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase145_macrophage_provenance_trace_remote_v1')
a=pd.read_csv(O/'legacy_object/metadata.tsv.gz',sep='\t',low_memory=False)
p='/mnt/BioDisk/ludingyu/project26210/object/01_annotation_adata_annotated.h5ad'
with h5py.File(p,'r') as h:b=read_elem(h['obs'])
print('OLD_PARENT',len(b),list(b.columns))
inter=set(a.cell_id)&set(b.index);c=b.loc[list(inter)]
rows={'legacy_h5ad':p,'parent_cells':len(b),'subtype_cells':len(a),'exact_parent_overlap':len(inter),'old_parent_cell_type_counts':b.cell_type.value_counts().to_dict() if 'cell_type' in b else {},'parent_macrophage_not_retained':int(((b.cell_type=='Macrophage')&~b.index.isin(a.cell_id)).sum()) if 'cell_type' in b else None}
if 'cell_type' in b:
 rows['old_subtype_parent_celltype']=c.cell_type.value_counts().to_dict()
 m=b[(b.cell_type=='Macrophage')&~b.index.isin(a.cell_id)];m.to_csv(O/'identity/OLD_PARENT_MACRO_NOT_IN_FINAL.tsv.gz',sep='\t',compression='gzip')
(O/'identity/LEGACY_PARENT_STATUS.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2));print(json.dumps(rows,ensure_ascii=False))
