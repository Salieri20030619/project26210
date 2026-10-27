from pathlib import Path
import json
import numpy as np,pandas as pd,loompy
from scipy import sparse
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/part2_rebuild_v1');I=R/'08_regulons/input_v1';O=R/'08_regulons/looms_v1';O.mkdir(exist_ok=False)
x=sparse.load_npz(I/'GENES_BY_CELLS_RAW.npz').tocsc();m=pd.read_csv(I/'CELL_METADATA.tsv.gz',sep='\t',dtype={'candidate':str},low_memory=False);genes=pd.read_csv(I/'GENES.tsv',sep='\t').gene.to_numpy();sel=m.selected_discovery.to_numpy(bool);assert x.shape==(len(genes),len(m))
for name,idx in [('all_cells',np.arange(len(m))),('discovery',np.flatnonzero(sel))]:
 mm=m.iloc[idx];loompy.create(str(O/(name+'.loom')),x[:,idx].astype(np.float32),{'Gene':genes.astype(object)},{'CellID':mm.cell_id.to_numpy(dtype=object),'patient':mm.curated_patient_id.to_numpy(dtype=object),'study':mm.study.to_numpy(dtype=object),'candidate':mm.candidate.to_numpy(dtype=object)})
(O/'STATUS.json').write_text(json.dumps({'status':'COMPLETE_LOOMS','cells':len(m),'discovery':int(sel.sum()),'genes':len(genes),'figures_created':False},indent=2));print('LOOMS_COMPLETE')
