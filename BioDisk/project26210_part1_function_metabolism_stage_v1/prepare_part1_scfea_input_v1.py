#!/usr/bin/env python3
from pathlib import Path
import gzip, json
import pandas as pd
from scipy.io import mmread

project=Path("/mnt/BioDisk/ludingyu/project26210")
root=project/"results/part1_rebuild_v1"
inp=root/"07_pyscenic_v1/input"
out=root/"09_metabolism_trajectory/metabolism_v1"
source=project/"envs/part2_6_2_7_20260726_v1/sources/scFEA/data/module_gene_m168.csv"
(out/"scFEA_input").mkdir(parents=True,exist_ok=True); (out/"qc").mkdir(parents=True,exist_ok=True)
genes=pd.read_csv(inp/"genes_v1.tsv",sep="\t")["gene"].astype(str).tolist()
cells=pd.read_csv(inp/"cells_v1.tsv",sep="\t")["exact_cell_id"].astype(str).tolist()
modules=pd.read_csv(source,index_col=0)
needed=set(modules.stack().dropna().astype(str)); selected=[i for i,g in enumerate(genes) if g in needed]
with gzip.open(inp/"raw_counts_genes_by_cells_v1.mtx.gz","rb") as handle: matrix=mmread(handle).tocsr()
assert matrix.shape==(len(genes),len(cells)) and matrix.nnz==35713877 and len(cells)==11950
sub=matrix[selected,:]
frame=pd.DataFrame.sparse.from_spmatrix(sub,index=[genes[i] for i in selected],columns=cells)
dest=out/"scFEA_input/malignant_all_cells_scFEA_M168_raw_counts.csv"
frame.to_csv(dest)
qc={"schema_version":"part1_scfea_input_v1","n_cells":len(cells),"n_input_genes":len(genes),"n_selected_model_genes":len(selected),"n_nonzero":int(sub.nnz),"raw_integer_nonnegative":bool(sub.data.min()>=0 and (sub.data==sub.data.astype(int)).all())}
(out/"qc/scFEA_input_qc_v1.json").write_text(json.dumps(qc,indent=2),encoding="utf-8")
(out/"SCFEA_INPUT_COMPLETE").write_text("PASS\n",encoding="utf-8")
print(json.dumps(qc,indent=2))
