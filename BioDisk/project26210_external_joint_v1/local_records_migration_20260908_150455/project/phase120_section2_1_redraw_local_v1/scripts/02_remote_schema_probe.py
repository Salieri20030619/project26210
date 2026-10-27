from pathlib import Path
import subprocess
code="""from pathlib import Path
import tempfile,subprocess,pandas as pd
r=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1')
for p in (r/'03_malignancy_cnv/cnv_annotated_h5ad_v1').glob('*'):print('H5AD',p.name)
for p in (r/'03_malignancy_cnv/malignancy_cnv_summary_v1').glob('*.tsv*'):print('TABLE',p.name)
for p in (r/'04_subclones_tymp/cnv_clone_tymp_cnmf_association_v1').glob('*.tsv*'):print('STATS',p.name)
rscript=Path('/mnt/BioDisk/wangmingyuan/project26210_envs/phase11_infercnv_conda_v2/bin/Rscript')
print('RSCRIPT',rscript.exists())
if rscript.exists():
 with tempfile.TemporaryDirectory(prefix='26210_cnv_probe_') as t:
  script=Path(t)/'probe.R'
  script.write_text('o <- readRDS(\"'+str(r/'03_malignancy_cnv/infercnv_v2/patients/GSE156632__RCC2/infercnv_run/preliminary.infercnv_obj')+'\"); print(slotNames(o)); print(dim(o@expr.data)); print(names(o@reference_grouped_cell_indices)); print(sapply(o@reference_grouped_cell_indices,length)); print(sapply(o@observation_grouped_cell_indices,length)); print(range(o@expr.data)); print(head(o@gene_order))')
  p=subprocess.run([str(rscript),'--vanilla',str(script)],capture_output=True,text=True,timeout=120);print(p.stdout);print(p.stderr)
"""
p=subprocess.run(['ssh','-T','-o','BatchMode=yes','sgs10-codex','/mnt/BioDisk/wangmingyuan/project26210_envs/part4_advanced_py_v1/bin/python -'],input=code.encode(),capture_output=True,timeout=180)
print(p.stdout.decode());print(p.stderr.decode());p.check_returncode()
