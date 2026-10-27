from pathlib import Path
import pandas as pd,subprocess
r=Path('D:/project26210_external_joint_v1')
for p in (r/'phase32_part1_m03_cnv_clone_trees_local_v1/input_snapshot_local_v1').glob('*'):
 if p.is_file():
  print(p.name,p.stat().st_size)
  if '.tsv' in p.name:print(pd.read_csv(p,sep='\t',nrows=0).columns.tolist())
for p in [r/'phase43_part1_client_figures_local_v1/source_remote_v1/tumor_epithelial_malignancy_calls_v1.tsv.gz',r/'phase43_part1_client_figures_local_v1/source_remote_v1/infercnv_query_cell_calls_v1.tsv.gz']:
 print(p.name,pd.read_csv(p,sep='\t',nrows=0).columns.tolist())
code="""from pathlib import Path
import pandas as pd
r=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1')
for folder in ['03_cnv','04_subclones_tymp','05_cnmf']:
 p=r/folder
 print('FOLDER',str(p))
 for f in p.glob('*'):print(f.name)
for p in r.glob('04_subclones_tymp/*/*'):
 if p.is_file() and ('tymp' in p.name.lower() or 'metadata' in p.name or 'landscape' in p.name):
  print(str(p))
  if '.tsv' in p.name:print(pd.read_csv(p,sep='\\t',nrows=0).columns.tolist())
"""
p=subprocess.run(['ssh','-T','-o','BatchMode=yes','-o','ConnectTimeout=20','sgs10-codex','/mnt/BioDisk/wangmingyuan/project26210_envs/part4_advanced_py_v1/bin/python -'],input=code.encode(),capture_output=True,timeout=120)
print(p.stdout.decode());print(p.stderr.decode());p.check_returncode()
