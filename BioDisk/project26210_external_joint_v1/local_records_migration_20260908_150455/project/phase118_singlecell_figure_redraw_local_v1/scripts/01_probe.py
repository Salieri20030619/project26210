from pathlib import Path
import subprocess
root=Path(__file__).resolve().parents[1]
code='''import anndata as ad, json
p='/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/02_integration_annotation/cell_annotation_v1/part1_scvi_integrated_cell_annotated_v1.h5ad'
a=ad.read_h5ad(p,backed='r')
print(json.dumps({'shape':a.shape,'obs':list(a.obs.columns),'obsm':list(a.obsm),'layers':list(a.layers),'X':str(a.X),'uns':list(a.uns)},default=str))
'''
r=subprocess.run(['ssh','-T','-o','BatchMode=yes','-o','ConnectTimeout=20','sgs10-codex','/mnt/BioDisk/wangmingyuan/project26210_envs/part4_advanced_py_v1/bin/python -'],input=code.encode(),capture_output=True,timeout=120)
(root/'qc'/'source_probe.json').write_bytes(r.stdout)
print(r.stdout.decode());print(r.stderr.decode());r.check_returncode()
