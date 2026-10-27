from pathlib import Path
import subprocess,json,datetime,sys
P=Path(__file__).resolve().parents[1]
commands=[('prepare',['/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python',str(P/'code/prepare_external_inputs.py')]),('score',['/mnt/BioDisk/wangmingyuan/miniconda3/envs/r_seurat/bin/Rscript','--vanilla',str(P/'code/score_external_sections.R')]),('analyze',['/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python',str(P/'code/analyze_external_relations.py')])]
for name,cmd in commands:
 (P/'STATUS.json').write_text(json.dumps({'status':'RUNNING','stage':name,'time':datetime.datetime.now().astimezone().isoformat(),'reference_exact_model_calibration_pass':False,'evidence':'exploratory published-marker external transfer'},indent=2))
 print('START',name,flush=True)
 with (P/'logs'/f'{name}.log').open('w') as f:q=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
 (P/'logs'/f'{name}.exit').write_text(str(q.returncode));print('END',name,q.returncode,flush=True)
 if q.returncode:
  (P/'STATUS.json').write_text(json.dumps({'status':'TECHNICAL_FAILURE','stage':name,'exit':q.returncode},indent=2));sys.exit(q.returncode)
(P/'STATUS.json').write_text(json.dumps({'status':'ANALYSIS_DONE_VALIDATION_PENDING','time':datetime.datetime.now().astimezone().isoformat()},indent=2))
