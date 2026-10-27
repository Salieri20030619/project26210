from pathlib import Path
import subprocess,sys,json,datetime
P=Path(__file__).resolve().parents[1];name=sys.argv[1]
(P/'STATUS.json').write_text(json.dumps({'status':'RUNNING','stage':name,'time':datetime.datetime.now().astimezone().isoformat()},indent=2))
with (P/'logs'/f'{name}.log').open('w') as f:q=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python','-u',str(P/'code'/f'{name}.py')],stdout=f,stderr=subprocess.STDOUT)
(P/'logs'/f'{name}.exit').write_text(str(q.returncode));print(name,'EXIT',q.returncode);print((P/'logs'/f'{name}.log').read_text()[-6500:]);sys.exit(q.returncode)
