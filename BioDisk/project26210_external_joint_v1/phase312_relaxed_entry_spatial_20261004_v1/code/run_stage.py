from pathlib import Path
import subprocess,sys,datetime,json
P=Path(__file__).resolve().parents[1];name=sys.argv[1];assert name.isidentifier();code=P/'code'/(name+'.py');assert code.exists()
with (P/'logs'/(name+'.log')).open('w') as log:r=subprocess.run([sys.executable,'-u',str(code)],stdout=log,stderr=subprocess.STDOUT)
(P/'logs'/(name+'.exit')).write_text(str(r.returncode)+'\n');print('STAGE_EXIT',name,r.returncode,datetime.datetime.now().astimezone().isoformat(),flush=True);print((P/'logs'/(name+'.log')).read_text()[-5000:],flush=True);sys.exit(r.returncode)
