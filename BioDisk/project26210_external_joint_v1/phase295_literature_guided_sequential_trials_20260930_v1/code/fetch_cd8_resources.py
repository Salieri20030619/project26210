from pathlib import Path
import json,subprocess,hashlib
P=Path('phase295_literature_guided_sequential_trials_20260930_v1');j=json.loads((P/'resources/cd8ippred_tree.json').read_text());out=P/'resources/cd8ippred';out.mkdir(exist_ok=True);audit=[]
for f in j['tree']:
 n=f['path']
 if f['type']!='blob' or not(n.startswith('R/') or n.startswith('inst/extdata/') or n in ['DESCRIPTION','vignettes/cd8_ip_predict.Rmd']):continue
 dest=out/n;dest.parent.mkdir(exist_ok=True,parents=True);url='https://raw.githubusercontent.com/bedapub/cd8ippred/'+j['sha']+'/'+n
 a=subprocess.run(['curl','-L','--connect-timeout','10','--max-time','60','-fsS',url,'-o',str(dest)],capture_output=True,text=True)
 if a.returncode: audit.append(dict(file=n,status='DOWNLOAD_FAILED',error=a.stderr));continue
 b=dest.read_bytes();actual=hashlib.sha1(('blob '+str(len(b))+'\0').encode()+b).hexdigest();assert actual==f['sha']
 audit.append(dict(file=n,status='VERIFIED',bytes=len(b),git_blob=actual,sha256=hashlib.sha256(b).hexdigest(),url=url))
(P/'resources/CD8_RESOURCE_AUDIT.json').write_text(json.dumps(audit,indent=2));print('resources',len(audit),flush=True)
