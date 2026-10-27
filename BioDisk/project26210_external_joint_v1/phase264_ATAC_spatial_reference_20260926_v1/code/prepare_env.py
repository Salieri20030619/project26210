from pathlib import Path
import requests,json,subprocess,sys,os
P=Path(__file__).resolve().parents[1];S=P/'software';S.mkdir(exist_ok=True)
r=requests.get('https://pypi.org/pypi/RGT/json',proxies={'https':'http://127.0.0.1:17897'},timeout=(10,30));r.raise_for_status()
j=r.json();(P/'metadata/RGT_PYPI.json').write_text(json.dumps({k:j['info'].get(k) for k in ['name','version','requires_python','requires_dist','project_urls']},indent=2));print('RGT_VERSION',j['info']['version'],flush=True)
subprocess.run([sys.executable,'-m','venv',str(S/'venv')],check=True)
cmd=[str(S/'venv/bin/python'),'-m','pip','install','--disable-pip-version-check','--proxy','http://127.0.0.1:17897','RGT=='+j['info']['version']]
rc=subprocess.run(cmd).returncode
(P/'ENV_INSTALL_EXIT_STATUS').write_text(str(rc)+'\n')
if rc:raise SystemExit(rc)
subprocess.run([str(S/'venv/bin/python'),'-m','pip','freeze'],stdout=(P/'software/requirements.lock.txt').open('w'),check=True)
print('RGT_ENV_READY',flush=True)

