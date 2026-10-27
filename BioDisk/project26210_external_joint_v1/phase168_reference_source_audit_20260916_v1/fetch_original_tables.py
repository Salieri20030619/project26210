from pathlib import Path
import subprocess,concurrent.futures,json
O=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase168_reference_source_audit_20260916_v1')
urls={f'original_mmc{k}.xlsx':f'https://ars.els-cdn.com/content/image/1-s2.0-S1074761322000814-mmc{k}.xlsx' for k in [2,3]}
urls['original_article_hal_response']='https://inserm.hal.science/inserm-03653110v1/document'
def get(it):
 name,u=it
 r=subprocess.run(['curl','-L','-sS','--connect-timeout','8','--max-time','25','-w','%{http_code}',u,'-o',str(O/name)],capture_output=True,text=True,timeout=30)
 p=O/name
 return dict(file=name,url=u,http=r.stdout,exit=r.returncode,bytes=p.stat().st_size if p.exists() else 0,magic=p.read_bytes()[:8].hex() if p.exists() else '')
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:log=list(pool.map(get,urls.items()))
(O/'ORIGINAL_EXTRA_FETCH.json').write_text(json.dumps(log,indent=2));print(json.dumps(log))
