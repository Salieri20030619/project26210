from pathlib import Path
import subprocess,concurrent.futures,json,zipfile
O=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase169_GSE175540_treatment_audit_20260916_v1')
urls={
'reanalysis_supp.zip':'https://mdpi-res.com/d_attachment/ijms/ijms-25-03887/article_deploy/ijms-25-03887-s001.zip',
'reanalysis_supp_alt.zip':'https://www.mdpi.com/article/10.3390/ijms25073887/s1',
'reanalysis_repo.json':'https://api.github.com/repos/SongyunLi21/TLS-predictive-model-2.0/git/trees/main?recursive=1'}
def get(it):
 name,u=it;r=subprocess.run(['curl','-L','-sS','--connect-timeout','8','--max-time','25','-w','%{http_code}',u,'-o',str(O/name)],capture_output=True,text=True,timeout=30)
 p=O/name;z=zipfile.is_zipfile(p) if p.exists() else False
 if z:
  with zipfile.ZipFile(p) as f:
   print('ZIP',name,f.namelist())
   dest=O/(name+'.extracted');dest.mkdir(exist_ok=True)
   for x in f.infolist():
    assert (dest/x.filename).resolve().is_relative_to(dest.resolve());f.extract(x,dest)
 return dict(file=name,url=u,http=r.stdout,exit=r.returncode,zip=z)
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:log=list(pool.map(get,urls.items()))
(O/'SECONDARY_FETCH_LOG.json').write_text(json.dumps(log,indent=2));print(json.dumps(log))
