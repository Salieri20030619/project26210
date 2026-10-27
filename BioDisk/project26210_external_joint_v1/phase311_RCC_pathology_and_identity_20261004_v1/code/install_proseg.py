from pathlib import Path
import subprocess,json,hashlib,tarfile,os,traceback
P=Path(__file__).resolve().parents[1];rc=1
try:
 meta=json.loads((P/'metadata/proseg_conda.json').read_text())
 entry=next(x for x in meta['files'] if x['basename']=='linux-64/rust-proseg-2.0.6-h4349ce8_0.tar.bz2')
 url='https://conda.anaconda.org/bioconda/'+entry['basename'];dest=P/'downloads/rust-proseg-2.0.6-h4349ce8_0.tar.bz2'
 subprocess.run(['curl','--proxy','http://127.0.0.1:7897','--retry','2','--retry-all-errors','--connect-timeout','15','--max-time','300','-fLsS','-o',str(dest),url],check=True)
 assert dest.stat().st_size==entry['size']
 assert hashlib.sha256(dest.read_bytes()).hexdigest()==entry['sha256']
 tool=P/'resources/proseg206';tool.mkdir()
 with tarfile.open(dest,'r:bz2') as tar:
  for m in tar.getmembers():
   target=(tool/m.name).resolve();assert target.is_relative_to(tool.resolve())
   if m.isfile():
    target.parent.mkdir(parents=True,exist_ok=True)
    with tar.extractfile(m) as f:target.write_bytes(f.read())
    target.chmod(m.mode)
   elif m.isdir():target.mkdir(exist_ok=True,parents=True)
   else:raise ValueError('Unexpected archive link')
 exe=tool/'bin/proseg'
 ver=subprocess.check_output([str(exe),'--version'],text=True).strip()
 helptext=subprocess.check_output([str(exe),'--help'],text=True)
 (P/'metadata/proseg206_help.txt').write_text(helptext)
 rec={'status':'PASS','version':ver,'url':url,'package_sha256':entry['sha256'],'binary_sha256':hashlib.sha256(exe.read_bytes()).hexdigest(),'choice':'Pinned small self-contained Bioconda build; no global installation, existing environments unchanged.'}
 (P/'qc/PROSEG_INSTALL.json').write_text(json.dumps(rec,indent=2));print(rec,flush=True);print(helptext,flush=True);rc=0
except Exception:traceback.print_exc()
finally:(P/'logs/install_proseg.exit').write_text(str(rc)+'\n')
raise SystemExit(rc)
