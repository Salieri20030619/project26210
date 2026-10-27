from pathlib import Path
import json,subprocess,hashlib,datetime,zipfile
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase235_SpaMTP_TJRCC_pilot_20260921_v1')
files=[f for f in json.loads((P/'metadata/figshare_article_24599295_cached.json').read_text())['files'] if f['name'] in ['spatial_metabolites_annotation.xlsx','spatial_metabolomics.zip']]
results=[]
for f in files:
 dest=P/'downloads'/f['name'];part=dest.with_suffix(dest.suffix+'.part')
 with (P/'logs'/('download_'+str(f['id'])+'.log')).open('a') as log:
  r=subprocess.run(['curl','--fail','--location','--connect-timeout','10','--max-time','240','--retry','1','--output',str(part),f['download_url']],stdout=log,stderr=subprocess.STDOUT,timeout=520)
 row=dict(name=f['name'],url=f['download_url'],exit_code=r.returncode,expected_bytes=f['size'],expected_md5=f['computed_md5'],actual_bytes=part.stat().st_size if part.exists() else 0)
 if r.returncode==0:
  row['md5']=hashlib.md5(part.read_bytes()).hexdigest();row['verified']=row['md5']==f['computed_md5'] and row['actual_bytes']==f['size']
  if row['verified']:part.rename(dest)
 else:row['verified']=False
 results.append(row);(P/'DOWNLOAD_STATUS.json').write_text(json.dumps(results,indent=2));print(json.dumps(row),flush=True)
 if not row['verified']:print((P/'logs'/('download_'+str(f['id'])+'.log')).read_text()[-1200:],flush=True);break
(P/'DOWNLOAD_EXIT_STATUS').write_text('0\n' if len(results)==2 and all(r['verified'] for r in results) else '1\n')
