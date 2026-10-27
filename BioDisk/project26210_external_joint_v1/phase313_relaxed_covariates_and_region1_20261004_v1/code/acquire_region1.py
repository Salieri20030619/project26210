from pathlib import Path
import subprocess,json,hashlib,os,traceback,datetime
P=Path(__file__).resolve().parents[1];rc=1
try:
 url='https://api.figshare.com/v2/articles/25685961/files/45867546';meta=P/'metadata/region1_file.json'
 subprocess.run(['curl','--proxy','http://127.0.0.1:7897','--retry','2','--retry-all-errors','--connect-timeout','20','--max-time','60','-fLsS','-o',str(meta),url],check=True)
 info=json.loads(meta.read_text());assert info['size']==330540918 and info['computed_md5']=='5eb1e20b77592b173643da3ce887869a'
 dest=P/'downloads/tumor1_region1_transcripts.csv.gz';part=Path(str(dest)+'.partial')
 print('DOWNLOAD_START',datetime.datetime.now().astimezone().isoformat(),info['id'],info['size'],flush=True)
 cmd=['curl','--proxy','http://127.0.0.1:7897','--retry','3','--retry-all-errors','--connect-timeout','20','--max-time','1800','-fLsS','-C','-','-o',str(part),'https://api.figshare.com/v2/file/download/45867546']
 subprocess.run(cmd,check=True);assert part.stat().st_size==info['size']
 md5=hashlib.md5();sha=hashlib.sha256()
 with part.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):md5.update(b);sha.update(b)
 assert md5.hexdigest()==info['computed_md5'];os.replace(part,dest)
 rec={'status':'PASS','source':url,'download_url':'https://api.figshare.com/v2/file/download/45867546','file':str(dest),'bytes':dest.stat().st_size,'md5':md5.hexdigest(),'sha256':sha.hexdigest(),'route':'proxy; <2GB','completed':datetime.datetime.now().astimezone().isoformat()}
 (P/'qc/REGION1_DOWNLOAD_VALIDATION.json').write_text(json.dumps(rec,indent=2));print(json.dumps(rec,indent=2),flush=True);rc=0
except Exception:traceback.print_exc()
finally:(P/'logs/acquire_region1.exit').write_text(str(rc)+'\n')
raise SystemExit(rc)
