from pathlib import Path
import subprocess,json,datetime,hashlib,tarfile,traceback,concurrent.futures,os
P=Path(__file__).resolve().parents[1]
def now():return datetime.datetime.now().astimezone().isoformat()
def dump(f,o):Path(f).write_text(json.dumps(o,ensure_ascii=False,indent=2))
def curl(name,url,limit=35):
 f=P/'sources'/name
 r=subprocess.run(['curl','-q','-fL','--connect-timeout','12','--max-time',str(limit),'-o',str(f),url],stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,timeout=limit+5)
 return dict(file=name,url=url,exit=r.returncode,bytes=f.stat().st_size if f.exists() else 0,error=r.stderr[-500:] if r.returncode else '')
try:
 urls=[('filelist.txt','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE294nnn/GSE294006/suppl/filelist.txt')]
 urls += [(x+'.soft.txt','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+x+'&targ=self&form=text&view=full') for x in ['GSM8896252','GSM8896253','GSM8896254']]
 with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:out=list(ex.map(lambda x:curl(*x),urls))
 dump(P/'METADATA_ACCESS.json',out)
 dump(P/'STATUS.json',dict(status='DOWNLOADING',at=now()))
 url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE294nnn/GSE294006/suppl/GSE294006_RAW.tar'
 f=P/'downloads/GSE294006_RAW.tar.part'
 with open(P/'logs/download.log','w') as log:
  r=subprocess.run(['curl','-q','-fL','--connect-timeout','20','--max-time','600','--retry','2','--retry-delay','15','-C','-','-o',str(f),url],stdout=log,stderr=subprocess.STDOUT,timeout=1900)
 assert r.returncode==0,'download exit '+str(r.returncode)
 final=f.with_suffix('');os.rename(f,final)
 members=[]
 with tarfile.open(final) as tar:
  for member in tar.getmembers():
   if not member.isfile():continue
   target=(P/'extracted'/member.name).resolve()
   assert target.is_relative_to((P/'extracted').resolve())
   target.parent.mkdir(parents=True,exist_ok=True)
   b=tar.extractfile(member).read();assert len(b)==member.size
   target.write_bytes(b)
   members.append(dict(member=member.name,bytes=member.size,sha256=hashlib.sha256(b).hexdigest()))
 dump(P/'EXTRACTED_MANIFEST.json',members)
 dump(P/'DOWNLOAD_VALIDATION.json',dict(status='PASS',at=now(),archive_bytes=final.stat().st_size,archive_sha256=hashlib.sha256(final.read_bytes()).hexdigest(),members=len(members),source=url))
 dump(P/'STATUS.json',dict(status='DOWNLOAD_VERIFIED_ANALYSIS_PENDING',at=now()));(P/'DOWNLOAD_EXIT_STATUS').write_text('0\n')
 print(json.dumps(members),flush=True)
except BaseException:
 (P/'logs/ERROR_download.txt').write_text(traceback.format_exc());dump(P/'STATUS.json',dict(status='DOWNLOAD_FAILED',at=now(),error=traceback.format_exc()));(P/'DOWNLOAD_EXIT_STATUS').write_text('1\n');raise

