from pathlib import Path
import json,time,subprocess,datetime,traceback,os,zipfile,re,hashlib
P=Path(__file__).resolve().parents[1];PY='/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python'
def state(s,**kw):
 r=dict(status=s,time=datetime.datetime.now().astimezone().isoformat(),controller_pid=os.getpid(),**kw);(P/'STATUS.json').write_text(json.dumps(r,ensure_ascii=False,indent=2));print(s,flush=True)
def execute(script,log,*args):
 with (P/'logs'/f'{log}.log').open('ab') as f:r=subprocess.run([PY,str(P/'code'/script),*args],stdout=f,stderr=subprocess.STDOUT)
 (P/'logs'/f'{log}.exit').write_text(str(r.returncode));return r.returncode
try:
 pid=int((P/'logs/OBRADOVIC_download.pid').read_text());state('WAITING_OBRADOVIC_DIRECT_DOWNLOAD',download_pid=pid,monitor_interval_seconds=900,HCRN='queued_proxy_after_Obradovic_preparation; public endpoint currently HTTP403')
 # One persistent monitor checks completion every 15 minutes; no agent-level polling.
 while True:
  ex=P/'logs/OBRADOVIC_download.exit'
  if ex.exists():
   if ex.read_text().strip()!='0':raise RuntimeError('OBRADOVIC_DOWNLOAD_FAILED')
   v=json.loads((P/'qc/OBRADOVIC_DOWNLOAD_VERIFIED.json').read_text());j=json.loads((P/'resources/OBRADOVIC_DOWNLOAD_JOB.json').read_text());assert v['sha256']==j['sha256'] and v['bytes']==j['size'];assert (P/'downloads'/j['filename']).stat().st_size==j['size'];break
  proc=Path(f'/proc/{pid}/cmdline')
  if not proc.exists() or str(P/'code/download_direct_parallel.py').encode() not in proc.read_bytes():raise RuntimeError('Download process exited without successful completion record')
  time.sleep(900)
 state('PREPARING_OBRADOVIC');rc=execute('prepare_obradovic.py','OBRADOVIC_preparation')
 obr='READY_FOR_IDENTITY_REVIEW' if rc==0 and (P/'qc/OBRADOVIC_PREPARATION_COMPLETE.json').exists() else 'PREPARATION_FAILED_REVIEW_LOG'
 state('DOWNLOADING_HCRN_PROXY',Obradovic=obr);hr=execute('download_matrix.py','HCRN_download_wrapper',str(P/'resources/HCRN_DOWNLOAD_JOB.json'))
 if hr==0:
  assert (P/'logs/HCRN_download.exit').read_text().strip()=='0';entries=json.loads((P/'metadata/HCRN_ZIP_CONTENTS.json').read_text());names=[e['path'] for e in entries if '__MACOSX' not in e['path']];ids=sorted(set(re.findall(r'260-\d{4}', '\n'.join(names))))
  (P/'metadata/HCRN_ARCHIVE_PATIENT_IDENTIFIERS.json').write_text(json.dumps({'patient_ids':ids,'paths':names},indent=2));hstatus='DOWNLOADED_VERIFIED_NEEDS_SAMPLE_MAP_REVIEW'
 else:hstatus='DOWNLOAD_FAILED_SEE_LOG; public stream previously HTTP403, API HTTP401'
 state('INPUT_STAGE_REQUIRES_REVIEW',Obradovic=obr,HCRN=hstatus,analysis_complete=False);(P/'logs/controller.exit').write_text('0' if rc==0 and hr==0 else '2')
except BaseException:
 traceback.print_exc();state('PIPELINE_FAILED_REVIEW_LOG',analysis_complete=False);(P/'logs/controller.exit').write_text('1');raise
