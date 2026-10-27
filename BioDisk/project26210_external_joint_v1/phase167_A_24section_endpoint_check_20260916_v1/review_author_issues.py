from pathlib import Path
import json,os,signal,requests
O=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase167_A_24section_endpoint_check_20260916_v1')
target=str(O/'fetch_patient_sources.py')
stopped=[]
for p in Path('/proc').iterdir():
 if not p.name.isdigit():continue
 try:
  args=(p/'cmdline').read_bytes().split(b'\0')
  if target.encode() in args:
   os.kill(int(p.name),signal.SIGTERM);stopped.append(int(p.name))
 except (FileNotFoundError,PermissionError,ProcessLookupError):pass
(O/'RAW_DOWNLOAD_SUPERSEDED.json').write_text(json.dumps(dict(stopped_pids=stopped,reason='Slow raw.githubusercontent.com transfer; all 8 requested files obtained through official GitHub blob API at frozen commit'),indent=2))
issues=json.loads((O/'metadata_api/issues.json').read_text())
for x in issues:
 print('ISSUE',x['number'],x['title'],x.get('body'))
 if x['comments']:
  r=requests.get(x['comments_url'],timeout=15);r.raise_for_status()
  (O/'metadata_api'/('issue_'+str(x['number'])+'_comments.json')).write_bytes(r.content)
  for y in r.json():print('COMMENT',y['user']['login'],y['body'])
