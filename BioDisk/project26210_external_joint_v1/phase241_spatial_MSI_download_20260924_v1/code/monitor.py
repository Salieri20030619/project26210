import os,pathlib,json,datetime,select
P=pathlib.Path(__file__).resolve().parents[1]
pid=int((P/'RUNNER_PID').read_text())
def now():return datetime.datetime.now().astimezone().isoformat()
def snapshot(final=False):
 state=json.loads((P/'STATUS.json').read_text()) if (P/'STATUS.json').exists() else {}
 d={'recorded_at':now(),'runner_pid':pid,'state':state,'final':final,'interval_seconds':900}
 part=P/'downloads/spatial_metabolomics.zip.part'
 d['downloaded_partial_bytes']=part.stat().st_size if part.exists() else 0
 if final:
  d['exit_code']=(P/'EXIT_STATUS').read_text().strip() if (P/'EXIT_STATUS').exists() else None
  d['run_finished']=(P/'RUN_FINISHED').exists()
  d['validation']=json.loads((P/'FINAL_VALIDATION.json').read_text()) if (P/'FINAL_VALIDATION.json').exists() else {}
  d['final_acceptance']='PASS' if d['exit_code']=='0' and d['run_finished'] and d['validation'].get('status')=='PASS' and (P/'COMPLETE_VERIFIED').exists() else 'FAIL_OR_INCOMPLETE'
 (P/'monitoring/LATEST.json').write_text(json.dumps(d,ensure_ascii=False,indent=2))
 with open(P/'monitoring/HISTORY.jsonl','a') as f:f.write(json.dumps(d)+'\n')
 if final:
  (P/'monitoring/COMPLETION_VALIDATION.json').write_text(json.dumps(d,ensure_ascii=False,indent=2))
  (P/'monitoring/EXIT_STATUS').write_text('0\n' if d['final_acceptance']=='PASS' else '1\n')
snapshot()
try: fd=os.pidfd_open(pid)
except ProcessLookupError:snapshot(True);raise SystemExit(0)
while True:
 ready,_,_=select.select([fd],[],[],900)
 if ready:
  snapshot(True);break
 snapshot()
os.close(fd)
