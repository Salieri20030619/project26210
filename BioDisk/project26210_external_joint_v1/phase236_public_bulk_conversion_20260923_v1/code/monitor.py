from common import *
import select,traceback,csv
pid=int((P/'RUNNER_PID').read_text());root=P.parent
try:fd=os.pidfd_open(pid)
except ProcessLookupError:fd=None
def snapshot(alive):
 state=json.loads((P/'STATUS.json').read_text()) if (P/'STATUS.json').exists() else {}
 o=dict(recorded_at=now(),runner_pid=pid,runner_present=alive,interval_seconds=900,state=state,exit_code=(P/'EXIT_STATUS').read_text().strip() if (P/'EXIT_STATUS').exists() else None)
 save(P/'monitoring/LATEST.json',o)
 with (P/'monitoring/HISTORY.jsonl').open('a') as f:f.write(json.dumps(o)+'\n')
 return o
snapshot(fd is not None)
while fd is not None:
 ready,_,_=select.select([fd],[],[],900)
 if ready:os.close(fd);fd=None;break
 snapshot(True)
result=snapshot(False);passed=False
try:
 assert result['exit_code']=='0' and (P/'RUN_FINISHED').exists() and (P/'COMPLETE_VERIFIED').exists()
 qa=json.loads((P/'MATRIX_QA.json').read_text());assert qa['status']=='PASS' and qa['sample_libraries']==60
 for item in qa['outputs']:
  q=pathlib.Path(item['path']);assert q.exists() and sha(q)==item['sha256']
 assert len(qa['outputs'])==9
 passed=True
except BaseException:result['acceptance_error']=traceback.format_exc()
result['final_acceptance']='PASS' if passed else 'INCOMPLETE_REVIEW_REQUIRED';save(P/'monitoring/COMPLETION_VALIDATION.json',result)
(P/'monitoring/EXIT_STATUS').write_text('0\n' if passed else '1\n')
text='60 bulk RNA-seq libraries converted; raw fragment counts, TPM, FPKM and cohort-specific QC/matrices validated.' if passed else 'Background conversion ended with incomplete acceptance; preserve all outputs and inspect recorded errors before resuming.'
note='\n\n## '+now()+' Phase236 bulk conversion '+result['final_acceptance']+'\n'+text+'\nEvidence: '+str(P/'monitoring/COMPLETION_VALIDATION.json')+'; scope is technical conversion only, no treatment/clinical comparisons and no donor count inference.\n'
for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
 q=root/name;q.write_text(note+q.read_text())
with (root/'RESULTS_MANIFEST.tsv').open('a') as f:f.write('\n'+'\t'.join(['phase236_public_bulk_conversion_final','bulk_FASTQ_conversion','json',str(P/'monitoring/COMPLETION_VALIDATION.json'),now(),result['final_acceptance'],'60 libraries expected; cohorts kept separate','No biological comparisons; original FASTQ and BAM retained'])+'\n')
