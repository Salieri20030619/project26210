import concurrent.futures, datetime, json, os, subprocess
from pathlib import Path
P=Path(__file__).resolve().parents[1]
env=os.environ.copy(); env.pop('PYTHONPATH',None)
env.update(RGTDATA=str(P/'software/rgtdata'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
lineages=['Endothelial','Malignant','Myeloid','SMC_Pericyte']
out=P/'footprints'; out.mkdir(exist_ok=True)
def run(lineage):
    cmd=[str(P/'software/venv/bin/rgt-hint'),'footprinting','--organism=hg38','--atac-seq','--paired-end','--output-location='+str(out),'--output-prefix='+lineage,str(P/'atac_chr1_pilot'/f'{lineage}_chr1_fragment_endpoints.bam'),str(P/'tables/RCC100_chr1_all_peaks.bed')]
    start=datetime.datetime.now().astimezone().isoformat()
    with (P/'logs'/f'footprint_{lineage}.log').open('w') as log:
        log.write('COMMAND '+json.dumps(cmd)+'\n'); log.flush()
        rc=subprocess.call(cmd,env=env,stdout=log,stderr=subprocess.STDOUT)
        log.write(f'\nFOOTPRINT_EXIT_STATUS={rc}\n')
    record=dict(lineage=lineage,exit_status=rc,started=start,ended=datetime.datetime.now().astimezone().isoformat(),command=cmd)
    (P/'qc'/f'FOOTPRINT_{lineage}_STATUS.json').write_text(json.dumps(record,indent=2))
    return record
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
    results=list(pool.map(run,lineages))
status=dict(status='COMMANDS_FINISHED' if all(r['exit_status']==0 for r in results) else 'COMMAND_FAILED',results=results)
(P/'qc/FOOTPRINT_BATCH_STATUS.json').write_text(json.dumps(status,indent=2))
rc=0 if all(r['exit_status']==0 for r in results) else 1
print('FOOTPRINT_BATCH_EXIT_STATUS='+str(rc),flush=True)
raise SystemExit(rc)
