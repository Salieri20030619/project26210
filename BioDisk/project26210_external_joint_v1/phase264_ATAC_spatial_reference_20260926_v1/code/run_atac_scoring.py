from pathlib import Path
import os,subprocess,json,datetime,collections,traceback
P=Path(__file__).resolve().parents[1]
L=['Endothelial','Malignant','Myeloid','SMC_Pericyte']
env=os.environ.copy();env.pop('PYTHONPATH',None)
env.update(RGTDATA=str(P/'software/rgtdata'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',MPLBACKEND='Agg')
def run(cmd,name):
    with (P/'logs'/f'{name}.log').open('w') as log:
        log.write('COMMAND '+json.dumps(cmd)+'\n');log.flush()
        result=subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT)
        log.write(f'\nEXIT_STATUS={result.returncode}\n')
    (P/'qc'/f'{name.upper()}_STATUS.json').write_text(json.dumps(dict(exit_status=result.returncode,ended=datetime.datetime.now().astimezone().isoformat(),command=cmd),indent=2))
    if result.returncode:raise RuntimeError(f'{name} exited {result.returncode}')
try:
    run([str(P/'software/venv/bin/rgt-motifanalysis'),'matching','--organism','hg38','--motif-dbs',str(P/'reference/four_fixed_motifs'),'--output-location',str(P/'motif_matching'),'--input-files']+[str(P/'footprints'/f'{x}.bed') for x in L],'motif_matching')
    unique={};counts={}
    for x in L:
        file=P/'motif_matching'/f'{x}_mpbs.bed'
        rows=[s.strip().split('\t') for s in file.open() if s.strip() and not s.startswith('#')]
        counts[x]=len(rows)
        for r in rows:
            key=(r[0],int(r[1]),int(r[2]),r[3],r[5])
            if key not in unique or float(r[4])>float(unique[key][4]):unique[key]=r
    union=P/'motif_matching/UNION_UNIQUE_mpbs.bed'
    union.write_text(''.join('\t'.join(unique[k])+'\n' for k in sorted(unique)))
    support=dict(collections.Counter(k[3] for k in unique))
    (P/'qc/MOTIF_UNION_QC.json').write_text(json.dumps(dict(per_lineage_rows=counts,union_unique_sites=len(unique),motif_support=support,deduplication='chr,start,end,motif,strand; common site universe used for every lineage'),indent=2))
    if not unique:raise RuntimeError('No motif matches')
    run([str(P/'software/venv/bin/rgt-hint'),'differential','--organism=hg38','--bc','--nc=2','--mpbs-files='+str(union),'--reads-files='+','.join(str(P/'atac_chr1_pilot'/f'{x}_chr1_fragment_endpoints.bam') for x in L),'--conditions='+','.join(L),'--output-location='+str(P/'differential'),'--output-prefix=RCC100_chr1'],'differential')
    (P/'qc/ATAC_SCORING_BATCH_STATUS.json').write_text(json.dumps(dict(status='COMMANDS_FINISHED',exit_status=0,ended=datetime.datetime.now().astimezone().isoformat()),indent=2))
    print('ATAC_SCORING_EXIT_STATUS=0',flush=True)
except Exception:
    traceback.print_exc()
    (P/'qc/ATAC_SCORING_BATCH_STATUS.json').write_text(json.dumps(dict(status='FAILED',exit_status=1,ended=datetime.datetime.now().astimezone().isoformat()),indent=2))
    raise
