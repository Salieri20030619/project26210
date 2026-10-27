from pathlib import Path
import os,sys,subprocess,json,datetime,hashlib,shutil,time,signal
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/part2_rebuild_v1')
O=R/'08_regulons/network_v2_20260915'
Q=Path('/mnt/BioDisk/ludingyu/project26210')
PY=Q/'.venv-pyscenic/bin/python'
WRAP=Q/'run_pyscenic_compat.py'
REF=Q/'references/pyscenic/hg38'
GATE=Q/'check_memory_gate.sh'
accept=R/'provenance/pre_scenic_acceptance_20260915_v2/STATUS.json'
assert json.loads(accept.read_text())['status']=='PASS'
def stamp():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def save(stage,status,**kw):
 obj=dict(time_utc=stamp(),status=status,stage=stage,pid=os.getpid(),output=str(O),acceptance=str(accept),figures_created=False,**kw)
 tmp=O/'STATUS.tmp';tmp.write_text(json.dumps(obj,indent=2));tmp.replace(O/'STATUS.json')
 print(json.dumps(obj),flush=True)
def gate():
 with (O/'RESOURCE_GATE.log').open('a') as f:subprocess.run([str(GATE)],stdout=f,stderr=subprocess.STDOUT,check=True)
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
if '--launch' in sys.argv:
 O.mkdir(exist_ok=False)
 with (O/'launcher.log').open('ab') as f:
  p=subprocess.Popen([sys.executable,str(Path(__file__).resolve()),'--run'],stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,start_new_session=True,close_fds=True)
 (O/'PID.txt').write_text(str(p.pid)+'\n')
 print('LAUNCHED',p.pid,str(O));sys.exit(0)
env=os.environ.copy()
env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONUNBUFFERED='1')
stage='preflight'
try:
 save(stage,'RUNNING');gate()
 import numpy as np
 for a,b in {'object':object,'float':float,'int':int,'bool':bool}.items():
  if a not in np.__dict__:setattr(np,a,b)
 import loompy,pandas as pd
 paths=[R/'08_regulons/looms_v1/discovery.loom',R/'08_regulons/looms_v1/all_cells.loom',REF/'allTFs_hg38.txt',REF/'hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather',REF/'hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather',REF/'motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl',WRAP,accept]
 for p in paths:assert p.is_file() and p.stat().st_size>0,p
 with loompy.connect(str(paths[0]),mode='r') as a,loompy.connect(str(paths[1]),mode='r') as b:
  assert a.shape==(9479,14210) and b.shape==(9479,33098)
  assert len(set(a.ca.CellID))==14210 and len(set(b.ca.CellID))==33098
  assert set(a.ca.CellID)<=set(b.ca.CellID) and np.array_equal(a.ra.Gene,b.ra.Gene)
  genes=set(map(str,a.ra.Gene));ids=set(map(str,b.ca.CellID))
  md=pd.read_csv(R/'08_regulons/input_v1/CELL_METADATA.tsv.gz',sep='\t',index_col=0)
  assert ids==set(map(str,md.index))
  tf=set(paths[2].read_text().splitlines())&genes;assert len(tf)>100
 manifest=[dict(path=str(p),bytes=p.stat().st_size,sha256=digest(p)) for p in paths]
 (O/'INPUT_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
 (O/'PREFLIGHT.json').write_text(json.dumps(dict(status='PASS',genes=9479,discovery_cells=14210,all_cells=33098,TF_overlap=len(tf),workers=8,seed=262148,python=sys.version),indent=2))
 base=[str(PY),str(WRAP)]
 commands=[
 ('GRNBoost2',base+['grn',str(paths[0]),str(paths[2]),'--method','grnboost2','--num_workers','8','--seed','262148','--cell_id_attribute','CellID','--gene_attribute','Gene','--sparse','-o',str(O/'adjacencies.tsv')]),
 ('motif_pruning',base+['ctx',str(O/'adjacencies.tsv'),str(paths[3]),str(paths[4]),'--annotations_fname',str(paths[5]),'--expression_mtx_fname',str(paths[0]),'--num_workers','8','--cell_id_attribute','CellID','--gene_attribute','Gene','-o',str(O/'regulons.dat')]),
 ('AUCell',base+['aucell',str(paths[1]),str(O/'regulons.dat'),'--num_workers','8','--seed','262148','--cell_id_attribute','CellID','--gene_attribute','Gene','-o',str(O/'AUC_ALL_CELLS.csv')])]
 (O/'COMMANDS.json').write_text(json.dumps(commands,indent=2))
 for stage,cmd in commands:
  gate();save(stage,'RUNNING')
  with (O/(stage+'.log')).open('wb') as f:
   p=subprocess.Popen(cmd,stdout=f,stderr=subprocess.STDOUT,env=env,stdin=subprocess.DEVNULL)
   (O/'CHILD_PID.txt').write_text(str(p.pid)+'\n')
   code=p.wait()
  (O/(stage+'_EXIT.json')).write_text(json.dumps(dict(stage=stage,returncode=code,time_utc=stamp())))
  if code!=0:raise RuntimeError(stage+' returncode='+str(code))
  gate()
  if stage=='GRNBoost2':
   d=pd.read_csv(O/'adjacencies.tsv',sep='\t')
   assert len(d)>0 and {'TF','target','importance'}<=set(d)
   assert np.isfinite(d.importance).all() and set(d.TF)<=tf and set(d.target)<=genes
   assert not d.duplicated(['TF','target']).any()
  elif stage=='motif_pruning':
   assert (O/'regulons.dat').stat().st_size>0
   import pickle
   with (O/'regulons.dat').open('rb') as f:regs=pickle.load(f)
   assert len(regs)>0 and all(len(r.genes)>0 for r in regs)
  else:
   d=pd.read_csv(O/'AUC_ALL_CELLS.csv',index_col=0)
   assert d.index.is_unique and set(d.index)==ids and d.shape[1]>0 and np.isfinite(d.values).all()
   assert (d.values>=0).all() and (d.values<=1+1e-6).all()
  (O/(stage+'_VERIFIED.json')).write_text(json.dumps(dict(status='PASS',time_utc=stamp())))
 save('complete','COMPLETE_SCENIC_CORE_VERIFIED',next_action='Patient-level regulon integration and independent biological review remain; no figures.')
 (O/'RUN_COMPLETE').write_text(stamp())
except BaseException as e:
 save(stage,'FAILED',error=repr(e))
 raise

