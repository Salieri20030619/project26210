from pathlib import Path
import json,datetime,hashlib,subprocess
import numpy as np,pandas as pd
from scipy import io,stats
P=Path(__file__).resolve().parents[1];ROOT=P.parent;EXT=ROOT/'phase247_external_spatial_cases_20260924_v1'
def now():return datetime.datetime.now().astimezone().isoformat()
def sha(f):return hashlib.sha256(Path(f).read_bytes()).hexdigest()
assert (P/'EXIT_STATUS').read_text().strip()=='0'
plan=json.loads((P/'PLAN_BEFORE_ANALYSIS.json').read_text());ef=pd.read_csv(P/'tables/EXTERNAL_CASE_MAINLINE_EFFECTS.tsv',sep='\t')
check=[];scorechecks=[]
for sample in plan['external_cases']:
 r=subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True)
 with open(P/'logs/validation_memory.log','a') as f:f.write(sample+'\n'+r.stdout+r.stderr)
 assert r.returncode==0
 d=pd.read_csv(P/'inputs'/f'{sample}_MAINLINE_INPUTS.tsv.gz',sep='\t')
 info=pd.read_csv(EXT/'tables/EXTERNAL_INPUT_QC.tsv',sep='\t').set_index('sample').loc[sample]
 prefix=info.accession+'_'+sample+'_ST.'
 mat=io.mmread(EXT/'extracted'/(prefix+'matrix.mtx.gz')).tocsr().astype(float)
 genes=pd.read_csv(EXT/'extracted'/(prefix+'features.tsv.gz'),sep='\t',header=None)[1]
 b=pd.read_csv(EXT/'extracted'/(prefix+'barcodes.tsv.gz'),sep='\t',header=None)[0]
 ix=pd.Index(b).get_indexer(d.barcode);assert (ix>=0).all();mat=mat[:,ix].T.tocsr();tot=np.asarray(mat.sum(1)).ravel()
 for name,gs in plan['programs'].items():
  used=[g for g in gs if genes.eq(g).any()]
  raw=np.column_stack([np.asarray(mat[:,np.flatnonzero(genes.eq(g))].sum(1)).ravel() for g in used])
  z=stats.zscore(np.log1p(raw/tot[:,None]*10000),axis=0,ddof=0,nan_policy='omit')
  score=np.nan_to_num(z,nan=0).mean(1)
  if d[name].notna().all():
   delta=float(np.max(np.abs(score-d[name])));assert delta<1e-10
   scorechecks.append(dict(sample=sample,program=name,max_difference=delta))
 for row in ef[(ef['sample']==sample)&(ef.model=='A0')&ef.scale.isin(['whole_tissue','block4x4'])].itertuples():
  g=d[[row.exposure,row.target]]
  if row.scale=='block4x4':
   g=d.groupby('block4')[[row.exposure,row.target]].median();g=g[d.groupby('block4').size()>=10]
  good=g.dropna()
  if pd.notna(row.effect):
   rr=stats.spearmanr(good.iloc[:,0],good.iloc[:,1]).statistic;assert abs(rr-row.effect)<1e-12
   check.append(dict(sample=sample,endpoint=row.exposure+'__'+row.target,scale=row.scale,rho=rr,status='PASS'))
pd.DataFrame(check).to_csv(P/'tables/INDEPENDENT_CASE_RHO_CHECKS.tsv',sep='\t',index=False)
pd.DataFrame(scorechecks).to_csv(P/'tables/INDEPENDENT_PROGRAM_SCORE_CHECKS.tsv',sep='\t',index=False)
for item in json.loads((P/'INPUT_HASHES.json').read_text()):assert sha(item['path'])==item['sha256']
note='固定4x4主尺度中，代谢转录背景—SPP1 RNA为3/3病例正向（中位rho0.5588，范围0.1618–0.7688）；去SPP1状态程序对应为2正1负，不能称全部病例重复。TYMP在两例的块中位数变化不足，TYMP相关主尺度端点仅1例可估计；该例TYMP—SPP1 rho0.7491、TYMP—代谢背景rho0.4253。不可估计不当作阴性，也不降低区域汇总资格以追求更多阳性。'
f=P/'REPORT_CN.md';old=f.read_text();f.write_text(old.replace('\n## Methods','\n\n'+note+'\n\n## Methods'))
(P/'INDEPENDENT_VALIDATION.json').write_text(json.dumps(dict(status='PASS',at=now(),independent_program_score_checks=len(scorechecks),independent_case_rho_checks=len(check),original_input_hashes=12,no_new_P=True),indent=2))
print('VALIDATED',len(scorechecks),len(check))

