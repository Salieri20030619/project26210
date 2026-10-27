from pathlib import Path
import os,sys,json,datetime,hashlib,subprocess,traceback,gc,ast
import numpy as np,pandas as pd
from scipy import sparse,stats,io
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import dijkstra
P=Path(__file__).resolve().parents[1];ROOT=P.parent;EXT=ROOT/'phase247_external_spatial_cases_20260924_v1';OLD=ROOT/'phase234_part4_replan_execution_20260921_v1'
def now():return datetime.datetime.now().astimezone().isoformat()
def dump(f,o):Path(f).write_text(json.dumps(o,ensure_ascii=False,indent=2,default=lambda x:x.item() if hasattr(x,'item') else str(x)))
def tab(d,f):d.to_csv(P/'tables'/f,sep='\t',index=False)
def sha(f):return hashlib.sha256(Path(f).read_bytes()).hexdigest()
def gate(stage):
 r=subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True)
 with open(P/'logs/memory_gate.log','a') as h:h.write(stage+'\n'+r.stdout+r.stderr)
 if r.returncode:raise RuntimeError('memory '+stage)
tree=ast.parse((OLD/'scripts/common.py').read_text())
exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='fit'],type_ignores=[]),'frozen_fit','exec'),globals())
def main():
 gate('freeze')
 prog=json.loads((OLD/'PROGRAM_DEFINITIONS.json').read_text())['programs']
 names=['polyamine_mta_context','spp1_state_no_spp1','immunosuppression_core','invasion_emt']
 programs={k:prog[k] for k in names}
 ep=[('TYMP_RNA','polyamine_mta_context','mainline'),('TYMP_RNA','SPP1_RNA','mainline'),('polyamine_mta_context','SPP1_RNA','mainline'),('polyamine_mta_context','spp1_state_no_spp1','mainline'),('TYMP_RNA','spp1_state_no_spp1','mainline'),('spp1_state_no_spp1','immunosuppression_core','functional_supplement'),('spp1_state_no_spp1','invasion_emt','functional_supplement')]
 plan=dict(at=now(),priority='independent spatial case controls of TYMP-metabolic context-SPP1 mainline; CD8 only inherited small supplement',external_cases=['sc06T','sc24T','sc27T'],programs=programs,endpoints=ep,primary_scale='fixed4x4 original coordinate-rank blocks, median within block, >=10spots/block and>=8blocks',parallel_scales=['whole_tissue','fixed400um_graph_neighborhood_excludes_self'],score='same Phase234 within-section mean gene Z of log1p CP10k; require >=half genes and>=3 genes,>=20 spots with>=2 genes detected',covariate_sensitivity='external whole-tissue only: log1pUMI, macrophage-marker score, epithelial/malignant-context score; technical adjustment not tumor-cell-specific attribution',old_cohort_comparison='same RNA-level endpoints calculated descriptively for all38 prior sections; do not replace old q05 or current C4 endpoints',statistics='effects/directions only for independent3; no P/FDR or new spot inference; prior cohorts kept separate',no_new_figures=True,known_previous_results=True)
 dump(P/'PLAN_BEFORE_ANALYSIS.json',plan)
 qcs=[];rows=[];spotout=[];inputhash=[]
 def estimates(d,sid,dataset,patient,D=None,cov=None):
  for exposure,target,role in ep:
   for scale in ['whole_tissue','block4x4','neighborhood400um']:
    x=d[exposure].to_numpy();y=d[target].to_numpy();n=len(d);reason=''
    if scale=='block4x4':
     group=d.groupby('block4')[[exposure,target]].median();counts=d.groupby('block4').size();group=group[counts>=10];x=group[exposure].to_numpy();y=group[target].to_numpy()
     rho,n,reason=fit(x,y,min_n=8)
    elif scale=='neighborhood400um':
     nns=np.asarray(D.sum(1)).ravel();W=sparse.diags(1/np.maximum(nns,1))@D
     rho,n,reason=fit(x,W@y,mask=nns>=3,min_n=50)
    else:rho,n,reason=fit(x,y,min_n=50)
    rows.append(dict(sample=sid,dataset=dataset,patient=patient,role=role,exposure=exposure,target=target,scale=scale,model='A0',effect=rho,eligible_spots_or_blocks=n,reason=reason))
   if cov is not None:
    rho,n,reason=fit(d[exposure],d[target],cov=cov,min_n=50)
    rows.append(dict(sample=sid,dataset=dataset,patient=patient,role=role,exposure=exposure,target=target,scale='whole_tissue',model='depth_lineage_adjusted',effect=rho,eligible_spots_or_blocks=n,reason=reason))
 externalqc=pd.read_csv(EXT/'tables/EXTERNAL_INPUT_QC.tsv',sep='\t')
 for r in externalqc.itertuples():
  gate('external_'+r.sample)
  prefix=r.accession+'_'+r.sample+'_ST.'
  files=[EXT/'extracted'/(prefix+s+'.gz') for s in ['matrix.mtx','features.tsv','barcodes.tsv','tissue_positions_list.csv']]
  inputhash.extend(dict(path=str(f),sha256=sha(f)) for f in files)
  mat=io.mmread(files[0]).tocsr().astype(float);g=pd.read_csv(files[1],sep='\t',header=None);b=pd.read_csv(files[2],sep='\t',header=None)[0]
  previous=pd.read_csv(EXT/'tables'/f'{r.sample}_SPOT_MARKERS.tsv.gz',sep='\t')
  ix=pd.Index(b).get_indexer(previous.barcode);assert (ix>=0).all();mat=mat[:,ix].T.tocsr();tot=np.asarray(mat.sum(1)).ravel()
  d=previous.copy();assert len(d)==r.QC_spots
  req=sorted(set(sum(programs.values(),[]))|{'TYMP','SPP1'});raw={}
  for gene in req:
   ids=np.flatnonzero(g[1].eq(gene).to_numpy());raw[gene]=np.asarray(mat[:,ids].sum(1)).ravel() if len(ids) else np.zeros(len(d))
  logs={gene:np.log1p(v/tot*10000) for gene,v in raw.items()}
  for gene in ['TYMP','SPP1']:
   assert np.allclose(d[gene+'_RNA'],logs[gene],atol=1e-12)
  for name,gs in programs.items():
   used=[gene for gene in gs if g[1].eq(gene).any()];arr=np.column_stack([logs[gene] for gene in used]);sd=arr.std(0)
   score=((arr-arr.mean(0))/np.where(sd>0,sd,1)).mean(1)
   detected=(np.column_stack([raw[gene]>0 for gene in used]).sum(1)>=min(2,len(used))).sum()
   good=len(used)>=max(3,int(np.ceil(len(gs)/2))) and detected>=20 and np.std(score)>0
   d[name]=score if good else np.nan
   qcs.append(dict(sample=r.sample,program=name,genes=';'.join(used),missing=';'.join(gene for gene in gs if gene not in used),detected_spots=int(detected),eligible=bool(good)))
  d['block4']=pd.qcut(d.x_um.rank(method='first'),4,labels=False).to_numpy()*4+pd.qcut(d.y_um.rank(method='first'),4,labels=False).to_numpy()
  uv=d[['x_um','y_um']].to_numpy();pairs=cKDTree(uv).query_pairs(100.01,output_type='ndarray')
  G=sparse.coo_matrix((np.full(2*len(pairs),100.),(np.r_[pairs[:,0],pairs[:,1]],np.r_[pairs[:,1],pairs[:,0]])),shape=(len(d),len(d))).tocsr()
  rr=[];cc=[]
  for start in range(0,len(d),256):
   dst=dijkstra(G,directed=False,indices=np.arange(start,min(start+256,len(d))),limit=400);a,c=np.where((dst>0)&np.isfinite(dst));rr.extend(a+start);cc.extend(c)
  D=sparse.coo_matrix((np.ones(len(rr)),(rr,cc)),shape=G.shape).tocsr()
  cov=np.c_[np.log1p(tot),d.Macrophage,d.Malignant_context]
  estimates(d,r.sample,'GSE294006',r.accession,D,cov)
  d.to_csv(P/'inputs'/f'{r.sample}_MAINLINE_INPUTS.tsv.gz',sep='\t',index=False)
  print('EXTERNAL_DONE',r.sample,flush=True);del mat;gc.collect()
 gate('matched_definition_comparison')
 regs=pd.read_csv(OLD/'tables/INPUT_SECTION_REGISTER.tsv',sep='\t')
 for r in regs.itertuples():
  gate('reference_'+r.sample)
  d=pd.read_csv(OLD/'inputs'/f'{r.sample}_INPUTS.tsv.gz',sep='\t')
  d['block4']=d.block_4
  D=sparse.load_npz(OLD/'graphs'/f'{r.sample}_distance800.npz').tocsr();D.data=(D.data<=400).astype(float);D.eliminate_zeros()
  estimates(d,r.sample,r.dataset,r.patient,D)
 table=pd.DataFrame(rows);tab(table,'ALL_MATCHED_DEFINITION_SECTION_EFFECTS.tsv')
 ext=table[table.dataset=='GSE294006'];tab(ext,'EXTERNAL_CASE_MAINLINE_EFFECTS.tsv')
 keys=['dataset','role','exposure','target','scale','model'];patients=[]
 for key,g in table.groupby(keys+['patient'],sort=True):
  v=g.effect.dropna().to_numpy();e=float(np.tanh(np.arctanh(np.clip(v,-.999999999,.999999999)).mean())) if len(v) else np.nan
  patients.append(dict(zip(keys+['patient'],key),effect=e,n_sections=len(v)))
 pat=pd.DataFrame(patients);tab(pat,'PATIENT_EFFECTS.tsv');summary=[]
 for key,g in pat.groupby(keys,sort=True):
  v=g.effect.dropna().to_numpy();summary.append(dict(zip(keys,key),n_patients=len(v),positive=int((v>0).sum()),negative=int((v<0).sum()),zero=int((v==0).sum()),median_rho=float(np.median(v)) if len(v) else np.nan,min_rho=float(v.min()) if len(v) else np.nan,max_rho=float(v.max()) if len(v) else np.nan))
 sm=pd.DataFrame(summary);tab(sm,'COHORT_DIRECTION_COMPARISON.tsv');tab(pd.DataFrame(qcs),'EXTERNAL_PROGRAM_ELIGIBILITY.tsv')
 dump(P/'INPUT_HASHES.json',inputhash)
 for f in inputhash:assert sha(f['path'])==f['sha256']
 assert len(ext)==84 and len(table)==882 and not any(k in table.columns for k in ['P','p','FDR'])
 # Export CD8 only as previously computed companion, no additional thresholds or tests.
 pd.read_csv(EXT/'tables/EXTERNAL_CASE_EFFECTS.tsv',sep='\t').query("endpoint.str.contains('CD8')",engine='python').to_csv(P/'tables/CD8_INHERITED_SMALL_SUPPLEMENT.tsv',sep='\t',index=False)
 lines=['# 独立空间病例对照：TYMP—代谢转录背景—SPP1主线','', '## Results','本轮以新增GSE294006三登记病例为独立来源对照，旧38片仅作同定义方向比较；细胞/spot/区域不扩增患者数。','']
 for r in sm[(sm.dataset=='GSE294006')&(sm.scale=='block4x4')&(sm.role=='mainline')].itertuples():
  lines.append(f'- {r.exposure}—{r.target}: {r.n_patients}例，{r.positive}正/{r.negative}负，中位rho={r.median_rho:.6f}，范围[{r.min_rho:.6f},{r.max_rho:.6f}]。')
 lines+=['','## Methods','沿用三病例官方QC和原始counts，以Phase234固定程序基因和同一log1p每万计数/片内平均Z算法计算。主病例对照为固定4x4块中位数相关，另保留全组织和固定400um图邻域；独立病例全组织再报告深度/谱系程序调整。病例数仅3，不报告推断P/FDR。','主线固定五个两端关联，功能程序两项补充；全部三例、全部尺度和调整结果均输出。旧GLMF/GSE175540/TJ-RCC按相同RNA端点计算方向比较，来源分别汇总，同一患者不同制备Fisher-z合并，不重复计数。','','## Discussion','这是混合spot的RNA区域对应。代谢背景不等同实测MTA，SPP1 RNA或去SPP1程序不等同C4或绝对巨噬丰度；TYMP也未被限定为恶性细胞来源。三例均治疗后，不能据此检验治疗效应。对照变量与旧q05端点不同，不把同源RNA代理当成同一细胞丰度端点的独立复现。','CD8仅引用Phase247已有结果作少量补充，没有新增CD8检验或为改善其结果扩大样本筛选。','', '来源：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294006']
 (P/'REPORT_CN.md').write_text('\n'.join(lines))
 dump(P/'FINAL_VALIDATION.json',dict(status='PASS_DESCRIPTIVE_MAINLINE_CASES',at=now(),external_cases=3,external_effect_records=len(ext),old_sections_compared=38,total_effect_records=len(table),original_input_hashes_unchanged=len(inputhash),P_tests=0,new_figures=0))
 dump(P/'STATUS.json',dict(status='COMPLETE_VERIFIED',at=now()));(P/'EXIT_STATUS').write_text('0\n');(P/'RUN_FINISHED').write_text(now())
 print(sm[(sm.dataset=='GSE294006')&(sm.scale=='block4x4')].to_string(index=False),flush=True)
try:main()
except BaseException:
 (P/'logs/ERROR.txt').write_text(traceback.format_exc());dump(P/'STATUS.json',dict(status='FAILED',error=traceback.format_exc(),at=now()));(P/'EXIT_STATUS').write_text('1\n');raise

