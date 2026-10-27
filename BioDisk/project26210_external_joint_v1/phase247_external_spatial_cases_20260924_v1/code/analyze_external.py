from pathlib import Path
import json,datetime,subprocess,traceback,hashlib,gzip,ast,gc,sys
import numpy as np,pandas as pd,h5py
from scipy import sparse,stats,io
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import dijkstra,connected_components
from PIL import Image
P=Path(__file__).resolve().parents[1];ROOT=P.parent
def now():return datetime.datetime.now().astimezone().isoformat()
def dump(f,o):Path(f).write_text(json.dumps(o,ensure_ascii=False,indent=2,default=lambda x:x.item() if hasattr(x,'item') else str(x)))
def tab(d,f):d.to_csv(P/'tables'/f,sep='\t',index=False)
def sha(f):return hashlib.sha256(Path(f).read_bytes()).hexdigest()
def gate(stage):
 r=subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True)
 with open(P/'logs/memory_gate.log','a') as f:f.write(stage+'\n'+r.stdout+r.stderr)
 if r.returncode:raise RuntimeError('memory '+stage)
tree=ast.parse((ROOT/'phase234_part4_replan_execution_20260921_v1/scripts/common.py').read_text())
exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='fit'],type_ignores=[]),'frozen_fit','exec'),globals())
def main():
 gate('preflight')
 assert (P/'DOWNLOAD_EXIT_STATUS').read_text().strip()=='0'
 dl=json.loads((P/'DOWNLOAD_VALIDATION.json').read_text());assert dl['status']=='PASS' and dl['archive_bytes']==111564800
 files=pd.read_csv(P/'sources/filelist.txt',sep='\t');check=files[files.iloc[:,0]=='File']
 manifest=json.loads((P/'EXTRACTED_MANIFEST.json').read_text());assert len(check)==len(manifest)==27
 expected=dict(zip(check.Name,check.Size.astype(int)))
 for item in manifest:
  f=P/'extracted'/item['member'];assert f.stat().st_size==expected[f.name] and sha(f)==item['sha256']
  with gzip.open(f,'rb') as h:
   while h.read(1024*1024):pass
 cases=[('GSM8896252','sc06T'),('GSM8896253','sc24T'),('GSM8896254','sc27T')]
 markers=json.loads((ROOT/'phase246_spatial_marker_structure_20260924_v1/PLAN_BEFORE_ANALYSIS.json').read_text())['markers']
 markers['Malignant_context']=['EPCAM','KRT8','KRT18','CA9']
 plan=dict(at=now(),source='GSE294006 official processed count and spatial files',cases=cases,markers=markers,QC='in_tissue==1 and total counts>=1000, per author metadata; barcodes one-to-one and original raw counts integer',source_definition='continuous SPP1 log1p CP10k in spots detecting>=2 macrophage marker genes; RNA coexpression proxy, not old V3 or C4',target_definition='same fixed lineage sets and co-detection as Phase246',graph='Visium array horizontal pitch100um, first ring 100um, shortest path400um excluding source; no bridging holes/components',models='A0; A3 source+neighbor log1p total counts, macrophage marker score, malignant-context score, and source target marker',distance='high/low20% source SPP1 RNA, same marker top20% target eligibility and >=20 each reachable group as Phase246',biological_inference='three registered cases, effects and directions only; no spot-level or patient-level P/FDR; no treatment group comparison',new_figures=0,original_inputs='read only',known_prior_results=True)
 dump(P/'PLAN_BEFORE_ANALYSIS.json',plan)
 results=[];qcs=[];geneqc=[];identities=[]
 for accession,sample in cases:
  gate(sample);dump(P/'STATUS.json',dict(status='CASE_ANALYSIS',at=now(),sample=sample))
  prefix=accession+'_'+sample+'_ST.'
  def file(suffix):return P/'extracted'/(prefix+suffix+'.gz')
  m=io.mmread(file('matrix.mtx')).tocsr().astype(float)
  genes=pd.read_csv(file('features.tsv'),sep='\t',header=None,dtype=str)
  bc=pd.read_csv(file('barcodes.tsv'),sep='\t',header=None)[0].astype(str)
  pos=pd.read_csv(file('tissue_positions_list.csv'),header=None,names=['barcode','in_tissue','array_row','array_col','px_row','px_col'])
  if pos.iloc[0].barcode=='barcode':pos=pos.iloc[1:].copy()
  pos=pos.set_index('barcode');assert pos.index.is_unique and bc.is_unique and m.shape==(len(genes),len(bc))
  pos=pos.loc[bc].apply(pd.to_numeric)
  assert np.isfinite(m.data).all() and (m.data>=0).all() and np.allclose(m.data,np.round(m.data))
  tot=np.asarray(m.sum(0)).ravel();keep=(tot>=1000)&(pos.in_tissue.to_numpy()==1)
  m=m[:,keep].T.tocsr();tot=tot[keep];pos=pos.loc[bc[keep]].copy();b=bc[keep].to_numpy()
  assert len(b)>=100
  req=sorted(set(sum(markers.values(),[]))|{'SPP1','TYMP'});raw={}
  for g in req:
   idx=np.flatnonzero(genes.iloc[:,1].to_numpy()==g)
   raw[g]=np.asarray(m[:,idx].sum(1)).ravel() if len(idx) else np.zeros(len(b))
  logged={g:np.log1p(v/tot*10000) for g,v in raw.items()}
  score={};det={};eligible={}
  for name,gs in markers.items():
   missing=[g for g in gs if g not in set(genes.iloc[:,1])]
   a=np.column_stack([logged[g] for g in gs]);sd=a.std(0);s=((a-a.mean(0))/np.where(sd>0,sd,1)).mean(1)
   detected=np.column_stack([raw[g]>0 for g in gs]).sum(1)>=2
   if name=='CD8lineage':detected &= (raw['CD8B']>0)&np.column_stack([raw[g]>0 for g in markers['Tcell']]).any(1)
   ok=not missing and detected.sum()>=20 and s.std()>0
   score[name]=s if ok else np.full(len(b),np.nan);det[name]=detected;eligible[name]=ok
   geneqc.append(dict(accession=accession,sample=sample,program=name,eligible=ok,detected_spots=int(detected.sum()),missing=';'.join(missing)))
  uv=np.c_[pos.array_col.to_numpy()*50,pos.array_row.to_numpy()*np.sqrt(3)*50]
  pairs=cKDTree(uv).query_pairs(100.01,output_type='ndarray')
  assert len(pairs)>0
  lengths=np.linalg.norm(uv[pairs[:,0]]-uv[pairs[:,1]],axis=1);assert np.allclose(lengths,100)
  G=sparse.coo_matrix((np.full(len(pairs)*2,100.),(np.r_[pairs[:,0],pairs[:,1]],np.r_[pairs[:,1],pairs[:,0]])),shape=(len(b),len(b))).tocsr()
  comps,labels=connected_components(G,directed=False)
  rr=[];cc=[]
  for start in range(0,len(b),256):
   distances=dijkstra(G,directed=False,indices=np.arange(start,min(start+256,len(b))),limit=400)
   ii,jj=np.where((distances>0)&np.isfinite(distances));rr.extend(ii+start);cc.extend(jj)
  H=sparse.coo_matrix((np.ones(len(rr)),(rr,cc)),shape=G.shape).tocsr();nn=np.asarray(H.sum(1)).ravel();W=sparse.diags(1/np.maximum(nn,1))@H
  src=det['Macrophage']&(nn>=3)&eligible['Macrophage'];x=logged['SPP1'];depth=np.log1p(tot)
  for target in ['CD8lineage','Endothelial','Fibroblast']:
   y=score[target]
   for model in ['A0','A3']:
    cov=None if model=='A0' else np.column_stack([depth,W@depth,score['Macrophage'],W@score['Macrophage'],score['Malignant_context'],W@score['Malignant_context'],y])
    effect,n,reason=fit(x,W@y,cov=cov,mask=src)
    results.append(dict(accession=accession,sample=sample,endpoint='SPP1_RNA__'+target+'_400um',model=model,effect_type='rho',effect=effect,eligible_spots=n,reason=reason,n_registered_cases_for_effect=1))
   tm=det[target]&np.isfinite(y);mask=det['Macrophage']&eligible['Macrophage'];delta=np.nan;reason='';nt=nh=nl=0
   if mask.sum()<100 or tm.sum()<20:reason='insufficient_qualified_source_or_target'
   else:
    hi,lo=np.quantile(x[mask],[.8,.2]);ids=np.flatnonzero(tm&(y>=np.quantile(y[tm],.8)));nt=len(ids)
    if nt<20 or hi<=lo:reason='insufficient_rich_targets_or_source_variation'
    else:
     ds=dijkstra(G,directed=False,indices=ids,min_only=True)
     highs=mask&(x>=hi)&np.isfinite(ds);lows=mask&(x<=lo)&np.isfinite(ds);nh=int(highs.sum());nl=int(lows.sum())
     if min(nh,nl)<20:reason='insufficient_reachable_source_groups'
     else:delta=float(np.median(ds[highs])-np.median(ds[lows]))
   results.append(dict(accession=accession,sample=sample,endpoint='SPP1_high_minus_low_distance_to_'+target,model='fixed20percent',effect_type='difference_um',effect=delta,eligible_spots=nh+nl,reason=reason,target_spots=nt,high_sources=nh,low_sources=nl,n_registered_cases_for_effect=1))
  sf=json.loads(gzip.decompress(file('scalefactors_json.json').read_bytes()))
  from io import BytesIO
  image=Image.open(BytesIO(gzip.decompress(file('tissue_hires_image.png').read_bytes())))
  xy=pos[['px_col','px_row']].to_numpy()*sf['tissue_hires_scalef'];bounds=(xy[:,0]>=0)&(xy[:,1]>=0)&(xy[:,0]<image.width)&(xy[:,1]<image.height);assert bounds.all()
  qcs.append(dict(accession=accession,sample=sample,input_genes=len(genes),input_spots=len(bc),QC_spots=len(b),UMI_min=float(tot.min()),components=int(comps),first_ring_edges=len(pairs),median_edge_um=float(np.median(lengths)),image_coordinates_in_bounds=int(bounds.sum()),raw_counts_integer=True,macro_marker_source_spots=int(src.sum())))
  txt=(P/'sources'/f'{accession}.soft.txt').read_text()
  attrs=[s.split(' = ',1)[1] for s in txt.splitlines() if s.startswith('!Sample_characteristics_ch1')]
  identities.append(dict(accession=accession,sample=sample,attributes='; '.join(attrs),identity='distinct registered sample aliases; no genotypic patient identity verification',no_treatment_comparison=True))
  pd.DataFrame(dict(barcode=b,x_um=uv[:,0],y_um=uv[:,1],SPP1_RNA=x,TYMP_RNA=logged['TYMP'],**score)).to_csv(P/'tables'/f'{sample}_SPOT_MARKERS.tsv.gz',sep='\t',index=False)
  print('CASE_COMPLETE',sample,len(b),flush=True)
  del m;gc.collect()
 table=pd.DataFrame(results);tab(table,'EXTERNAL_CASE_EFFECTS.tsv');tab(pd.DataFrame(qcs),'EXTERNAL_INPUT_QC.tsv');tab(pd.DataFrame(geneqc),'EXTERNAL_MARKER_ELIGIBILITY.tsv');tab(pd.DataFrame(identities),'OFFICIAL_CASE_METADATA.tsv')
 assert len(table)==27 and table[['sample','endpoint','model']].duplicated().sum()==0
 assert not any(c in table.columns for c in ['P','p','FDR'])
 lines=['# GSE294006三例治疗后肾癌空间表达对照','', '## Results','三个官方登记病例的处理计数、空间坐标和图像通过输入核验。按已冻结定义计算表达邻域和区域距离，全部效应与不可估计原因见表；仅报告病例方向，不进行人群或治疗组检验。','']
 for r in results:
  if r['model']=='A0':lines.append(f"- {r['sample']} {r['endpoint']}: rho={r['effect']}, 合格来源spot={r['eligible_spots']}, 原因={r['reason']}。")
 lines+=['','## Methods','官方GSE294006处理归档111564800字节、27成员按官方filelist逐个大小核对，gzip完整解码与归档SHA256登记。保留in_tissue且UMI≥1000的spot，沿用作者QC阈值；原始计数按每万计数log1p后形成冻结谱系标记分数。图使用Visium阵列100um第一圈及400um最短路径邻域。','来源采用检测到至少两个巨噬标记的spot上的SPP1连续RNA，目标和距离资格沿用Phase246。输出A0及深度/组成程序/来源目标调整的A3效果。所有病例和端点均报告，没有按相关方向筛选。','','## Discussion','这是新增来源的三例表达空间对照，来源变量是RNA共表达代理，与原SPP1巨噬内比例不同，不称同一端点的直接独立复现或C4识别成功。样本为治疗后组织，官方病例包含不同治疗及病理反应；不能据三例比较药物效果，也不能以spot代替独立患者数。','Visium混合spot与H&E不提供独立CD8/CAF细胞分割或蛋白身份；本轮不能补成已确认的物理接触/排斥。所有方向与资格不足结果均保留。','', '来源：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294006']
 (P/'REPORT_CN.md').write_text('\n'.join(lines))
 dump(P/'FINAL_VALIDATION.json',dict(status='PASS_DESCRIPTIVE_EXTERNAL_CASES',at=now(),official_file_sizes_checked=27,gzip_integrity_checked=27,cases=3,section_effect_records=len(table),biological_P_tests=0,new_figures=0,all_image_coordinates_checked=True))
 dump(P/'STATUS.json',dict(status='COMPLETE_VERIFIED',at=now()));(P/'EXIT_STATUS').write_text('0\n');(P/'RUN_FINISHED').write_text(now())
 print(table.to_string(index=False),flush=True)
try:main()
except BaseException:
 (P/'logs/ERROR_analysis.txt').write_text(traceback.format_exc());dump(P/'STATUS.json',dict(status='FAILED',error=traceback.format_exc(),at=now()));(P/'EXIT_STATUS').write_text('1\n');raise

