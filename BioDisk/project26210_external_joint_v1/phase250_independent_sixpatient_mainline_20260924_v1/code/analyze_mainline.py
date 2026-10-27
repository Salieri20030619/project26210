from pathlib import Path
import json,datetime,hashlib,subprocess,traceback,ast,gc,sys
import numpy as np,pandas as pd,h5py
from scipy import sparse,stats
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import dijkstra,connected_components
from PIL import Image
P=Path(__file__).resolve().parents[1];ROOT=P.parent;SRC=ROOT/'phase249_independent_spatial_sources_20260924_v1';OLD=ROOT/'phase234_part4_replan_execution_20260921_v1'
def now():return datetime.datetime.now().astimezone().isoformat()
def dump(f,o):Path(f).write_text(json.dumps(o,ensure_ascii=False,indent=2,default=lambda x:x.item() if hasattr(x,'item') else str(x)))
def tab(d,f):d.to_csv(P/'tables'/f,sep='\t',index=False)
def sha(f):return hashlib.sha256(Path(f).read_bytes()).hexdigest()
def gate(stage):
 r=subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True)
 with open(P/'logs/memory_gate.log','a') as f:f.write(stage+'\n'+r.stdout+r.stderr)
 if r.returncode:raise RuntimeError('memory '+stage)
tree=ast.parse((OLD/'scripts/common.py').read_text())
exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='fit'],type_ignores=[]),'frozen_fit','exec'),globals())
def main():
 gate('input_identity_freeze')
 assert (SRC/'DOWNLOAD_EXIT_STATUS').read_text().strip()=='0'
 manifest=json.loads((SRC/'DOWNLOADED_MANIFEST.json').read_text());assert len(manifest)==38 and all(r['status']=='PASS' for r in manifest)
 for r in manifest:
  f=SRC/'downloads'/r['name'];assert f.stat().st_size==r['bytes'] and sha(f)==r['sha256']
 meta=pd.read_csv(SRC/'sources/sdrf_fire.txt',sep='\t').drop_duplicates('Source Name')
 registry=pd.DataFrame(dict(sample=meta['Source Name'],patient='EMTAB12767:'+meta['Characteristics[individual]'].str.replace(' ','',regex=False),disease=meta['Characteristics[disease]'],clinical_history=meta['Characteristics[clinical information]']))
 registry['tissue']=np.where(registry.disease.str.contains('adjacent',case=False),'adjacent','tumor')
 assert len(registry)==9 and registry.patient.nunique()==6 and (registry.tissue=='tumor').sum()==7
 registry['treatment_context']=np.where(registry.clinical_history.str.contains('nivolumab|ipilimumab',case=False),'ICB_EXPOSURE_RECORDED','NOT_EXPLICITLY_RECORDED')
 tab(registry,'CASE_IDENTITY_FROZEN.tsv')
 inherited=json.loads((ROOT/'phase248_independent_cases_mainline_20260924_v1/PLAN_BEFORE_ANALYSIS.json').read_text())
 programs=inherited['programs'];endpoints=inherited['endpoints']
 marker={'Macrophage':['CD68','CSF1R','FCER1G','AIF1'],'Malignant_context':['EPCAM','KRT8','KRT18','CA9'],'CD8lineage':['CD3D','CD3E','TRAC','CD8A','CD8B']}
 plan=dict(at=now(),source='E-MTAB-12767 public official processed data',patients=6,capture_arrays=9,tumor_arrays=7,adjacent_arrays=2,sampling='all registered patients; HG1 multiple tumor preparations pooled within patient; LG2 capture includes primary TME and thrombus, not split without annotation',QC='in_tissue, total UMI>=500, genes_detected>500, mitochondrial fraction<0.30, corresponding to publication spot criteria; no outcome selection',programs=programs,markers=marker,mainline_endpoints=endpoints,primary='fixed4x4 block medians,>=10spots/block,>=8blocks; same per-gene Z log1pCP10k scoring as Phase248',parallel='whole-tissue and fixed400um graph; whole-tissue depth/macrophage/malignant-context adjusted sensitivity',paired_control='HG1 andHG3 adjacent versus own tumor, mean log1pCP10k signals, no comparison of within-section Z means; two donors, descriptive',CD8='single400um RNA-program endpoint A0 plus depth/lineage/source-target adjusted; optional descriptive only, no distances or expanded pathway search',statistical_unit='independent author-registered individual; all outputs descriptive because six patients; no P/FDR and no pooled9 inference',no_figures=True,no_new_mapping=True,prior_mainline_three_case_results_known=True)
 dump(P/'PLAN_BEFORE_ANALYSIS.json',plan)
 rows=[];means=[];quality=[];geneqc=[]
 for r in registry.itertuples():
  gate('sample_'+r.sample);dump(P/'STATUS.json',dict(status='CASE_MAINLINE_ANALYSIS',at=now(),sample=r.sample))
  f=SRC/'downloads'/f'{r.sample}_filtered_feature_bc_matrix.h5'
  with h5py.File(f,'r') as h:
   a=h['matrix'];m=sparse.csc_matrix((a['data'][:],a['indices'][:],a['indptr'][:]),shape=tuple(a['shape'][:])).tocsr()
   genes=np.array([x.decode() for x in a['features/name'][:]]);barcodes=np.array([x.decode() for x in a['barcodes'][:]])
  assert m.shape==(len(genes),len(barcodes)) and len(set(barcodes))==len(barcodes) and np.isfinite(m.data).all() and (m.data>=0).all() and np.allclose(m.data,np.round(m.data))
  pos=pd.read_csv(SRC/'downloads'/f'{r.sample}_tissue_positions_list.csv',header=None,names=['barcode','in_tissue','array_row','array_col','px_row','px_col'])
  if pos.barcode.iloc[0]=='barcode':pos=pos.iloc[1:]
  pos=pos.set_index('barcode').loc[barcodes].apply(pd.to_numeric)
  total=np.asarray(m.sum(0)).ravel();ng=np.asarray((m>0).sum(0)).ravel();mt=np.asarray(m[np.char.startswith(genes,'MT-')].sum(0)).ravel()/np.maximum(total,1)
  keep=(pos.in_tissue.to_numpy()==1)&(total>=500)&(ng>500)&(mt<.30)
  original=len(barcodes);m=m[:,keep].T.tocsr().astype(float);total=total[keep];pos=pos.iloc[np.flatnonzero(keep)];barcodes=barcodes[keep];n=len(barcodes);assert n>=100
  req=sorted(set(sum(programs.values(),[])+sum(marker.values(),[])+['TYMP','SPP1']))
  raw={g:np.asarray(m[:,np.flatnonzero(genes==g)].sum(1)).ravel() if (genes==g).any() else np.zeros(n) for g in req}
  log={g:np.log1p(v/total*10000) for g,v in raw.items()}
  uv=np.c_[pos.array_col.to_numpy()*50,pos.array_row.to_numpy()*np.sqrt(3)*50]
  d=pd.DataFrame(dict(barcode=barcodes,x_um=uv[:,0],y_um=uv[:,1],TYMP_RNA=log['TYMP'],SPP1_RNA=log['SPP1'],log1p_UMI=np.log1p(total)))
  detect={}
  for name,gs in {**programs,**marker}.items():
   used=[g for g in gs if (genes==g).any()]
   vv=np.column_stack([log[g] for g in used]);sd=vv.std(0);score=((vv-vv.mean(0))/np.where(sd>0,sd,1)).mean(1)
   detected=np.column_stack([raw[g]>0 for g in used]).sum(1)>=min(2,len(used))
   if name=='CD8lineage':detected &= (raw['CD8B']>0)&np.column_stack([raw[g]>0 for g in ['CD3D','CD3E','TRAC']]).any(1)
   good=len(used)>=max(3,int(np.ceil(len(gs)/2))) and detected.sum()>=20 and score.std()>0
   if name in marker:good &= len(used)==len(gs)
   d[name]=score if good else np.nan;detect[name]=detected
   geneqc.append(dict(sample=r.sample,program=name,used=';'.join(used),missing=';'.join(g for g in gs if g not in used),eligible=bool(good),detected_spots=int(detected.sum())))
   if name in programs:means.append(dict(sample=r.sample,patient=r.patient,tissue=r.tissue,variable=name,mean_logCP10k=float(vv.mean()) if len(used)==len(gs) else np.nan,full_genes=len(used)==len(gs)))
  for g in ['TYMP','SPP1']:means.append(dict(sample=r.sample,patient=r.patient,tissue=r.tissue,variable=g+'_RNA',mean_logCP10k=float(log[g].mean()),full_genes=bool((genes==g).any())))
  d['block4']=pd.qcut(d.x_um.rank(method='first'),4,labels=False).to_numpy()*4+pd.qcut(d.y_um.rank(method='first'),4,labels=False).to_numpy()
  pairs=cKDTree(uv).query_pairs(100.01,output_type='ndarray');assert np.allclose(np.linalg.norm(uv[pairs[:,0]]-uv[pairs[:,1]],axis=1),100)
  G=sparse.coo_matrix((np.full(len(pairs)*2,100.),(np.r_[pairs[:,0],pairs[:,1]],np.r_[pairs[:,1],pairs[:,0]])),shape=(n,n)).tocsr()
  rr=[];cc=[]
  for start in range(0,n,256):
   dist=dijkstra(G,directed=False,indices=np.arange(start,min(start+256,n)),limit=400);a,b=np.where((dist>0)&np.isfinite(dist));rr.extend(a+start);cc.extend(b)
  H=sparse.coo_matrix((np.ones(len(rr)),(rr,cc)),shape=G.shape).tocsr();nn=np.asarray(H.sum(1)).ravel();W=sparse.diags(1/np.maximum(nn,1))@H
  cov=np.c_[d.log1p_UMI,d.Macrophage,d.Malignant_context]
  def add(x,y,role,scale,model,e,ne,why):rows.append(dict(sample=r.sample,patient=r.patient,tissue=r.tissue,role=role,exposure=x,target=y,scale=scale,model=model,effect=e,eligible_spots_or_blocks=ne,reason=why))
  for x,y,role in endpoints:
   for scale in ['whole_tissue','block4x4','neighborhood400um']:
    if scale=='whole_tissue':e,ne,why=fit(d[x],d[y])
    elif scale=='block4x4':
     z=d.groupby('block4')[[x,y]].median();cnt=d.groupby('block4').size();z=z[cnt>=10];e,ne,why=fit(z[x],z[y],min_n=8)
    else:e,ne,why=fit(d[x],W@d[y].to_numpy(),mask=nn>=3)
    add(x,y,role,scale,'A0',e,ne,why)
   e,ne,why=fit(d[x],d[y],cov=cov);add(x,y,role,'whole_tissue','depth_lineage_adjusted',e,ne,why)
  mask=detect['Macrophage']&(nn>=3)
  for model in ['A0','A3']:
   c=None if model=='A0' else np.c_[cov,d.CD8lineage,W@d.log1p_UMI.to_numpy(),W@d.Macrophage.to_numpy(),W@d.Malignant_context.to_numpy()]
   e,ne,why=fit(d.SPP1_RNA,W@d.CD8lineage.to_numpy(),cov=c,mask=mask)
   add('SPP1_RNA','CD8lineage','CD8_small_supplement','neighborhood400um',model,e,ne,why)
  sf=json.loads((SRC/'downloads'/f'{r.sample}_scalefactors_json.json').read_text());image=Image.open(SRC/'downloads'/f'{r.sample}_tissue_hires_image.png')
  xy=pos[['px_col','px_row']].to_numpy()*sf['tissue_hires_scalef'];inbounds=(xy[:,0]>=0)&(xy[:,1]>=0)&(xy[:,0]<image.width)&(xy[:,1]<image.height);assert inbounds.all()
  quality.append(dict(sample=r.sample,patient=r.patient,tissue=r.tissue,input_spots=original,QC_spots=n,genes=len(genes),mitochondrial_fraction_max=float(mt[keep].max()),components=int(connected_components(G,directed=False,return_labels=False)),image_coordinate_pass=True,TYMP_positive_spots=int((raw['TYMP']>0).sum()),SPP1_positive_spots=int((raw['SPP1']>0).sum())))
  d.to_csv(P/'inputs'/f'{r.sample}_MAINLINE_INPUTS.tsv.gz',sep='\t',index=False)
  print('CASE_DONE',r.sample,n,flush=True);del m;gc.collect()
 gate('patient_summary')
 effects=pd.DataFrame(rows);tab(effects,'ALL_SECTION_EFFECTS.tsv');tab(pd.DataFrame(quality),'INPUT_QC.tsv');tab(pd.DataFrame(geneqc),'PROGRAM_ELIGIBILITY.tsv')
 keys=['tissue','role','exposure','target','scale','model'];patients=[]
 for key,g in effects.groupby(keys+['patient'],sort=True):
  v=g.effect.dropna().to_numpy();e=float(np.tanh(np.arctanh(np.clip(v,-.999999999,.999999999)).mean())) if len(v) else np.nan
  patients.append(dict(zip(keys+['patient'],key),effect=e,n_sections=len(v)))
 patient=pd.DataFrame(patients);tab(patient,'PATIENT_EFFECTS.tsv');summary=[]
 for key,g in patient.groupby(keys,sort=True):
  v=g.effect.dropna().to_numpy();summary.append(dict(zip(keys,key),n_patients=len(v),positive=int((v>0).sum()),negative=int((v<0).sum()),zero=int((v==0).sum()),median_rho=float(np.median(v)) if len(v) else np.nan,min_rho=float(v.min()) if len(v) else np.nan,max_rho=float(v.max()) if len(v) else np.nan))
 sm=pd.DataFrame(summary);tab(sm,'PATIENT_DIRECTION_SUMMARY.tsv')
 means=pd.DataFrame(means);tab(means,'SECTION_MEAN_LOG_SIGNALS.tsv');means=means.groupby(['patient','tissue','variable']).mean(numeric_only=True).reset_index()
 paired=means.pivot(index=['patient','variable'],columns='tissue',values='mean_logCP10k').dropna().reset_index();paired['tumor_minus_adjacent']=paired.tumor-paired.adjacent
 tab(paired,'PAIRED_TUMOR_ADJACENT_DESCRIPTIVE.tsv');assert paired.patient.nunique()==2
 assert len(effects)==270 and patient[patient.tissue=='tumor'].patient.nunique()==6
 assert not any(c in effects.columns for c in ['P','p','FDR'])
 lines=['# E-MTAB-12767独立六患者主线空间对照','', '## Results','九张捕获阵列来自六个官方登记个体；七张肿瘤阵列按六患者汇总，两张旁组织用于同患者描述性对照。HG1两张肿瘤阵列不重复计算患者数。','']
 for r in sm[(sm.tissue=='tumor')&(sm.role=='mainline')&(sm.scale=='block4x4')].itertuples():
  lines.append(f'- {r.exposure}—{r.target}：有效{r.n_patients}患者，{r.positive}正/{r.negative}负，中位rho={r.median_rho:.6f}，范围[{r.min_rho:.6f},{r.max_rho:.6f}]。')
 lines+=['','## Methods','官方38文件按公开大小及本地SHA256完整核验。采用文中spot质量门槛：in_tissue、UMI≥500、检出基因>500、线粒体比例<30%。固定五条主线与两条功能补充；主尺度为4x4区域中位数，全部病例同时保留全组织和400um图邻域。基因集和评分沿用Phase248，整套定义在计算前保存。','同患者肿瘤切片效应Fisher-z等权合并。旁组织比较使用每万计数log1p信号均值，不比较片内Z分数均值；两配对供者仅描述。六患者均按方向/效应呈现，不作P/FDR或治疗组比较。CD8仅两种固定400um口径作小型补充。','','## Discussion','来源独立依据官方登记、机构与研究样本表，未做跨研究基因型鉴定。LG2同一阵列包含原发肿瘤与肿瘤栓，缺少独立区域标注时不人为拆成两患者；HG1有ICB暴露记录，其余未明确记录不等同确证未治疗。','混合spot的TYMP、SPP1及代谢转录程序不能替代细胞特异来源、实测MTA、C4丰度或因果中介。旁组织也不能自动作为完全健康肾组织。保留阴性、尺度差异、调整减弱及缺失资格。','', '来源：https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-12767 ；https://www.nature.com/articles/s41698-023-00441-5']
 (P/'REPORT_CN.md').write_text('\n'.join(lines))
 dump(P/'FINAL_VALIDATION.json',dict(status='PASS_DESCRIPTIVE_INDEPENDENT_CASES',at=now(),official_files_checked=38,capture_arrays=9,patients=6,tumor_arrays=7,adjacent_paired_patients=2,section_effects=len(effects),QC_spots=int(pd.DataFrame(quality).QC_spots.sum()),new_figures=0,new_models_trained=0,P_tests=0))
 dump(P/'STATUS.json',dict(status='COMPLETE_VERIFIED',at=now()));(P/'EXIT_STATUS').write_text('0\n');(P/'RUN_FINISHED').write_text(now())
 print(sm[(sm.tissue=='tumor')&(sm.role=='mainline')&(sm.scale=='block4x4')].to_string(index=False),flush=True)
try:main()
except BaseException:
 (P/'logs/ERROR.txt').write_text(traceback.format_exc());dump(P/'STATUS.json',dict(status='FAILED',at=now(),error=traceback.format_exc()));(P/'EXIT_STATUS').write_text('1\n');raise

