from common import *
import h5py
from anndata.io import read_elem
from scipy import sparse,stats
from scipy.optimize import nnls
from scipy.spatial import cKDTree
CLASSES=['Fibroblast_like','Mural','Endothelial','Epithelial','Macrophage','Tcell']
def train(pb,v,features,held=None):
 mask=np.ones(len(pb),bool) if held is None else pb.study.ne(held).to_numpy()
 f=pb[mask&pb.label.eq('Fibroblast_like')]
 if held is not None and (len(f)<2 or f.study.nunique()<2):return None
 cols=[]
 for cl in CLASSES:
  z=mask&pb.label.eq(cl)
  if not z.any():return None
  cols.append(np.mean([v[z&pb.study.eq(st)][:,features].mean(0) for st in pb.loc[z,'study'].unique()],axis=0))
 mat=np.array(cols).T;scale=np.sqrt(np.maximum(mat.mean(1),1.))
 return mat/scale[:,None],scale
def fitcoef(model,y):
 mat,scale=model;coef,res=nnls(mat,y/scale)
 prop=coef/coef.sum() if coef.sum()>0 else coef
 return coef,prop,res
def feature_gate(pb,v,jj):
 pred=[];mix=[]
 for st in pb.study.unique():
  model=train(pb,v,jj,st)
  if model is None:continue
  test=pb.study.eq(st).to_numpy()
  for i in np.flatnonzero(test):
   _,prop,_=fitcoef(model,v[i,jj]);pred.append((pb.label.iloc[i],CLASSES[int(prop.argmax())]))
  bg=[v[test&pb.label.eq(cl)][:,jj].mean(0) for cl in CLASSES[1:] if (test&pb.label.eq(cl)).any()]
  if len(bg)!=5:continue
  bg=np.mean(bg,axis=0)
  for i in np.flatnonzero(test&pb.label.eq('Fibroblast_like')):
   for w in [0.,.1,.25,.5]:
    _,prop,_=fitcoef(model,w*v[i,jj]+(1-w)*bg);mix.append((w,float(prop[0])))
 pos=[x==y for x,y in pred if x=='Fibroblast_like'];neg=[y=='Fibroblast_like' for x,y in pred if x!='Fibroblast_like']
 recall=float(np.mean(pos)) if pos else np.nan;fp=float(np.mean(neg)) if neg else np.nan
 mae=float(np.mean([abs(x-y) for x,y in mix])) if mix else np.nan;null95=float(np.quantile([y for x,y in mix if x==0],.95)) if mix else np.nan
 return dict(recall=recall,false_positive=fp,mixture_MAE=mae,null95=null95,heldout_fibro=len(pos),passed=bool(recall>=.8 and fp<=.1 and mae<=.15 and null95<=.15))
def partial(x,y,covs,mask):
 ok=mask&np.isfinite(x)&np.isfinite(y)&np.isfinite(covs).all(1)
 n=int(ok.sum())
 if n<100:return np.nan,n,'fewer_than_100_spots'
 xx=stats.rankdata(x[ok]);yy=stats.rankdata(y[ok])
 if np.std(xx)<1e-9 or np.std(yy)<1e-9:return np.nan,n,'constant_endpoint'
 cc=np.column_stack([stats.rankdata(covs[ok,j]) for j in range(covs.shape[1])])
 cc=cc[:,np.std(cc,axis=0)>1e-9];cc=(cc-cc.mean(0))/cc.std(0) if cc.shape[1] else cc
 design=np.column_stack([np.ones(n),cc]);rx=xx-design@np.linalg.lstsq(design,xx,rcond=None)[0];ry=yy-design@np.linalg.lstsq(design,yy,rcond=None)[0]
 if np.std(rx)<1e-8 or np.std(ry)<1e-8:return np.nan,n,'constant_residual'
 return float(np.corrcoef(rx,ry)[0,1]),n,'OK'
def main():
 gate('spatial_pilot')
 pb=read(P/'tables/REFERENCE_PSEUDOBULK_REGISTER.tsv')
 vv=read(P/'inputs/REFERENCE_PSEUDOBULK_CPM.tsv.gz');features=vv.columns.tolist();v=vv.to_numpy()
 audit=R/'phase286_CAF_identity_boundary_audit_20260929_v1';a=read(audit/'inputs/CELL_IDENTITY_AUDIT.tsv.gz');z=np.load(audit/'inputs/IDENTITY_PANEL_COUNTS.npz')
 ix=[z['genes'].tolist().index(g) for g in features];counts=z['counts'][:,ix];total=z['total_counts']
 neg=[]
 unseen=['B cell','NK cell','Plasma cell','Monocyte','Dendritic cell','Mast cell','Neutrophil']
 for (st,pat,label),g in a[a.qc_pass&a.cell_type_broad_v1.isin(unseen)].groupby(['study','curated_patient_id','cell_type_broad_v1']):
  if len(g)<20:continue
  model=train(pb,v,np.arange(len(features)),st)
  if model is None:continue
  ii=g.index.to_numpy();y=counts[ii].sum(0)/total[ii].sum()*1e6
  coef,prop,res=fitcoef(model,y)
  neg.append(dict(study=st,patient=pat,original_label=label,cells=len(g),predicted=CLASSES[int(prop.argmax())],fibro_weight=float(prop[0]),residual=res))
 nd=pd.DataFrame(neg);tab(nd,'UNMODELED_LINEAGE_CHALLENGE.tsv')
 q95=float(nd.fibro_weight.quantile(.95));fp=float(nd.predicted.eq('Fibroblast_like').mean());negpass=q95<=.15 and fp<=.1
 dump(P/'qc/UNMODELED_BACKGROUND_GATE.json',dict(status='PASS' if negpass else 'NOT_PASSED',pseudobulks=len(nd),fibro_q95=q95,false_positive=fp))
 print('BACKGROUND_GATE',negpass,q95,fp,flush=True)
 if not negpass:
  dump(P/'qc/SPATIAL_PILOT_COMPLETE.json',dict(status='STOPPED_AT_BACKGROUND_GATE',new_associations=0,reason='unmodeled_lineage_false_positive_risk'));return
 reg=read(R/'phase234_part4_replan_execution_20260921_v1/tables/INPUT_SECTION_REGISTER.tsv')
 quals=[];effects=[];cache={};manifest=[]
 for r in reg.itertuples():
  inp=R/'phase234_part4_replan_execution_20260921_v1/inputs'/f'{r.sample}_INPUTS.tsv.gz';d=read(inp)
  f=Path(r.input_path);st=f.stat();manifest.append(dict(path=str(f),bytes=st.st_size,mtime_ns=st.st_mtime_ns,input_table_sha256=sha(inp)))
  with h5py.File(f,'r') as h:
   gs=np.asarray(read_elem(h['var'][h['var'].attrs['_index']])).astype(str);ids=np.asarray(read_elem(h['obs'][h['obs'].attrs['_index']])).astype(str)
   if 'gene_symbol_original' in h['var']:gs=np.asarray(read_elem(h['var/gene_symbol_original'])).astype(str)
   assert len(set(gs))==len(gs)
   inds=pd.Index(ids).get_indexer(d.spot_id.astype(str));assert (inds>=0).all()
   jj=tuple(i for i,g in enumerate(features) if g in set(gs))
   present=[features[i] for i in jj]
   basic=len(jj)>=.8*len(features) and len(set(present)&{'DCN','LUM','PDGFRA'})>=2 and len(set(present)&{'COL1A1','COL1A2','COL3A1'})>=2
   if jj not in cache:cache[jj]=feature_gate(pb,v,list(jj)) if basic else dict(passed=False,reason='feature_coverage')
   tech=cache[jj];qualified=basic and tech['passed']
   rec=dict(sample=r.sample,dataset=r.dataset,patient=r.patient,present=len(jj),total=len(features),missing=';'.join(g for g in features if g not in present),qualified=qualified,**tech)
   if not qualified:quals.append(rec);print('SPATIAL_SKIP',r.sample,flush=True);continue
   node=h['layers/counts'] if 'layers/counts' in h else h['X']
   raw=read_elem(node).tocsr();assert raw.shape==(len(ids),len(gs)) and (raw.data>=0).all() and np.allclose(raw.data,np.round(raw.data))
   sums=np.asarray(raw.sum(1)).ravel()
   totals=np.asarray(read_elem(h['obs/total_counts'])) if 'total_counts' in h['obs'] else sums
   assert (totals>0).all() and (totals+1e-6>=sums).all()
   yy=raw[:,[np.flatnonzero(gs==g)[0] for g in present]].toarray()[inds]/totals[inds,None]*1e6
  model=train(pb,v,list(jj));out=[fitcoef(model,y) for y in yy]
  coef=np.array([x[0] for x in out]);prop=np.array([x[1] for x in out]);res=np.array([x[2] for x in out])
  d['fibro_coefficient']=coef[:,0];d['fibro_RNA_weight']=prop[:,0];d['fit_residual']=res
  rec.update(spots=len(d),nonzero_fibro=int((coef[:,0]>0).sum()),median_fibro_RNA_weight=float(np.median(prop[:,0])),median_fit_residual=float(np.median(res)))
  quals.append(rec);d[['spot_id','fibro_coefficient','fibro_RNA_weight','fit_residual']].to_csv(P/'tables'/f'{r.sample}_FIBRO_PILOT.tsv.gz',sep='\t',index=False)
  np.savez_compressed(P/'inputs'/f'{r.sample}_SPATIAL_PANEL.npz',spot_id=d.spot_id.to_numpy(str),genes=np.array(present),CPM=yy,coefficients=coef,weights=prop,residual=res,totals=totals[inds])
  if not bool(r.physical_eligible):continue
  xy=d[['physical_x_um','physical_y_um']].to_numpy();tree=cKDTree(xy)
  for radius in [200,400]:
   lists=tree.query_ball_point(xy,radius);ii=[];jjj=[]
   for k,lst in enumerate(lists):
    neigh=[j for j in lst if j!=k];ii.extend([k]*len(neigh));jjj.extend(neigh)
   A=sparse.csr_matrix((np.ones(len(ii)),(ii,jjj)),shape=(len(d),len(d)));num=np.asarray(A.sum(1)).ravel()
   W=sparse.diags(1/np.maximum(num,1))@A
   cols=['log1p_UMI','total_macro_fraction','malignant_fraction','endothelial_fraction','pericyte_fraction','Tcell_fraction','SPP1_RNA','fibro_coefficient']
   neigh=W@d[cols].to_numpy()
   zz=d[['spot_id','TYMP_RNA','SPP1_RNA','CD8_expression']].copy()
   for k,c in enumerate(cols):zz['near_'+c]=neigh[:,k]
   good=num>=3;zz['neighbors']=num;zz['eligible']=good
   rf=stats.rankdata(zz.near_fibro_coefficient)/(len(d)+1);rs=stats.rankdata(zz.near_SPP1_RNA)/(len(d)+1);joint=rf*rs
   zz['joint_neighborhood_rank']=joint
   core=np.column_stack([d[['log1p_UMI','total_macro_fraction','malignant_fraction']].to_numpy(),zz[['near_log1p_UMI','near_total_macro_fraction','near_malignant_fraction']].to_numpy()])
   vascular=np.column_stack([core,d[['endothelial_fraction','pericyte_fraction']].to_numpy(),zz[['near_endothelial_fraction','near_pericyte_fraction']].to_numpy()])
   for mod,cv in [('core',core),('vascular_sensitivity',vascular)]:
    for end,x,y,fam in [('SPP1_near_fibro',d.SPP1_RNA.to_numpy(),zz.near_fibro_coefficient.to_numpy(),'organization'),('TYMP_near_fibro',d.TYMP_RNA.to_numpy(),zz.near_fibro_coefficient.to_numpy(),'organization'),('joint_neighborhood_CD8',joint,d.CD8_expression.to_numpy(),'immune_distribution')]:
     cov=np.column_stack([cv,zz.near_Tcell_fraction]) if fam=='immune_distribution' else cv
     eff,n,reason=partial(x,y,cov,good)
     effects.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,radius=radius,model=mod,family=fam,endpoint=end,effect=eff,spots=n,reason=reason))
   zz.to_csv(P/'inputs'/f'{r.sample}_NEIGHBORHOOD_{radius}.tsv.gz',sep='\t',index=False)
  print('SPATIAL_DONE',r.sample,flush=True)
 tab(quals,'SPATIAL_REFERENCE_QUALIFICATION.tsv');ef=pd.DataFrame(effects);tab(ef,'NEIGHBORHOOD_SECTION_EFFECTS.tsv')
 keys=['dataset','radius','model','family','endpoint'];patients=[]
 for key,g in ef.groupby(keys+['patient']):
  values=g.effect.dropna().to_numpy();value=float(np.tanh(np.mean(np.arctanh(np.clip(values,-.999999,.999999))))) if len(values) else np.nan
  patients.append(dict(zip(keys+['patient'],key),effect=value,sections=';'.join(g.loc[g.effect.notna(),'sample']),n_sections=len(values)))
 pt=pd.DataFrame(patients);tab(pt,'NEIGHBORHOOD_PATIENT_EFFECTS.tsv')
 summaries=[]
 for key,g in pt.groupby(keys):
  vals=g.effect.dropna().to_numpy();n=len(vals);p=float(stats.wilcoxon(vals,alternative='two-sided').pvalue) if n>=10 and np.any(vals!=0) else (1. if n>=10 else np.nan)
  summaries.append(dict(zip(keys,key),patients=n,positive=int((vals>0).sum()),negative=int((vals<0).sum()),median_effect=float(np.median(vals)) if n else np.nan,P=p))
 su=pd.DataFrame(summaries);su['FDR']=np.nan
 for _,g in su.groupby(['dataset','radius','model','family']):
  valid=g.P.dropna();order=valid.sort_values();q=np.minimum.accumulate((order.to_numpy()*len(order)/np.arange(1,len(order)+1))[::-1])[::-1] if len(order) else []
  su.loc[order.index,'FDR']=np.minimum(q,1.)
 tab(su,'NEIGHBORHOOD_COHORT_SUMMARY.tsv')
 loo=[]
 for key,g in pt.groupby(keys):
  g=g[g.effect.notna()]
  for patient in g.patient:
   vals=g.loc[g.patient!=patient,'effect'];loo.append(dict(zip(keys,key),omitted_patient=patient,patients=len(vals),median_effect=vals.median(),positive=int((vals>0).sum())))
 tab(loo,'NEIGHBORHOOD_LEAVE_PATIENT.tsv');dump(P/'qc/SPATIAL_INPUT_MANIFEST.json',manifest)
 dump(P/'qc/SPATIAL_PILOT_COMPLETE.json',dict(status='COMPUTED',sections_qualified=sum(x['qualified'] for x in quals),section_effects=len(ef),patient_effects=len(pt),new_tumor_stroma_boundary=False,interpretation='exploratory_neighborhood_RNA_association_not_barrier_crossing',figures=0))
 print('SPATIAL_PILOT_COMPLETE',flush=True)
if __name__=='__main__':main()
