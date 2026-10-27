from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,subprocess,datetime
P=Path(__file__).resolve().parents[1];R=P.parent
def gate(s):
 q=subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True)
 with (P/'logs/qualification_memory.log').open('a') as f:f.write(s+'\n'+q.stdout+q.stderr)
 if q.returncode:raise RuntimeError('MEMORY_GATE')
def sha(f):
 h=hashlib.sha256()
 with open(f,'rb') as s:
  for b in iter(lambda:s.read(1048576),b''):h.update(b)
 return h.hexdigest()
def read(f):return pd.read_csv(f,sep='\t',float_precision='round_trip')
def table(v,n):pd.DataFrame(v).to_csv(P/'tables'/n,sep='\t',index=False)
def main():
 assert (P/'qc/EXTRACTION_COMPLETE').exists()
 ids=json.loads((R/'phase286_CAF_identity_boundary_audit_20260929_v1/resources/IDENTITY_PANELS_FROZEN.json').read_text())
 panels=json.loads((R/'phase292_barrier_execution_and_external_feasibility_20260930_v1/PANELS_FROZEN.json').read_text())
 table([dict(program=k,gene=g) for k,v in panels.items() for g in v],'FROZEN_FUNCTION_GENES.tsv')
 rows=[];ident=[];gene_rows=[];lineage=[];manifest=[]
 for ds in ['GSE264586','GSE294109']:
  gate(ds)
  fs=[P/'inputs'/f'{ds}_METADATA.tsv.gz',P/'inputs'/f'{ds}_PANEL_COUNTS.tsv.gz']
  manifest += [dict(path=str(f),sha256=sha(f)) for f in fs]
  m=read(fs[0]);x=read(fs[1]);assert x.cell_id.is_unique and m.cell_id.tolist()==x.cell_id.tolist()
  x=x.drop(columns='cell_id');assert np.isfinite(x.to_numpy()).all() and (x.to_numpy()>=0).all() and (x.to_numpy()==np.floor(x.to_numpy())).all()
  assert set(sum(panels.values(),[])).issubset(x)
  det={k:(x[[g for g in gs if g in x]]>0).sum(axis=1) for k,gs in ids.items()}
  broad=m.phase7_broad_annotation.eq('Fibroblast')
  core=(det['fibro_core']>=2)&(det['fibro_support']>=1)
  competing=(det['immune']>=2)|(det['endothelial']>=2)|(det['epithelial']>=2)
  mural=(det['pericyte']>=2)|((x[['MYH11','CNN1']]>0).any(axis=1)&(x[['ACTA2','TAGLN']]>0).any(axis=1))
  strict=broad&core&~competing&~mural
  label=np.where(~broad,'not_old_fibro',np.where(competing,'competing_lineage',np.where(mural,'mural_overlap',np.where(core,'strict_fibro_candidate','insufficient_marker'))))
  mm=m[['cell_id','patient_key','sample_id','phase7_broad_annotation','phase7_fine_annotation']].copy()
  mm['qualification_label']=label;mm['strict_fibro']=strict;mm.to_csv(P/'inputs'/f'{ds}_IDENTITY_FREEZE.tsv.gz',sep='\t',index=False)
  for (sample,cl),g in mm.groupby(['sample_id','qualification_label']):
   ident.append(dict(dataset=ds,sample_id=sample,label=cl,cells=len(g)))
  # marker evidence per existing broad/fine group, without re-clustering
  for (b,fi),g in m.groupby(['phase7_broad_annotation','phase7_fine_annotation']):
   for gn in sorted(set(sum(ids.values(),[]))):
    if gn in x:
     lineage.append(dict(dataset=ds,broad=b,fine=fi,gene=gn,cells=len(g),detected=int((x.loc[g.index,gn]>0).sum()),mean_raw=float(x.loc[g.index,gn].mean())))
  for (pat,sample),g in m.groupby(['patient_key','sample_id']):
   idx=g.index
   treatments=g.treatment.fillna('Unknown').unique();times=g.timepoint.fillna('Unknown').unique()
   assert len(treatments)==len(times)==1
   stratum=times[0] if ds=='GSE264586' else treatments[0]
   for receiver,mask in [('strict_fibro',strict),('old_fibro_pool',broad)]:
    fi=idx[mask.loc[idx]];ep=idx[m.loc[idx,'phase7_fine_annotation'].eq('Epithelial_tumor_like')];ma=idx[m.loc[idx,'phase7_fine_annotation'].isin(['Macrophage_C1QC','Macrophage_SPP1'])]
    out=dict(dataset=ds,patient=pat,library=sample,treatment=treatments[0],timepoint=times[0],stratum=stratum,receiver=receiver,receiver_cells=len(fi),tumor_like_cells=len(ep),macro_cells=len(ma))
    total=float(m.loc[fi,'total_counts'].sum());v=x.loc[fi].sum(axis=0).astype(float);score=np.log1p(v/max(total,1)*1e4)
    out['receiver_depth']=float(np.log1p(total/len(fi))) if len(fi) else np.nan
    out['receiver_total_UMI']=total
    for src,ii in [('tumor_like',ep),('macro',ma)]:
     tt=float(m.loc[ii,'total_counts'].sum());vv=x.loc[ii,['TYMP','SPP1']].sum().astype(float)
     for gn in ['TYMP','SPP1']:out[src+'_'+gn]=float(np.log1p(vv[gn]/tt*1e4)) if tt else np.nan
    for name,gs in panels.items():
     detected=int((v[gs]>0).sum());out[name+'_detected_genes']=detected
     out[name]=float(score[gs].mean()) if detected>=np.ceil(len(gs)/2) and total else np.nan
     for gn in gs:gene_rows.append(dict(dataset=ds,patient=pat,library=sample,receiver=receiver,program=name,gene=gn,counts=float(v[gn]),logCP10K=float(score[gn]),total_UMI=total))
    out['qualified']=len(ep)>=50 and len(ma)>=100 and len(fi)>=(20 if receiver=='strict_fibro' else 50)
    out['qualification_reason']='OK' if out['qualified'] else ';'.join(n for n,test in [('tumor_like_lt50',len(ep)>=50),('macro_lt100',len(ma)>=100),('receiver_below_gate',len(fi)>=(20 if receiver=='strict_fibro' else 50))] if not test)
    rows.append(out)
 table(rows,'EXTERNAL_LIBRARY_INPUTS.tsv');table(ident,'IDENTITY_CAPACITY.tsv');table(gene_rows,'EXTERNAL_GENE_PSEUDOBULK.tsv');table(lineage,'ORIGINAL_LINEAGE_MARKER_EVIDENCE.tsv')
 df=pd.DataFrame(rows);pa=[]
 for (ds,pat,stratum,recv),g in df[df.qualified].groupby(['dataset','patient','stratum','receiver']):
  for program in panels:
   z=g[g[program].notna()]
   if z.empty:continue
   pa.append(dict(dataset=ds,patient=pat,stratum=stratum,receiver=recv,endpoint=program,libraries=';'.join(z.library),n_libraries=len(z),score=z[program].mean(),receiver_depth=z.receiver_depth.mean(),tumor_like_TYMP=z.tumor_like_TYMP.mean(),macro_TYMP=z.macro_TYMP.mean(),macro_SPP1=z.macro_SPP1.mean()))
 table(pa,'EXTERNAL_PATIENT_FREEZE.tsv')
 manifest += [dict(path=str(P/'EXTERNAL_IDENTITY_AND_STATISTICS_PROTOCOL.md'),sha256=sha(P/'EXTERNAL_IDENTITY_AND_STATISTICS_PROTOCOL.md'))]
 (P/'qc/EXTERNAL_INPUT_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
 (P/'qc/PATIENT_FREEZE.json').write_text(json.dumps(dict(at=datetime.datetime.now().astimezone().isoformat(),hash=sha(P/'tables/EXTERNAL_PATIENT_FREEZE.tsv'),association_not_yet_run=True),indent=2))
 print(df.groupby(['dataset','receiver']).agg(libraries=('library','size'),qualified=('qualified','sum')).to_string())
 if pa:print(pd.DataFrame(pa).groupby(['dataset','stratum','receiver']).patient.nunique().to_string())
 print('QUALIFICATION_COMPLETE')
if __name__=='__main__':main()
