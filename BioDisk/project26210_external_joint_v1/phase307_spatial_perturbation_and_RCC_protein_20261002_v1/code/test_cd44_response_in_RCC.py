from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='4'
import json,datetime,hashlib,math,subprocess,gc
import numpy as np,pandas as pd,anndata as ad
from scipy.stats import rankdata,t as td,false_discovery_control
from scipy.linalg import orth
P=Path(__file__).resolve().parents[1];Q=P.parent/'phase305_larger_immune_sequential_20261002_v1'
proto={'created':datetime.datetime.now().astimezone().isoformat(),'source':'Cell2026 SPAC-seq official Supplementary Table6 (mmc6.xlsx), CD44_KO_vs_Ctrl_Coculture and CD44_KO_vs_Ctrl_Static. Author DESeq2 gene-level effects/P retained, not recalculated independent biological treatment tests. Original small replicate experimental context.',
'frozen_gene_rule':'author pvalue<0.05 AND abs(log2FoldChange)>=0.5 AND baseMean>=50. Require MGI homology class with exactly1 mouse symbol and1 human symbol, unique mapping globally. Exclude CD44,SPP1,TYMP and recipient identity genes CD3D/CD3E/TRAC/CD8A/CD8B to avoid defining marker effects. Do not require FDR for selection. All rejected genes/reasons preserved.',
'score':'CD44-loss response=mean log1pCP10K over KO-up orthologs minus mean over KO-down orthologs. Positive is KO-like, not validated SPP1 response, exhaustion score, or treatment prediction. Co-culture main; static separate sensitivity. No target-driven weights or RCC-based gene selection.',
'RCC_qualification':'Phase305 frozen qualified CD8>=100/patient/timepoint. >=80%genes measured in every contributing library for each direction and >=50%full directional panel detected in>=5%CD8. Both directional panels required, >=10genes per direction. Macro source TYMP/SPP1 qualification unchanged.',
'models':['batch_stage_adjusted','batch_stage_CD8depth_adjusted'],'adjustment':'OBRADOVIC authorbatch andpT; depth sensitivity medianCD8UMI. HCRN BL/Progression separate, no pooled repeat donors. Expected SPP1-to-KO-like score negative; TYMP secondary.',
'inference':'Independent patient rank partial correlations; n>=10 anddf>=5 two-sidedt approximation. Family4 per cohort/timepoint/source:2programsx2models, complete rawP/FDR. n<10effects only. Cell/bins never n. All reverse/zero/notqualified results retained.',
'qualification_relation':'Separate from spatial bin analysis. Uses pure-CD8 external experimental gene list in preexisting RCC scRNA. Same RCC donors asPhase305/306: exploratory reanalysis, not independent cohort replication. CD44 has multiple ligands; response is not SPP1-specific.'}
pf=P/'resources/CD44_RESPONSE_PROTOCOL.json';assert not pf.exists();pf.write_text(json.dumps(proto,indent=2))
checks=[]
def ck(name,ok):
 checks.append({'check':name,'passed':bool(ok)})
 if not ok:raise AssertionError(name)
def eq(a,b):return np.allclose(a,b,atol=1e-10,equal_nan=True)
h=pd.read_csv(P/'resources/HOM_MouseHumanSequence.rpt',sep='\t',dtype=str);maps=[]
for k,g in h.groupby('DB Class Key'):
 mm=g.loc[g['NCBI Taxon ID'].eq('10090'),'Symbol'].drop_duplicates();hh=g.loc[g['NCBI Taxon ID'].eq('9606'),'Symbol'].drop_duplicates()
 if len(mm)==len(hh)==1:maps.append({'homology_class':k,'mouse_gene':mm.iloc[0],'human_gene':hh.iloc[0]})
orths=pd.DataFrame(maps);orths=orths[~orths.mouse_gene.duplicated(False)&~orths.human_gene.duplicated(False)];orths.to_csv(P/'tables/MGI_ONE_TO_ONE_ORTHOLOGS.tsv',sep='\t',index=False);mapping=orths.set_index('mouse_gene').human_gene.to_dict();exclude={'CD44','SPP1','TYMP','CD3D','CD3E','TRAC','CD8A','CD8B'};panels={};selection=[]
for name,sheet in [('Coculture','CD44_KO_vs_Ctrl_Coculture'),('Static','CD44_KO_vs_Ctrl_Static')]:
 d=pd.read_excel(P/'downloads/SPAC_mmc6.xlsx',sheet_name=sheet);ck(name+'_unique_genes',d.Gene_name.is_unique)
 for rr in d.to_dict('records'):
  human=mapping.get(rr['Gene_name']);reasons=[]
  if not rr['pvalue']<.05:reasons.append('AUTHOR_RAWP')
  if not abs(rr['log2FoldChange'])>=.5:reasons.append('EFFECT_SIZE')
  if not rr['baseMean']>=50:reasons.append('EXTERNAL_MEASUREMENT')
  if human is None:reasons.append('NO_UNIQUE_ONE_TO_ONE_ORTHOLOG')
  if human in exclude:reasons.append('TARGET_OR_IDENTITY_EXCLUSION')
  selection.append({'program':name,**rr,'human_gene':human,'selected':not reasons,'direction':'up' if rr['log2FoldChange']>0 else 'down','reason':';'.join(reasons)})
 selected=pd.DataFrame(selection);sub=selected[selected.program.eq(name)&selected.selected];panels[name]={sign:sub.loc[sub.direction.eq(sign),'human_gene'].tolist() for sign in ['up','down']};ck(name+'_both_directions',all(len(v)>=10 for v in panels[name].values()));print('FROZEN_PANEL',name,{k:len(v) for k,v in panels[name].items()},flush=True)
pd.DataFrame(selection).to_csv(P/'tables/EXTERNAL_CD44_GENE_SELECTION_ALL.tsv',sep='\t',index=False);(P/'resources/CD44_LOSS_PANELS_FROZEN.json').write_text(json.dumps(panels,indent=2));interest=sorted(set(sum([v for d in panels.values() for v in d.values()],[])));acc={};meta_records={}
for ds in ['OBRADOVIC','HCRN']:
 m=pd.read_csv(Q/f'tables/{ds}_QUALIFIED_CELL_METADATA.tsv.gz',sep='\t',index_col=0)
 for lib in m.loc[m.lineage_CD8,'library'].unique():
  subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],check=True,stdout=subprocess.DEVNULL);a=ad.read_h5ad(Q/f'inputs/{ds}_{lib}_qc.h5ad');mm=m.loc[a.obs_names];mask=mm.lineage_CD8.to_numpy(bool);xx=a.X[mask].tocsr();n=xx.shape[0];unit=mm.analysis_unit.iloc[0];idx=a.var_names.get_indexer(interest)
  if unit not in acc:acc[unit]={'n':0,'total':0.,'sum':np.zeros(len(interest)),'det':np.zeros(len(interest)),'available':set(interest),'umi':[]};meta_records[unit]={'dataset':ds,'analysis_unit':unit,'patient':mm.patient.iloc[0],'stratum':mm.stratum.iloc[0],'technical_batch':mm.technical_batch.iloc[0]}
  q=acc[unit];q['n']+=n;q['total']+=float(xx.sum());q['sum'][idx>=0]+=np.asarray(xx[:,idx[idx>=0]].sum(0)).ravel();q['det'][idx>=0]+=np.asarray((xx[:,idx[idx>=0]]>0).sum(0)).ravel();q['available']&=set(a.var_names);q['umi'].extend(np.asarray(xx.sum(1)).ravel().tolist());del a,xx;gc.collect()
 print('RCC_AGGREGATED',ds,flush=True)
rows=[];cover=[];index={g:i for i,g in enumerate(interest)}
for unit,q in acc.items():
 base=meta_records[unit]|{'cells':q['n'],'CD8_median_UMI':float(np.median(q['umi']))};pb=np.log1p(q['sum']*1e4/max(q['total'],1))
 for program,parts in panels.items():
  values={};reason=[]
  if q['n']<100:reason.append('CD8_CELL_GATE')
  for sign,gs in parts.items():
   measured=[g for g in gs if g in q['available']];detect=sum(q['det'][index[g]]>=.05*q['n'] for g in measured)
   if len(measured)<.8*len(gs):reason.append(sign+'_MEASUREMENT')
   if detect<.5*len(gs):reason.append(sign+'_DETECTION')
   values[sign]=np.mean([pb[index[g]] for g in measured]) if measured else np.nan
   calc=math.fsum(math.log1p(q['sum'][index[g]]*1e4/max(q['total'],1)) for g in measured)/len(measured) if measured else np.nan;ck(unit+program+sign+'_independent_score',eq(calc,values[sign]))
   for gene in gs:cover.append(base|{'program':program,'direction':sign,'gene':gene,'measured_all_libraries':gene in q['available'],'detected_cells':q['det'][index[gene]],'total_gene_counts':q['sum'][index[gene]],'total_RNA_counts':q['total'],'logCP10K':pb[index[gene]]})
  rows.append(base|{'program':program,'up_score':values['up'],'down_score':values['down'],'value':values['up']-values['down'],'eligible':not reason,'reason':';'.join(reason)})
d=pd.DataFrame(rows);d.to_csv(P/'tables/CD44_RESPONSE_PATIENT_SCORES.tsv',sep='\t',index=False);pd.DataFrame(cover).to_csv(P/'tables/CD44_RESPONSE_GENE_COVERAGE.tsv',sep='\t',index=False)
old=pd.read_csv(Q/'tables/PATIENT_MEASURES.tsv',sep='\t');sources=old[old.lineage.eq('Macrophage')&old.eligible];effects=[];pairs=[]
for ds,st in [('OBRADOVIC','Untreated'),('HCRN','Baseline'),('HCRN','Progression')]:
 for src in ['SPP1','TYMP']:
  ss=sources[sources.dataset.eq(ds)&sources.stratum.eq(st)&sources.measure.eq(src)][['analysis_unit','value','pT_category']].rename(columns={'value':'source_value'})
  for program in panels:
   g=d[d.dataset.eq(ds)&d.stratum.eq(st)&d.program.eq(program)&d.eligible].merge(ss,on='analysis_unit',validate='one_to_one');n=len(g);ck(ds+st+src+program+'_unique_patients',g.patient.nunique()==n)
   for rr in g.to_dict('records'):pairs.append(rr|{'source':src})
   for model in proto['models']:
    cols=[np.ones(n)]
    if ds=='OBRADOVIC':cols.extend(pd.get_dummies(g[['technical_batch','pT_category']],drop_first=True,dtype=float).to_numpy().T)
    if model.endswith('CD8depth_adjusted') and n and np.ptp(g.CD8_median_UMI)>0:cols.append(rankdata(g.CD8_median_UMI))
    D=np.column_stack(cols);rank=int(np.linalg.matrix_rank(D)) if n else 0;df=n-rank-1;rho=np.nan;pval=np.nan
    if n>=3:
     xy=np.column_stack([rankdata(g.source_value),rankdata(g.value)]);res=xy-D@np.linalg.lstsq(D,xy,rcond=None)[0]
     if min(np.std(res,axis=0))>1e-9:
      rho=float(np.corrcoef(res.T)[0,1]);base=orth(D);test=xy-base@base.T@xy;ck(ds+st+src+program+model+'_orthogonal_rho',eq(rho,test[:,0]@test[:,1]/np.sqrt((test[:,0]@test[:,0])*(test[:,1]@test[:,1]))))
     if n>=10 and df>=5 and np.isfinite(rho):pval=float(2*td.sf(abs(rho)*np.sqrt(df/max(1-rho*rho,1e-15)),df))
    effects.append({'dataset':ds,'stratum':st,'source':src,'program':program,'model':model,'patients':n,'df':df,'rho':rho,'P':pval,'family_planned_n':4,'expected_direction':'negative','reason':'' if np.isfinite(pval) else 'SMALL_N_OR_NOT_ESTIMABLE'})
e=pd.DataFrame(effects);e['FDR']=np.nan;e['family_tested_n']=0
for key,ids in e.groupby(['dataset','stratum','source']).groups.items():
 valid=[i for i in ids if pd.notna(e.at[i,'P'])];e.loc[list(ids),'family_tested_n']=len(valid)
 if valid:
  pv=e.loc[valid,'P'].to_numpy();q=false_discovery_control(pv);order=np.argsort(pv);manual=np.empty(len(pv));manual[order]=np.minimum(np.minimum.accumulate((pv[order]*len(pv)/np.arange(1,len(pv)+1))[::-1])[::-1],1);ck(str(key)+'_BH',eq(q,manual));e.loc[valid,'FDR']=q
e.to_csv(P/'tables/CD44_RESPONSE_ASSOCIATIONS_ALL.tsv',sep='\t',index=False);pd.DataFrame(pairs).to_csv(P/'tables/CD44_RESPONSE_ASSOCIATION_PAIRS.tsv',sep='\t',index=False);ck('complete24models',len(e)==24);ck('noHCRN_P',e.loc[e.dataset.eq('HCRN'),'P'].isna().all());(P/'qc/CD44_RESPONSE_VALIDATION_CHECKS.json').write_text(json.dumps(checks,indent=2));(P/'qc/CD44_RESPONSE_COMPLETE.json').write_text(json.dumps({'status':'PASS','checks':len(checks),'time':datetime.datetime.now().astimezone().isoformat(),'protocol_sha256':hashlib.sha256(pf.read_bytes()).hexdigest(),'source_xlsx_sha256':hashlib.sha256((P/'downloads/SPAC_mmc6.xlsx').read_bytes()).hexdigest(),'orthology_sha256':hashlib.sha256((P/'resources/HOM_MouseHumanSequence.rpt').read_bytes()).hexdigest()},indent=2));print('RCC_RESPONSE_RESULTS',e[e.dataset.eq('OBRADOVIC')].to_string(index=False),flush=True)
