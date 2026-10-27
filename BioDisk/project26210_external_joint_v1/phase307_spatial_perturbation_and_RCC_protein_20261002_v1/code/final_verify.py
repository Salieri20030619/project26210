from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='4'
import numpy as np,pandas as pd,anndata as ad,json,datetime,subprocess,openpyxl
from scipy.stats import t as td
P=Path(__file__).resolve().parents[1];Q=P.parent/'phase305_larger_immune_sequential_20261002_v1';checks=[]
def ck(name,ok):
 checks.append({'check':name,'passed':bool(ok)})
 if not ok:raise AssertionError(name)
for name in ['download_ot1','spatial_analysis','spatial_validation','guide_localization','cd44_response','cd44_pseudobulk_v2']:ck(name+'_exit',(P/f'logs/{name}.exit').read_text().strip()=='0')
panels=json.loads((P/'resources/CD44_LOSS_PANELS_FROZEN.json').read_text());cv=pd.read_csv(P/'tables/CD44_RESPONSE_GENE_COVERAGE.tsv',sep='\t');scores=pd.read_csv(P/'tables/CD44_PSEUDOBULK_PATIENT_SCORES_V2.tsv',sep='\t');interest=sorted(cv.gene.unique());gi={g:i for i,g in enumerate(interest)}
for ds in ['OBRADOVIC','HCRN']:
 subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],check=True,stdout=subprocess.DEVNULL)
 a=ad.read_h5ad(Q/f'inputs/{ds}_counts_qc.h5ad',backed='r');m=pd.read_csv(Q/f'tables/{ds}_QUALIFIED_CELL_METADATA.tsv.gz',sep='\t',index_col=0);ii=np.flatnonzero(m.lineage_CD8.to_numpy(bool));x=a[ii,:].X.tocsr();meta=m.iloc[ii];gene_idx=a.var_names.get_indexer(interest);present=gene_idx>=0;missing=set(np.array(interest)[~present]);ck(ds+'_missing_not_claimed_measured',not cv.loc[cv.dataset.eq(ds)&cv.measured_all_libraries,'gene'].isin(missing).any());a.file.close();tot=np.asarray(x.sum(1)).ravel();target=x[:,gene_idx[present]]
 for unit,g in meta.groupby('analysis_unit'):
  mask=meta.analysis_unit.eq(unit).to_numpy();sums=np.zeros(len(interest));detect=np.zeros(len(interest));sums[present]=np.asarray(target[mask].sum(0)).ravel();detect[present]=np.asarray((target[mask]>0).sum(0)).ravel();total=tot[mask].sum();rows=cv[cv.analysis_unit.eq(unit)]
  ids=np.array([gi[z] for z in rows.gene]);ck(unit+'_joined_vs_original_library_gene_sums',np.array_equal(sums[ids],rows.total_gene_counts.to_numpy()));ck(unit+'_joined_vs_original_library_detection',np.array_equal(detect[ids],rows.detected_cells.to_numpy()));ck(unit+'_total_RNA',np.allclose(rows.total_RNA_counts,total))
  for rr in scores[scores.analysis_unit.eq(unit)].itertuples():
   vals={};ok=mask.sum()>=100
   for sign in ['up','down']:
    q=rows[rows.program.eq(rr.program)&rows.direction.eq(sign)];measured=q[q.measured_all_libraries];jj=[gi[z] for z in measured.gene];vals[sign]=np.log1p(sums[jj]*1e4/total).mean();ok=ok and len(measured)>=.8*len(q) and sum(detect[gi[z]]>=5 for z in measured.gene)>=.5*len(q)
   ck(unit+rr.program+'_pseudobulk_score',np.isclose(vals['up']-vals['down'],rr.value,atol=1e-10));ck(unit+rr.program+'_qualification',bool(ok)==rr.eligible)
 print('FINAL_INPUT_VERIFIED',ds,flush=True)
e=pd.read_csv(P/'tables/CD44_PSEUDOBULK_ASSOCIATIONS_V2.tsv',sep='\t');pairs=pd.read_csv(P/'tables/CD44_PSEUDOBULK_PAIRS_V2.tsv',sep='\t')
for i,r in e.iterrows():
 g=pairs[pairs.dataset.eq(r.dataset)&pairs.stratum.eq(r.stratum)&pairs.source.eq(r.source)&pairs.program.eq(r.program)];ck(f'model{i}_N',len(g)==g.patient.nunique()==r.patients)
 if pd.notna(r.P):ck(f'model{i}_P',np.isclose(2*td.sf(np.sqrt(r.df)*abs(r.rho)/np.sqrt(1-r.rho*r.rho),r.df),r.P,atol=1e-12) and r.patients>=10 and r.df>=5)
for key,g in e.groupby(['dataset','stratum','source']):
 v=g[g.P.notna()].sort_values('P');n=len(v)
 if n:ck(str(key)+'_FDR',np.allclose(np.minimum(np.minimum.accumulate((v.P.to_numpy()*n/np.arange(1,n+1))[::-1])[::-1],1),v.FDR,atol=1e-12))
w=openpyxl.load_workbook(P/'downloads/RCC_SourceData.xlsx',read_only=False,data_only=True);ck('RCC_1G_has_ID','ID' in [c.value for c in w['1G'][1]]);ck('RCC_5C_lacks_ID','ID' not in [c.value for c in w['5C'][1]] and w['5C'].max_column==4);ck('RCC_5C_no_hidden_ID_columns',not any(v.hidden for v in w['5C'].column_dimensions.values()))
all_checks={}
for name,key in [('FINAL_VALIDATION.json','checks'),('LOCALIZATION_COMPLETE.json','checks'),('CD44_RESPONSE_COMPLETE.json','checks'),('CD44_PSEUDOBULK_VALIDATION_V2.json','check_n')]:
 j=json.loads((P/'qc'/name).read_text());ck(name+'_PASS',j['status']=='PASS');all_checks[name]=j[key]
all_checks['independent_final']=len(checks);total=sum(all_checks.values());(P/'qc/FINAL_PHASE307_CHECKS.json').write_text(json.dumps(checks,indent=2));(P/'qc/PHASE307_FINAL_VALIDATION.json').write_text(json.dumps({'status':'PASS','checks':total,'components':all_checks,'time':datetime.datetime.now().astimezone().isoformat(),'scope':'24officialMD5 + original spatial H5 + independent spatial geometry + independently joined RCC counts + complete qualification and all statistics'},indent=2));print('PHASE307_VALIDATED',total,flush=True)
