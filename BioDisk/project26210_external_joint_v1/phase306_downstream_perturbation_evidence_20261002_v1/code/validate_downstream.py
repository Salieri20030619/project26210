from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='4'
import json,subprocess,hashlib,datetime
import numpy as np,pandas as pd,anndata as ad,openpyxl
from scipy.linalg import orth
from scipy.stats import t as td
P=Path(__file__).resolve().parents[1];Q=P.parent/'phase305_larger_immune_sequential_20261002_v1';checks=[]
def ck(n,ok):
 checks.append({'check':n,'passed':bool(ok)})
 if not ok:raise AssertionError(n)
def eq(a,b):return np.allclose(a,b,atol=1e-10,equal_nan=True)
def bh(x,m=None):
 x=np.asarray(x);m=len(x) if m is None else m;o=np.argsort(x);v=np.minimum.accumulate((x[o]*m/np.arange(1,len(x)+1))[::-1])[::-1];out=np.empty(len(x));out[o]=np.minimum(v,1);return out
ck('analysis_exit',(P/'logs/analyze_downstream.exit').read_text().strip()=='0')
ck('source_sha',hashlib.sha256((P/'downloads/Trehan_SourceData.xlsx').read_bytes()).hexdigest()=='0d6e1b8b9726e8e457331a20a32e9e19713b4623290659b86bb4c28e928fbaf7')
d=pd.read_csv(P/'tables/RECIPIENT_PATIENT_ENDPOINTS.tsv',sep='\t');co=pd.read_csv(P/'tables/RECIPIENT_GENE_COVERAGE.tsv',sep='\t');acc={};genes=['PDCD1','ENTPD1','TOX','IFNG','GZMB','CD44']
for ds in ['OBRADOVIC','HCRN']:
 m=pd.read_csv(Q/f'tables/{ds}_QUALIFIED_CELL_METADATA.tsv.gz',sep='\t',index_col=0)
 for lib in m.loc[m.lineage_CD8,'library'].unique():
  subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],check=True,stdout=subprocess.DEVNULL)
  a=ad.read_h5ad(Q/f'inputs/{ds}_{lib}_qc.h5ad');mm=m.loc[a.obs_names];xx=a.X[mm.lineage_CD8.to_numpy(bool)].tocsr();n=xx.shape[0];unit=mm.analysis_unit.iloc[0];jj=a.var_names.get_indexer(genes);ck(lib+'_all_genes_measured',np.all(jj>=0));z=xx[:,jj].toarray()
  if unit not in acc:acc[unit]={'n':0,'tot':0,'sum':np.zeros(6),'det':np.zeros(6),'double':0,'umi':[]}
  q=acc[unit];q['n']+=n;q['tot']+=float(xx.sum());q['sum']+=z.sum(0);q['det']+=(z>0).sum(0);q['double']+=int(((z[:,0]>0)&(z[:,1]>0)).sum());q['umi'].extend(np.asarray(xx.sum(1)).ravel().tolist())
for i,r in d.iterrows():
 q=acc[r.analysis_unit];vals={'PDCD1_ENTPD1_double_detected_fraction':q['double']/q['n'],'TOX_detected_fraction':q['det'][2]/q['n'],'IFNG_GZMB_effector_mean_logCP10K':np.log1p(q['sum'][3:5]*1e4/q['tot']).mean()};ck(f'endpoint{i}_original_counts',eq(vals[r.endpoint],r.value));ck(f'endpoint{i}_n',q['n']==r.cells);ck(f'endpoint{i}_depth',eq(np.median(q['umi']),r.CD8_median_UMI));eligible=q['n']>=100 and (not r.endpoint.startswith('IFNG') or (q['det'][3:5]>=.05*q['n']).sum()>=1);ck(f'endpoint{i}_eligibility',eligible==r.eligible)
for i,r in co.iterrows():
 q=acc[r.analysis_unit];j=genes.index(r.gene);ck(f'gene{i}_det',q['det'][j]==r.detected_cells);ck(f'gene{i}_pseudo',eq(np.log1p(q['sum'][j]*1e4/q['tot']),r.pseudobulk_logCP10K))
e=pd.read_csv(P/'tables/RECIPIENT_ASSOCIATIONS_ALL.tsv',sep='\t');pairs=pd.read_csv(P/'tables/RECIPIENT_ASSOCIATION_PAIRS.tsv',sep='\t');source=pd.read_csv(Q/'tables/PATIENT_MEASURES.tsv',sep='\t');source=source[source.lineage.eq('Macrophage')&source.eligible]
for i,r in e.iterrows():
 g=pairs[pairs.dataset.eq(r.dataset)&pairs.stratum.eq(r.stratum)&pairs.source.eq(r.source)&pairs.endpoint.eq(r.endpoint)];ck(f'assoc{i}_patients',g.patient.nunique()==len(g)==r.patients)
 ss=source[source.dataset.eq(r.dataset)&source.stratum.eq(r.stratum)&source.measure.eq(r.source)].set_index('analysis_unit');dd=d[d.endpoint.eq(r.endpoint)&d.eligible].set_index('analysis_unit');ck(f'assoc{i}_source_values',eq(g.source_value,ss.loc[g.analysis_unit,'value']));ck(f'assoc{i}_recipient_values',eq(g.value,dd.loc[g.analysis_unit,'value']))
 if pd.isna(r.rho):ck(f'assoc{i}_noP',pd.isna(r.P));continue
 cols=[np.ones(len(g))]
 if r.dataset=='OBRADOVIC':cols.extend(pd.get_dummies(g[['technical_batch','pT_category']],drop_first=True,dtype=float).to_numpy().T)
 if r.model.endswith('CD8depth_adjusted'):cols.append(g.CD8_median_UMI.rank().to_numpy())
 D=np.column_stack(cols);basis=orth(D);xy=g[['source_value','value']].rank().to_numpy();rr=xy-basis@basis.T@xy;rho=rr[:,0]@rr[:,1]/np.sqrt((rr[:,0]@rr[:,0])*(rr[:,1]@rr[:,1]));ck(f'assoc{i}_rho',eq(rho,r.rho));ck(f'assoc{i}_df',len(g)-basis.shape[1]-1==r.df)
 if pd.notna(r.P):ck(f'assoc{i}_P',eq(2*td.sf(abs(rho)*np.sqrt(r.df/(1-rho*rho)),r.df),r.P) and len(g)>=10 and r.df>=5)
for key,g in e.groupby(['dataset','stratum','source']):
 v=g[g.P.notna()];ck(str(key)+'_family',len(g)==6 and g.family_tested_n.eq(len(v)).all())
 if len(v):ck(str(key)+'_BH',eq(bh(v.P),v.FDR))
w=openpyxl.load_workbook(P/'downloads/Trehan_SourceData.xlsx',read_only=True,data_only=True)
sp=pd.read_csv(P/'tables/EXTERNAL_SPATIAL_EFFECTS.tsv',sep='\t');sv=pd.read_csv(P/'tables/EXTERNAL_SPATIAL_SOURCE_VALUES.tsv',sep='\t')
column={'6F':3,'6G':7,'6H':11,'6I':15,'6J':19,'6K':23,'6L':3,'6M':7,'6N':11,'6O':15}
for i,r in sv.iterrows():
 c=column[r.panel];ck(f'spatial_value{i}',eq([w['Figure 6'].cell(int(r.source_row),c).value,w['Figure 6'].cell(int(r.source_row),c+1).value],[r.alphaSMA_positive,r.alphaSMA_negative]))
for i,r in sp.iterrows():
 g=sv[sv.panel.eq(r.panel)];delta=g.alphaSMA_positive-g.alphaSMA_negative;ck(f'spatial{i}_effect',eq(delta.mean(),r.mean_paired_difference));ck(f'spatial{i}_directions',(delta<0).sum()==r.negative_direction and (delta>0).sum()==r.positive_direction)
 if pd.notna(r.P):ck(f'spatial{i}_P',eq(2*td.sf(abs(delta.mean()/delta.std(ddof=1)*np.sqrt(len(g))),len(g)-1),r.P))
for key,g in sp.groupby('family'):
 v=g[g.P.notna()];ck('spatial_BH_'+key,eq(bh(v.P,int(g.family_planned_n.iloc[0])),v.FDR))
pr=pd.read_csv(P/'tables/EXTERNAL_STIMULATION_SOURCE_VALUES.tsv',sep='\t');pe=pd.read_csv(P/'tables/EXTERNAL_STIMULATION_EFFECTS.tsv',sep='\t');pc={'IgG':3,'Spp1':4,'Spp1_aCD44':5}
for i,r in pr.iterrows():ck(f'stim_source{i}',eq(w['Figure 4'].cell(int(r.source_row),pc[r.group]).value,r.score))
for i,r in pe.iterrows():
 g=pr[pr.panel.eq(r.panel)].groupby('group').score.mean();ck(f'stim{i}_delta',eq(g.Spp1-g.IgG,r.stimulation_difference));ck(f'stim{i}_attenuation',eq((g.Spp1-g.Spp1_aCD44)/(g.Spp1-g.IgG),r.attenuation_fraction));ck(f'stim{i}_noP',pd.isna(r.P) and pd.isna(r.FDR))
kr=pd.read_csv(P/'tables/EXTERNAL_KO_SOURCE_VALUES.tsv',sep='\t');ke=pd.read_csv(P/'tables/EXTERNAL_KO_EFFECTS.tsv',sep='\t')
for i,r in kr.iterrows():ck(f'KOsource{i}',eq(w['Figure 5'].cell(int(r.source_row),4 if r.population=='TRP2_positive' else 5).value,r.percent))
for i,r in ke.iterrows():
 g=kr[kr.panel.eq(r.panel)&kr.population.eq(r.population)];means=g.groupby('group').percent.mean();ck(f'KOeffect{i}',eq(means.KO-means.WT,r.difference_percentage_points));ck(f'KOnoP{i}',pd.isna(r.new_P) and pd.isna(r.new_FDR))
ck('all36_models',len(e)==36);ck('HCRN_no_inference',e.loc[e.dataset.eq('HCRN'),'P'].isna().all());ck('no_adjusted_positive_evidence',not ((e.P<.05)&((e.endpoint.str.startswith('IFNG')&(e.rho<0))|(~e.endpoint.str.startswith('IFNG')&(e.rho>0)))).any())
(P/'qc/VALIDATION_CHECKS.json').write_text(json.dumps(checks,indent=2));(P/'qc/FINAL_VALIDATION.json').write_text(json.dumps({'status':'PASS','checks':len(checks),'time':datetime.datetime.now().astimezone().isoformat(),'method':'Independent per-library raw-count reaggregation; orthogonal rank residuals; manual BH and paired t formulas; workbook cell provenance'},indent=2));print('VALIDATION_PASS',len(checks),flush=True)
