from pathlib import Path
import pandas as pd,numpy as np,json,subprocess,hashlib,datetime
from scipy.stats import pearsonr,spearmanr,hypergeom,fisher_exact
P=Path(__file__).resolve().parents[1]
protocol=json.loads((P/'resources/REPRODUCTION_PROTOCOL_FROZEN.json').read_text())
def gate():subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],check=True)
def table(s):return pd.read_csv(P/'tables'/('SOURCE_'+s.replace(' ','_')+'.tsv'),sep='\t',float_precision='round_trip')
def dump(n,o):(P/'qc'/n).write_text(json.dumps(o,indent=2,default=str))
gate();checks=[];corr=[];immune=[];merged={};mapping=[]
key=['imagerow','imagecol']
def check(name,ok,detail):checks.append({'check':name,'pass':bool(ok),'detail':detail})
for patient,m in protocol['identity'].items():
 gate();labels=table(m['labels']);check(patient+'_rows',len(labels)==protocol['units'][patient],len(labels))
 check(patient+'_unique_coord',not labels.duplicated(key).any(),len(labels))
 if 'scores' in m:
  sc=table(m['scores']);sc=sc.rename(columns={c:'FAP' if c.startswith('FAP') else 'SPP1' for c in sc if c not in key})
 else:
  a=table(m['FAP']);b=table(m['SPP1']);sc=a.merge(b,on=key,validate='one_to_one')
  sc=sc.rename(columns={c:'FAP' if c.startswith('FAP') else 'SPP1' for c in sc if c not in key})
 z=labels.merge(sc,on=key,how='outer',indicator=True,validate='one_to_one')
 check(patient+'_score_join',z['_merge'].eq('both').all(),z['_merge'].value_counts().to_dict());z=z.drop(columns='_merge')
 check(patient+'_scores_finite',np.isfinite(z[['FAP','SPP1']]).all().all(),'all scores')
 z['joint']=z.DefineTypes.eq(m['joint_label']);z['malignant']=z.DefineTypes.eq('Malignant epithelial cells')
 for scope,b in [('author_joint_cluster',z[z.joint]),('whole_section_diagnostic',z)]:
  rr=pearsonr(b.FAP,b.SPP1);sr=spearmanr(b.FAP,b.SPP1)
  row={'patient':patient,'scope':scope,'n_spots':len(b),'biological_n':1,'pearson_r':float(rr.statistic),'paper_spot_model_P':float(rr.pvalue),'spearman_r_diagnostic':float(sr.statistic),'P_type':'source spot-level model, not biological-patient inference','paper_rounded_R':{'P6':.58,'P7':.49}.get(patient) if scope=='author_joint_cluster' else None}
  corr.append(row)
 ims=table(m['immune']);q=z.merge(ims,on=key,how='outer',indicator=True,validate='one_to_one')
 counts=q['_merge'].value_counts().to_dict();ok=q['_merge'].eq('both').all();check(patient+'_immune_coordinate_match',ok,counts)
 mapping.append({'patient':patient,'label_sheet':m['labels'],'immune_sheet':m['immune'],'label_rows':len(z),'immune_rows':len(ims),'joined':counts,'pass':bool(ok)})
 if ok:
  q=q.drop(columns='_merge');q.to_csv(P/'tables'/f'JOINED_{patient}.tsv',sep='\t',index=False);merged[patient]=q
  for sensitivity in ['all_nonmalignant','exclude_unknown']:
   a=q[q.malignant];b=q[~q.malignant]
   if sensitivity=='exclude_unknown':b=b[~b.DefineTypes.str.lower().isin(['unknown','unkown'])]
   for gene in ['CD3D','CD8A','CD4','CD19','MS4A1','GZMA','PRF1','GNLY']:
    immune.append({'patient':patient,'comparison':sensitivity,'gene':gene,'n_malignant_spots':len(a),'n_other_spots':len(b),'mean_malignant':a[gene].mean(),'mean_other':b[gene].mean(),'delta_mean':a[gene].mean()-b[gene].mean(),'detection_malignant':a[gene].gt(0).mean(),'detection_other':b[gene].gt(0).mean(),'delta_detection':a[gene].gt(0).mean()-b[gene].gt(0).mean(),'biological_P':None})
 else:
  z.to_csv(P/'tables'/f'JOINED_{patient}_scores_only.tsv',sep='\t',index=False)
# Audit exact author plot source subset identity; no outcome-driven selection.
for pat,score_sheet,label_sheet,supp in [('P6','Figure 6j','Figure 6b','Supplementary Figure 9h'),('P7','Figure 6k','Figure 6e','Supplementary Figure 9g')]:
 a=table(score_sheet);allscore=table('Supplementary Figure 9k') if pat=='P6' else table('Figure 6g').merge(table('Figure 6h'),on=key)
 vals=allscore[[c for c in allscore if c not in key]].to_numpy();av=a.to_numpy();lab=table(label_sheet)
 intended=vals[lab.DefineTypes.eq('MARCO_FAP').to_numpy()]
 if pat=='P6':check('P6_main_source_allspots',av.shape==vals.shape and np.array_equal(av,vals),{'main_rows':len(a),'caption_joint_rows':len(intended)})
 else:check('P7_main_source_joint_subset',av.shape==intended.shape and np.array_equal(av,intended),{'main_rows':len(a),'caption_joint_rows':len(intended)})
 sb=table(supp);sv=sb[[c for c in sb if c.startswith('FAP') or c.startswith('SPP1')]].copy();sv=sv[[c for c in sv if c.startswith('FAP')]+[c for c in sv if c.startswith('SPP1')]].to_numpy()
 check(pat+'_MARCO_SPP1_vector_agreement',sv.shape==vals.shape and np.max(abs(sv-vals))<=5e-5,{'max_abs_difference':float(np.max(abs(sv-vals)))})
gate();go=table('Figure 6l');grows=[]
for _,a in go.iterrows():
 k,n=map(int,a.GeneRatio.split('/'));K,N=map(int,a.BgRatio.split('/'))
 hp=float(hypergeom.sf(k-1,N,K,n));fp=float(fisher_exact([[k,n-k],[K-k,N-K-n+k]],alternative='greater').pvalue)
 err=abs(np.log10(hp)-np.log10(float(a.pvalue)))
 grows.append({'patient':a.Patient,'GO':a.ID,'term':a.Description,'k':k,'n':n,'K':K,'N':N,'fold_enrichment':(k/n)/(K/N),'author_P':a.pvalue,'hypergeom_P':hp,'fisher_P':fp,'author_FDR':a['p.adjust'],'absolute_log10P_error':err,'selected_primary':a.ID in protocol['essential_R2']['terms']})
 check(a.Patient+'_'+a.ID+'_rawP',err<=1e-6,err)
 check(a.Patient+'_'+a.ID+'_gene_count',len(a.geneID.split('/'))==k and k==a.Count,k)
c=pd.DataFrame(corr);i=pd.DataFrame(immune);g=pd.DataFrame(grows)
c.to_csv(P/'tables/REPRODUCED_CORRELATIONS.tsv',sep='\t',index=False);i.to_csv(P/'tables/IMMUNE_MAP_QUANTIFICATION.tsv',sep='\t',index=False);g.to_csv(P/'tables/GO_STATISTICS_REPRODUCTION.tsv',sep='\t',index=False)
pd.DataFrame(checks).to_csv(P/'qc/REPRODUCTION_CHECKS.tsv',sep='\t',index=False);dump('SPATIAL_INPUT_MAPPING.json',mapping)
primary=c[c.scope.eq('author_joint_cluster')]
r1=bool(primary.pearson_r.gt(0).all() and all(abs(a.pearson_r-a.paper_rounded_R)<=.005 and a.paper_spot_model_P<2.2e-16 for a in primary[primary.patient.isin(['P6','P7'])].itertuples()))
gp=g[g.selected_primary];r2=bool(len(gp)==12 and gp.hypergeom_P.lt(.05).all() and gp.fold_enrichment.gt(1).all() and g.absolute_log10P_error.le(1e-6).all())
ip=i[i.comparison.eq('all_nonmalignant')&i.gene.eq('CD3D')];r3=bool(len(ip)==4 and ip.delta_mean.lt(0).all() and ip.delta_detection.lt(0).all())
result={'time':datetime.datetime.now().astimezone().isoformat(),'R1_pass':r1,'R2_pass':r2,'R3_pass':r3,'immune_patients_available':ip.patient.tolist(),'all_input_checks_pass':all(c['pass'] for c in checks),'success':bool(r1 and r2 and r3 and all(c['pass'] for c in checks)),'scope':'source-level reproduction','external_validation_started':False,'TYMP_positioning_started':False,'failed_checks':[a for a in checks if not a['pass']]}
dump('REPRODUCTION_DECISION.json',result)
print('DECISION',json.dumps(result,indent=2),flush=True)
print('CORRELATIONS',c.to_json(orient='records'),flush=True)
print('PRIMARY_IMMUNE',ip.to_json(orient='records'),flush=True)
print('GO_MAX_LOG_ERROR',g.absolute_log10P_error.max(),flush=True)
