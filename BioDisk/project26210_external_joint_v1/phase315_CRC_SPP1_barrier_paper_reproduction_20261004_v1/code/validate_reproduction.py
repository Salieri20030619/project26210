from pathlib import Path
import openpyxl,math,json,hashlib,subprocess,csv,datetime
import numpy as np,pandas as pd
from scipy.stats import rankdata,t
P=Path(__file__).resolve().parents[1]
subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],check=True)
w=openpyxl.load_workbook(P/'downloads/source.xlsx',read_only=True,data_only=True)
def raw(name):
 a=list(w[name].values);return pd.DataFrame(a[1:],columns=a[0])
checks=[]
def ck(name,ok,detail=None):checks.append({'check':name,'pass':bool(ok),'detail':detail})
def pc(a,b):
 a=list(map(float,a));b=list(map(float,b));n=len(a);am=math.fsum(a)/n;bm=math.fsum(b)/n
 return math.fsum((x-am)*(y-bm) for x,y in zip(a,b))/math.sqrt(math.fsum((x-am)**2 for x in a)*math.fsum((y-bm)**2 for y in b))
pr=json.loads((P/'resources/REPRODUCTION_PROTOCOL_FROZEN.json').read_text());results=pd.read_csv(P/'tables/REPRODUCED_CORRELATIONS.tsv',sep='\t');imm=pd.read_csv(P/'tables/IMMUNE_MAP_QUANTIFICATION.tsv',sep='\t');coords=['imagerow','imagecol']
ck('protocol_hash',hashlib.sha256((P/'resources/REPRODUCTION_PROTOCOL_FROZEN.json').read_bytes()).hexdigest()==(P/'resources/REPRODUCTION_PROTOCOL_FROZEN.sha256').read_text())
for f,sha in pr['original_files'].items():ck('input_SHA_'+f,hashlib.sha256((P/'downloads'/f).read_bytes()).hexdigest()==sha)
for patient,m in pr['identity'].items():
 subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],check=True)
 lab=raw(m['labels'])
 if 'scores' in m:sc=raw(m['scores'])
 else:sc=raw(m['FAP']).merge(raw(m['SPP1']),on=coords,validate='one_to_one')
 cols=[c for c in sc if c not in coords];sc=sc.rename(columns={cols[0]:'FAP',cols[1]:'SPP1'})
 z=lab.merge(sc,on=coords,validate='one_to_one');ck(patient+'_exact_n',len(z)==pr['units'][patient])
 for scope,a in [('author_joint_cluster',z[z.DefineTypes==m['joint_label']]),('whole_section_diagnostic',z)]:
  row=results[(results.patient==patient)&(results.scope==scope)].iloc[0]
  r=pc(a.FAP,a.SPP1);rho=pc(rankdata(a.FAP),rankdata(a.SPP1))
  ck(patient+scope+'_pearson_manual',abs(r-row.pearson_r)<1e-12,abs(r-row.pearson_r));ck(patient+scope+'_spearman_manual',abs(rho-row.spearman_r_diagnostic)<1e-12)
  if row.paper_spot_model_P>0:
   pp=2*t.sf(abs(r)*math.sqrt((len(a)-2)/(1-r*r)),len(a)-2)
   ck(patient+scope+'_P_t_formula',abs(math.log10(pp)-math.log10(row.paper_spot_model_P))<1e-8)
 ims=raw(m['immune']);coord_overlap=len(set(map(tuple,z[coords].to_numpy()))&set(map(tuple,ims[coords].to_numpy())))
 ck(patient+'_coordinate_overlap_record',coord_overlap==(0 if patient=='P7' else len(z)),coord_overlap)
 if patient!='P7':
  z=z.merge(ims,on=coords,validate='one_to_one')
  for row in imm[imm.patient==patient].itertuples():
   a=z[z.DefineTypes.eq('Malignant epithelial cells')];b=z[~z.DefineTypes.eq('Malignant epithelial cells')]
   if row.comparison=='exclude_unknown':b=b[~b.DefineTypes.str.lower().isin(['unknown','unkown'])]
   vals={'mean_malignant':math.fsum(a[row.gene])/len(a),'mean_other':math.fsum(b[row.gene])/len(b),'detection_malignant':sum(v>0 for v in a[row.gene])/len(a),'detection_other':sum(v>0 for v in b[row.gene])/len(b)}
   for k,v in vals.items():ck(patient+row.comparison+row.gene+'_'+k,abs(v-getattr(row,k))<1e-12)
a=raw('Figure 6m');b=raw('Supplementary Figure 9s');duplicate=bool(a.equals(b));ck('P7_candidate_immune_is_exact_P6_duplicate',duplicate,{'rows':len(a),'columns':list(a.columns)})
go=raw('Figure 6l');gt=pd.read_csv(P/'tables/GO_STATISTICS_REPRODUCTION.tsv',sep='\t');maxerr=0
for row in go.itertuples(index=False):
 k,n=map(int,row.GeneRatio.split('/'));K,N=map(int,row.BgRatio.split('/'))
 def lc(n,k):return math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1)
 terms=[lc(K,j)+lc(N-K,n-j)-lc(N,n) for j in range(k,min(K,n)+1)];mx=max(terms);logp=(mx+math.log(math.fsum(math.exp(v-mx) for v in terms)))/math.log(10)
 err=abs(logp-math.log10(row.pvalue));maxerr=max(maxerr,err);ck(row.Patient+'_'+row.ID+'_loggamma_independent',err<1e-8,err)
ck('no_external_or_TYMP_stage',not json.loads((P/'qc/REPRODUCTION_DECISION.json').read_text())['external_validation_started'])
# Read-only source inventory verifies no other published immune-coordinate worksheet is silently omitted.
wslist=[]
for ws in w:
 headers=next(ws.iter_rows(min_row=1,max_row=1,values_only=True))
 if 'CD3D' in headers and 'imagerow' in headers:wslist.append({'sheet':ws.title,'rows':ws.max_row-1,'cols':ws.max_column})
w.close();(P/'qc/IMMUNE_SOURCE_INVENTORY.json').write_text(json.dumps(wslist,indent=2))
pd.DataFrame(checks).to_csv(P/'qc/INDEPENDENT_VALIDATION_CHECKS.tsv',sep='\t',index=False)
out={'status':'PASS' if all(a['pass'] for a in checks) else 'FAIL','checks':len(checks),'failed':[a for a in checks if not a['pass']],'time':datetime.datetime.now().astimezone().isoformat(),'GO_max_log10P_error_independent':maxerr,'P6_duplicate_in_P7_candidate':duplicate,'meaning':'Verification of numerical results and recorded reproduction shortfalls; PASS does not mean full scientific reproduction success','method':'Original XLSX openpyxl extraction; manual centered correlation/rank correlation/t-tail; fsum means/detection; loggamma combinatorial GO tail; exact coordinates and hashes'}
(P/'qc/FINAL_VALIDATION.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2));print('IMMUNE_SHEETS',json.dumps(wslist));print('CD8',imm[(imm.comparison=='all_nonmalignant')&(imm.gene=='CD8A')].to_json(orient='records'))
