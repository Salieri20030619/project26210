from pathlib import Path
import zipfile,xml.etree.ElementTree as E,json,math,statistics,hashlib,subprocess
import pandas as pd,numpy as np
from scipy import stats,optimize
P=Path(__file__).resolve().parents[1];checks=[]
def check(name,ok,detail=None):
 checks.append(dict(check=name,pass_=bool(ok),detail=detail))
 if not ok:raise AssertionError((name,detail))
q=subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True);(P/'logs/validation_memory_gate.log').write_text(q.stdout+q.stderr);check('memory_gate',q.returncode==0)
ns={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'};rns='{http://schemas.openxmlformats.org/officeDocument/2006/relationships}'
raw={}
for idx in [9,10,20]:
 f=P/'downloads'/f'41588_2025_2380_MOESM{idx}_ESM.xlsx'
 with zipfile.ZipFile(f) as z:
  check(f'xlsx{idx}_crc',z.testzip() is None)
  rel={a.attrib['Id']:a.attrib['Target'] for a in E.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
  for sh in E.fromstring(z.read('xl/workbook.xml')).find('m:sheets',ns):
   target=rel[sh.attrib[rns+'id']];target=target.lstrip('/') if target.startswith('/') else 'xl/'+target
   for c in E.fromstring(z.read(target)).findall('.//m:c',ns):
    v=c.find('m:v',ns)
    if c.get('t') not in [None,'n'] or v is None or not v.text:continue
    raw[(idx,sh.attrib['name'],c.attrib['r'])]=float(v.text)
t=pd.read_csv(P/'tables/ALL_NUMERIC_SOURCE_CELLS.tsv',sep='\t');check('all_source_cells_count',len(t)==len(raw),len(raw))
for row in t.itertuples(index=False):check('source_'+str(row.workbook)+'_'+row.sheet+'_'+row.cell,math.isclose(row.value,raw[(row.workbook,row.sheet,row.cell)],rel_tol=1e-12,abs_tol=1e-12))
u=pd.read_csv(P/'tables/BIOLOGICAL_OR_DECLARED_UNIT_SUMMARY.tsv',sep='\t');s=pd.read_csv(P/'tables/ALL_GROUP_SUMMARIES.tsv',sep='\t');obs=pd.read_csv(P/'tables/ALL_OBSERVATIONS.tsv',sep='\t');keys=['panel','endpoint','group','unit','unit_kind']
for key,g in obs.groupby(keys):
 row=u
 for k,v in zip(keys,key):row=row[row[k]==v]
 check('unit_unique_'+str(key),len(row)==1)
 check('unit_mean_'+str(key),math.isclose(row.iloc[0]['mean'],statistics.mean(g.value),rel_tol=1e-12,abs_tol=1e-10))
for row in s.itertuples():
 a=u[(u.panel==row.panel)&(u.endpoint==row.endpoint)&(u.group==row.group)&(u.unit_kind==row.unit_kind)]['mean'].tolist()
 check('group_'+str((row.panel,row.endpoint,row.group)),len(a)==row.n_units and math.isclose(statistics.mean(a),row.mean,rel_tol=1e-12,abs_tol=1e-10))
models=pd.read_csv(P/'tables/QUALIFIED_MODEL_CHECKS.tsv',sep='\t');caf=models[models.panel=='Ext Fig10H'].iloc[0];wt=[raw[(20,'Ext Fig10H',f'A{i}')] for i in range(2,8)];ko=[raw[(20,'Ext Fig10H',f'B{i}')] for i in range(2,9)]
test=stats.ttest_ind(ko,wt,equal_var=True);check('CAF_OLS_P_independent_scipy',abs(test.pvalue-caf.P)<1e-12);check('CAF_delta_independent_source',abs(np.mean(ko)-np.mean(wt)-caf.estimate_KO_minus_WT)<1e-8)
# Independently fit a two-group random-intercept REML model via scalar profile likelihood.
groups=[]
for col in 'ABCDEFGHIJ':
 v=[val for (bi,sh,cell),val in raw.items() if bi==20 and sh=='Ext Fig10C' and cell.rstrip('0123456789')==col]
 groups.append(np.asarray(v))
n=np.array([len(a) for a in groups]);sy=np.array([a.sum() for a in groups]);sy2=np.array([(a*a).sum() for a in groups]);gg=np.array([0]*5+[1]*5);N=n.sum();k=2
check('CD8_animal_count',len(groups)==10 and all(n>0),{'animals':10,'observations':int(N)})
def profile(logr,details=False):
 rat=np.exp(logr);weights=n/(1+n*rat);A=np.array([weights[gg==g].sum() for g in [0,1]]);beta=np.array([(sy/(1+n*rat))[gg==g].sum()/A[g] for g in [0,1]]);b=beta[gg]
 rss=np.sum(sy2-2*b*sy+n*b*b-rat/(1+n*rat)*(sy-n*b)**2);sig=rss/(N-k)
 obj=(N-k)*np.log(sig)+np.log(1+n*rat).sum()+np.log(A).sum()
 return (beta,sig,rat*sig,A) if details else obj
fit=optimize.minimize_scalar(profile,bounds=(-15,15),method='bounded',options={'xatol':1e-10});beta,sigma2,tau2,A=profile(fit.x,True);effect=beta[1]-beta[0];se=np.sqrt(sigma2*(1/A).sum());pval=2*stats.t.sf(abs(effect/se),8);m=models[models.panel=='Ext Fig10C'].iloc[0]
check('CD8_REML_optimization',fit.success);check('CD8_REML_effect',abs(effect-m.estimate_KO_minus_WT)<1e-5);check('CD8_REML_SE',abs(se-m.SE)<1e-5);check('CD8_REML_P',abs(pval-m.P)<1e-8)
check('two_model_BH',np.allclose(stats.false_discovery_control(models.P),models.FDR_two_model_audit,atol=1e-12))
# Exact original workbook remains frozen by downloaded SHA256.
for a in json.loads((P/'resources/DOWNLOAD_MANIFEST.json').read_text()):
 if a['status']=='OK':check('download_hash_'+a['file'],hashlib.sha256((P/'downloads'/a['file']).read_bytes()).hexdigest()==a['sha256'])
# Fixed author Box-Cox transform retained as descriptive only for unidentified microscopy columns.
d=obs[obs.panel.isin(['Fig8J','Fig8K'])].copy();d['boxcox_lambda025']=4*(d.value**.25-1);v=d.groupby(keys)['boxcox_lambda025'].mean().reset_index();v.to_csv(P/'tables/COLLAGEN_ASMA_AUTHOR_TRANSFORM_DESCRIPTIVE.tsv',sep='\t',index=False)
(P/'qc/INDEPENDENT_MODEL_VERIFICATION.json').write_text(json.dumps({'CAF_P':test.pvalue,'CD8_effect':effect,'CD8_SE':se,'CD8_P':pval,'CD8_residual_variance':sigma2,'CD8_random_mouse_variance':tau2},indent=2))
(P/'qc/FINAL_VALIDATION.json').write_text(json.dumps({'status':'PASS','checks':len(checks),'numeric_source_cells':len(raw),'units':len(u),'groups':len(s),'checks_detail':checks},indent=2))
print('VALIDATION_PASS',len(checks),'SOURCE_CELLS',len(raw));print('INDEPENDENT_MODELS',test.pvalue,effect,se,pval)
