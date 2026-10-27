from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,re,zipfile,xml.etree.ElementTree as ET,subprocess
P=Path(__file__).resolve().parents[1];R=P.parent
q=subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True)
(P/'logs/perturbation_memory.log').write_text(q.stdout+q.stderr);assert q.returncode==0
mapping=pd.read_csv(R/'phase280_TYMP_source_niche_function_20260928_v1/tables/EXTERNAL_ONE_TO_ONE_MAPPING.tsv',sep='\t').dropna()
assert mapping.ENSEMBL.is_unique and mapping.SYMBOL.is_unique
mp=mapping.set_index('ENSEMBL').SYMBOL.to_dict()
meta=json.loads((R/'phase292_barrier_execution_and_external_feasibility_20260930_v1/resources/GSE204809_SAMPLE_DESIGN.json').read_text())
lookup={}
for row in meta:
 key=re.search(r'\[([^\]]+)\]',row['title'][0]).group(1).lower()
 d=dict(sample=row['sample'],title=row['title'][0])
 for v in row['characteristics_ch1']:
  k,val=v.split(': ',1);d[k]=val
 lookup[key]=d
panels=json.loads((R/'phase292_barrier_execution_and_external_feasibility_20260930_v1/PANELS_FROZEN.json').read_text())
targets=['TYMP','SPP1','CD44','ICAM1','VCAM1','CXCL9','CXCL10','CXCL11','CCL5']
effects=[];design=[];scores=[];coverage=[];checks=[];cache=[];fullgenes=None
def check(name,v):checks.append(dict(check=name,passed=bool(v)));assert v,name
for fn in sorted((P/'inputs').glob('GSE204809*.xlsx')):
 with zipfile.ZipFile(fn) as z:check('zip_'+fn.name,z.testzip() is None)
 df=pd.read_excel(fn,engine='openpyxl',index_col=0);df.index=df.index.astype(str)
 genes=df.index.str.replace(r'\.\d+$','',regex=True)
 check('unique_ensembl_'+fn.name,genes.is_unique)
 C=df.to_numpy(float);check('safe_values_'+fn.name,np.isfinite(C).all() and (C>=0).all())
 if fullgenes is None:fullgenes=list(genes)
 check('same_genes_'+fn.name,list(genes)==fullgenes)
 columns=[]
 for col in df.columns:
  key=re.sub('_R1_001$','',col).lower();check('sample_match_'+col,key in lookup)
  r=lookup[key];t=r['treatment'];rec=dict(file=fn.name,column=col,GSM=r['sample'],cell_line=r['cell line'],treatment=t,RIG=int(t.startswith('RIG')),TPI=int('Tipiracil' in t),library_date_label=re.search(r'RNA(\d+)SA',col).group(1),donor='NOT_DECLARED')
  design.append(rec);columns.append(rec)
 rig=np.array([r['RIG'] for r in columns],bool);tot=C.sum(0);check('nonzero_libraries_'+fn.name,(tot>0).all())
 uq=np.array([np.quantile(v[v>0],.75) for v in C.T])
 # UQ scaled to geometric mean total library size for comparable log CPM magnitude.
 uq_factors=uq/np.exp(np.log(uq).mean())*np.exp(np.log(tot).mean())
 logs={'CPM':np.log2(C/tot[None,:]*1e6+1),'upper_quartile':np.log2(C/uq_factors[None,:]*1e6+1)}
 for normalization,L in logs.items():
  delta=L[:,rig].mean(1)-L[:,~rig].mean(1)
  for j,g in enumerate(genes):
   effects.append(dict(file=fn.name,cell_line=columns[0]['cell_line'],TPI=columns[0]['TPI'],normalization=normalization,ENSEMBL=g,SYMBOL=mp.get(g,''),n_control=int((~rig).sum()),n_RIG=int(rig.sum()),control_mean=L[j,~rig].mean(),RIG_mean=L[j,rig].mean(),effect=delta[j],P=np.nan,FDR=np.nan))
 for program,gs in dict(panels,source_and_entry=targets).items():
  ids=[j for j,g in enumerate(genes) if mp.get(g) in gs]
  coverage.append(dict(file=fn.name,program=program,present=len(ids),requested=len(gs),missing=';'.join(sorted(set(gs)-set(mp.get(genes[j]) for j in ids)))))
  if program=='source_and_entry':continue
  for i,col in enumerate(df.columns):
   val=np.log1p(C[ids,i]/tot[i]*1e4).mean() if len(ids)>=.8*len(gs) else np.nan
   scores.append(dict(file=fn.name,column=col,GSM=columns[i]['GSM'],cell_line=columns[i]['cell_line'],TPI=columns[i]['TPI'],RIG=columns[i]['RIG'],program=program,score=val))
 # Independent raw sheet XML value replay, not pandas Excel parser.
 with zipfile.ZipFile(fn) as z:
  xml=ET.fromstring(z.read('xl/worksheets/sheet1.xml'));ns={'x':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
  B=np.zeros_like(C)
  for row in xml.findall('x:sheetData/x:row',ns):
   ri=int(row.attrib['r'])-2
   if ri<0:continue
   for c in row.findall('x:c',ns):
    ref=re.match(r'([A-Z]+)',c.attrib['r']).group(1);ci=0
    for ch in ref:ci=ci*26+ord(ch)-64
    ci-=2
    if ci<0:continue
    node=c.find('x:v',ns)
    if node is not None and c.attrib.get('t','n')=='n':B[ri,ci]=float(node.text)
  check('all_XML_values_'+fn.name,np.array_equal(C,B))
 np.savez_compressed(P/'inputs'/(fn.stem+'.npz'),counts=C,ensembl=np.array(genes),columns=np.array(df.columns))
 cache.append(dict(path=str(fn),sha256=hashlib.sha256(fn.read_bytes()).hexdigest(),rows=len(df),columns=len(df.columns),fractional_values=int(np.count_nonzero(C!=np.floor(C))),min=float(C.min()),max=float(C.max())))
e=pd.DataFrame(effects);e.to_csv(P/'tables/PERTURBATION_ALL_GENES.tsv.gz',sep='\t',index=False)
t=e[e.SYMBOL.isin(targets)];t.to_csv(P/'tables/PERTURBATION_TARGETS.tsv',sep='\t',index=False)
pd.DataFrame(scores).to_csv(P/'tables/PERTURBATION_PROGRAM_SCORES.tsv',sep='\t',index=False)
pd.DataFrame(coverage).to_csv(P/'tables/PERTURBATION_COVERAGE.tsv',sep='\t',index=False)
dd=pd.DataFrame(design);dd.to_csv(P/'tables/PERTURBATION_SAMPLE_MAP.tsv',sep='\t',index=False)
h=dd[dd.cell_line.eq('HUVEC')];D=np.column_stack([np.ones(len(h)),h.RIG,h.TPI,h.library_date_label.eq('210902').astype(int),h.RIG*h.TPI])
check('TPI_date_alias',np.array_equal(h.TPI,h.library_date_label.eq('210902').astype(int)))
check('rank_4_of_5',np.linalg.matrix_rank(D)==4 and D.shape[1]==5)
# Descriptive response-difference only; no direct drug effect.
pv=e[(e.cell_line=='HUVEC')].pivot(index=['normalization','ENSEMBL','SYMBOL'],columns='TPI',values='effect')
pv['RIG_response_difference_TPI_minus_noTPI']=pv[1]-pv[0];pv.to_csv(P/'tables/HUVEC_CROSS_BATCH_RESPONSE_DIFFERENCE.tsv.gz',sep='\t')
check('no_inferential_P',e.P.isna().all() and e.FDR.isna().all())
pd.DataFrame(checks).to_csv(P/'qc/PERTURBATION_CHECKS.tsv',sep='\t',index=False)
(P/'qc/PERTURBATION_VALIDATION.json').write_text(json.dumps(dict(status='PASS',checks=len(checks),inputs=cache,design_rank=4,design_columns=5,direct_TPI_effect='not_identifiable_separately_from_date_label',replicate_independence='unconfirmed',P_values_computed=False),indent=2))
print(t[(t.normalization=='CPM')&t.SYMBOL.isin(['TYMP','SPP1','ICAM1','VCAM1'])][['cell_line','TPI','SYMBOL','effect']].to_string(index=False))
print('PERTURBATION_CONTEXT_COMPLETE')
