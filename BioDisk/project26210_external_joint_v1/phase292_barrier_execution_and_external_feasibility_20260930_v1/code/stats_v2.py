from pathlib import Path
import importlib.util,json,hashlib
P=Path(__file__).resolve().parents[1]
s=(P/'code/analyze.py').read_text()
needle=" cols=[x,y]+cov"
assert needle in s
s=s.replace(needle," cov=[c for c in cov if d[c].nunique(dropna=True)>1]\n"+needle,1)
(P/'code/analysis_parameterization_v2.py').write_text(s)
spec=importlib.util.spec_from_file_location('a2',P/'code/analysis_parameterization_v2.py');a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a)
a.gate('model_parameterization_v2')
freeze=dict(addendum_sha256=a.sha(P/'MODEL_PARAMETERIZATION_ADDENDUM_v2.md'),original_freeze=json.loads((P/'FREEZE.json').read_text()),created_before_v2_statistics=True)
a.dump('FREEZE_v2.json',freeze)
d=a.read(P/'tables/PATIENT_INPUTS.tsv');panels=json.loads((P/'PANELS_FROZEN.json').read_text())
out=[];leave=[];genes=[]
models=['study_only','study_plus_perivascular_fraction','study_plus_perivascular_fraction_plus_receiver_depth']
gd=a.read(P/'tables/PATIENT_GENE_COUNTS.tsv').pivot(index='patient',columns='gene',values='logCP10K')
for src in ['malignant_TYMP','macro_TYMP','macro_SPP1']:
 for model in models:
  for end in panels:
   z=d[d.eligible&d[end].notna()].copy()
   dropped='perivascular_fraction' if model!='study_only' else ''
   equiv='study_only' if model=='study_plus_perivascular_fraction' else ('study_plus_receiver_depth' if model.endswith('receiver_depth') else 'study_only')
   res=a.cross(z,src,end,model);out.append(dict(source=src,endpoint=end,model=model,equivalent_model=equiv,dropped_constant_covariates=dropped,family=src+'__'+model,family_size=3,**res))
   for kind,vals in [('study',sorted(z.study.unique())),('patient',z.patient)]:
    for value in vals:
     w=z[z[kind].ne(value)];r=a.cross(w,src,end,model,False)
     leave.append(dict(source=src,endpoint=end,model=model,omission_type=kind,omitted=value,**r))
   if model=='study_plus_perivascular_fraction':
    for g in panels[end]:
     gs=[k for k in panels[end] if k in gd.columns and k!=g];w=z.copy();w['minus_gene_score']=gd.loc[w.patient,gs].mean(axis=1).to_numpy()
     genes.append(dict(source=src,endpoint=end,removed_gene=g,**a.cross(w,src,'minus_gene_score',model,False)))
df=a.pd.DataFrame(out);df['FDR']=a.np.nan
for _,ix in df.groupby('family').groups.items():
 q=a.stats.false_discovery_control(df.loc[ix,'P'].fillna(1).to_numpy(),method='bh')
 df.loc[ix,'FDR']=a.np.where(df.loc[ix,'P'].notna(),q,a.np.nan)
df['evidence']=a.np.where(df.P.isna(),'not_estimable',a.np.where(df.FDR.lt(.05),'FDR_significant',a.np.where(df.P.lt(.05),'nominal','not_supported')))
a.tab(df,'ASSOCIATIONS_v2.tsv');a.tab(leave,'LEAVE_OUT_DESCRIPTIVE_v2.tsv');a.tab(genes,'LEAVE_GENE_DESCRIPTIVE_v2.tsv')
print(df[['source','endpoint','equivalent_model','rho','P','FDR']].to_string(index=False))
print('STATISTICS_V2_COMPLETE')
