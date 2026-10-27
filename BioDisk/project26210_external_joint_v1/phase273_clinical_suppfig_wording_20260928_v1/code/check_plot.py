from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
from scipy.stats import norm
import pymupdf as fitz
p=Path(__file__).resolve().parents[1];r=p.parent
cm=pd.read_csv(r/'phase270_SPP1_literature_clinical_revalidation_20260927_v1/tables/CHECKMATE_ANALYSIS_INPUT.tsv',sep='\t')
cut=cm.Jiang5.median();d=cm[cm.arm=='NIVOLUMAB'].copy();d['group']=np.where(d.Jiang5>=cut,'High','Low')
km=pd.read_csv(p/'data/KM_POINTS.tsv',sep='\t');risk=pd.read_csv(p/'data/KM_RISK_TABLE.tsv',sep='\t');diff=[]
for group in ['Low','High']:
 z=d[d.group==group];S=1.;green=0.
 for tm in sorted(z.OS_time.unique()):
  nr=int((z.OS_time>=tm).sum());de=int(((z.OS_time==tm)&(z.OS_event==1)).sum())
  if de:
   S*=1-de/nr
   green+=de/(nr*(nr-de))
  got=km[(km.group==group)&np.isclose(km.time,tm,rtol=0,atol=1e-9)]
  assert len(got)==1
  row=got.iloc[0];se=green**.5
  lo=S*np.exp(-norm.ppf(.975)*se);hi=min(1.,S*np.exp(norm.ppf(.975)*se))
  err=max(abs(S-row.survival),abs(lo-row.lower),abs(hi-row.upper));assert err<1e-9
  diff.append(err)
for row in risk.itertuples():assert row.n_risk==sum((d.group==row.group)&(d.OS_time>=row.time))
md=pd.read_csv(p/'data/C4_TIME_MODEL.tsv',sep='\t').iloc[0];curve=pd.read_csv(p/'data/C4_TIME_HR_CURVE.tsv',sep='\t')
np.testing.assert_allclose(curve.HR,np.exp(md.beta_0+md.beta_log1p_month*np.log1p(curve.month)),atol=1e-12)
pts=pd.read_csv(r/'phase271_C4_state_clinical_extension_20260928_v1/tables/C4_PFS_TIME_SPECIFIC_HR.tsv',sep='\t')
for row in pts[pts.model_id==12].itertuples():
 got=curve[np.isclose(curve.month,row.month)].iloc[0]
 np.testing.assert_allclose(got[['HR','CI_low','CI_high']].astype(float),[row.HR,row.CI_low,row.CI_high],atol=1e-9)
pdf=fitz.open(p/'figures/Supplementary_Figure_S2.pdf')
assert len(pdf)==1 and abs(pdf[0].rect.width-595.2756)<.1 and abs(pdf[0].rect.height-841.8898)<.1
blocks=pdf[0].get_text('blocks')
assert all(b[0]>=0 and b[1]>=0 and b[2]<=pdf[0].rect.width+.1 and b[3]<=pdf[0].rect.height+.1 for b in blocks)
doc=fitz.open(p/'report/26210_manuscript_clinical_suppfig_v1.pdf');assert '再分析' not in ''.join(x.get_text() for x in doc)
old=r/'phase272_manuscript_clinical_integration_20260928_v1/qc/render'
changed=[]
for i in range(len(doc)):
 f=p/f'qc/render/page-{i+1:02}.png';o=old/f.name
 if not o.exists() or hashlib.sha256(f.read_bytes()).hexdigest()!=hashlib.sha256(o.read_bytes()).hexdigest():changed.append(i+1)
vis=p/'qc/visual';vis.mkdir(exist_ok=True)
for pg in changed:doc[pg-1].get_pixmap(matrix=fitz.Matrix(120/72,120/72)).save(vis/f'page-{pg:02}.png')
(p/'qc/INDEPENDENT_PLOT_CHECK.json').write_text(json.dumps({'KM_product_limit_and_log_CI_verified':True,'max_KM_error':max(diff),'risk_table_verified':True,'C4_curve_formula_and_frozen_HR_points_verified':True,'figure_A4_bounds_pass':True,'reanalysis_word_absent':True,'pages':len(doc),'changed_pages':changed,'unchanged_pages_pixel_identical':len(doc)-len(changed)},indent=2))
print('CHECK_PASS',len(doc),'changed pages',changed,'KM error',max(diff))
