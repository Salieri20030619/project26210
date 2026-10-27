from pathlib import Path
import json,subprocess
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Rectangle
w=Path(__file__).resolve().parents[1];q=w/'qc';out=w/'figures'
p=w.parent/'local_records_migration_20260908_150455/project/phase30_part1_enzyme_function_axis_local_v1'
s=pd.read_csv(q/'PRIOR_ENZYME_AXIS_COMPLETE_SUMMARY.tsv',sep='\t');d=pd.read_csv(p/'results_local_v1/enzyme_axis_patient_effects_v1.tsv',sep='\t')
keys=['feature','feature_class','endpoint','method'];checks=0
for r in s.itertuples():
 ix=np.ones(len(d),bool)
 for key in keys:ix &= d[key].values==getattr(r,key)
 vals=d.loc[ix,'effect'].dropna().values
 assert len(vals)==r.n_eligible_patients and np.isclose(np.median(vals),r.median_patient_effect,atol=1e-12)
 assert np.sum(vals>0)==r.n_positive and np.sum(vals<0)==r.n_negative and np.sum(vals==0)==r.n_zero
 checks+=1
for _,sub in s.groupby(['endpoint','method','feature_class'],dropna=False):
 pv=sub.wilcoxon_p_two_sided.values;order=np.argsort(pv);adj=np.empty(len(pv));adj[order]=np.minimum(1,np.minimum.accumulate((pv[order]*len(pv)/np.arange(1,len(pv)+1))[::-1])[::-1])
 assert np.allclose(adj,sub.BH_FDR,equal_nan=True,atol=1e-10)
axes=['TYMP_axis_primary','TYMP_axis_broad_sensitivity','RRM2B_axis_primary','RRM2B_axis_broad_sensitivity']
ends=[('F1_continuous','Spearman'),('F1_dominant_group','median_delta'),('F1_patient_quartile','median_delta')]+[('trajectory_terminal_shift',x) for x in ['Monocle3','Slingshot','scTour']]
mat=np.zeros((4,6));rows=[]
for i,axis in enumerate(axes):
 for j,(ep,method) in enumerate(ends):
  r=s[(s.feature==axis)&(s.endpoint==ep)&(s.method==method)].iloc[0];rows.append(r.to_dict());mat[i,j]=np.sign(r.median_patient_effect)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
fig,ax=plt.subplots(figsize=(17,7));fig.subplots_adjust(left=.15,right=.98,top=.81,bottom=.23)
ax.imshow(mat,cmap=ListedColormap(['#B7CDDF','#F2F3F4','#E9B6A7']),vmin=-1,vmax=1,aspect='auto')
for i,axis in enumerate(axes):
 for j,(ep,method) in enumerate(ends):
  r=s[(s.feature==axis)&(s.endpoint==ep)&(s.method==method)].iloc[0]
  ax.text(j,i,f'effect={r.median_patient_effect:+.3f}\n+ / - / 0 = {int(r.n_positive)}/{int(r.n_negative)}/{int(r.n_zero)}\nP={r.wilcoxon_p_two_sided:.3g}; q={r.BH_FDR:.3g}',ha='center',va='center',fontsize=8.5)
  if r.BH_FDR<.05:ax.add_patch(Rectangle((j-.48,i-.48),.96,.96,fill=False,ec='#363B41',lw=1.5))
  elif r.wilcoxon_p_two_sided<.05:ax.add_patch(Rectangle((j-.48,i-.48),.96,.96,fill=False,ec='#C48A2B',lw=2))
ax.set_yticks(range(4),['TYMP strict proxy','TYMP broad proxy','RRM2B strict proxy','RRM2B broad proxy'])
ax.set_xticks(range(6),['F1 correlation\nrho | n=20','F1-high clone\ndelta | n=14','F1 Q4 minus Q1\ndelta | n=20','Monocle3 terminal\ndelta | n=20','Slingshot terminal\ndelta | n=20','scTour terminal\ndelta | n=20'],fontsize=9);ax.tick_params(length=0)
fig.suptitle('Existing enzyme-function proxies contain directional and nominal evidence',x=.04,y=.964,ha='left',fontsize=16,fontweight='bold')
fig.text(.04,.91,'Recovery of the original Phase30 results under the current exploratory reporting goal; original P values and FDR unchanged.',fontsize=10,color='#52626D')
fig.text(.04,.055,'Fill indicates effect direction only: red positive, blue negative, grey zero. Black border: q<0.05; orange border: P<0.05 but q>=0.05.\nProxies combine RNA and scFEA components, not measured protein activity. + / - / 0 shows all patient directions, including ties.\nThe historical column name F1_dominant_group actually compares F1-high primary clones with other primary clones (verified in the original script).',fontsize=8.5,linespacing=1.45)
name='04_prior_enzyme_axis_evidence_recovered_v1'
for ext in ['png','pdf']:fig.savefig(out/(name+'.'+ext),dpi=360)
plt.close(fig);subprocess.run(['pdftoppm','-r','125','-singlefile','-png',str(out/(name+'.pdf')),str(q/(name+'_pdf'))],check=True)
pd.DataFrame(rows).to_csv(q/'PRIOR_AXIS_DISPLAY_SUMMARY.tsv',sep='\t',index=False)
(q/'PRIOR_AXIS_RECHECK.json').write_text(json.dumps({'status':'PASS','original_summary_effect_and_sign_rows_checked':checks,'all_original_BH_families_recomputed':True,'original_P_values_retained':True,'clone_group_name_corrected_in_display':True,'new_tests':False},indent=2))
print('PRIOR_AXIS_RESTORED',checks,flush=True)
