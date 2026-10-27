from pathlib import Path
import pandas as pd,numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
r=Path.cwd();d=r/'phase202_candidate_labels_comparison_20260918_v6'
a=pd.read_csv(r/'phase200_union_C4_C10_C15_reassessment_20260918_v1/bridge/SAME_LIBRARY_GROUP_METRICS.tsv',sep='\t')
b=pd.read_csv(r/'phase154_same_library_TYMP_C4_test_20260915_v1/SAME_LIBRARY_PATIENT_METRICS.tsv',sep='\t')
assert a.patient.nunique()==len(b)==16
assert set(a.patient)==set(b.patient)
from scipy import sparse
core=pd.read_csv(r/'phase200_union_C4_C10_C15_reassessment_20260918_v1/core/PATIENT_COMPOSITION_SPP1.tsv',sep='\t');core=core[core.scope=='all']
m=pd.read_csv(r/'phase200_union_C4_C10_C15_reassessment_20260918_v1/inputs/MACRO_METADATA.tsv.gz',sep='\t',low_memory=False)
x=sparse.load_npz(r/'phase200_union_C4_C10_C15_reassessment_20260918_v1/inputs/MACRO_COUNTS.npz')
genes=pd.read_csv(r/'phase200_union_C4_C10_C15_reassessment_20260918_v1/inputs/GENES.tsv',sep='\t').gene.tolist()
assert x.shape==(len(m),len(genes))
v=x[:,genes.index('SPP1')].toarray().ravel();m['raw_spp1']=v;m['det_spp1']=v>0
whole=m.groupby('curated_patient_id').agg(raw=('raw_spp1','sum'),lib=('library_counts','sum'),detection=('det_spp1','mean'),cells=('raw_spp1','size'))
whole=whole[whole.cells>=20];whole['expression']=np.log2(1+whole.raw/whole.lib*1e6)
whole.to_csv(d/'tables/WHOLE_DISCOVERY_PATIENT_VALUES.tsv',sep='\t')
core.to_csv(d/'tables/DISCOVERY_PATIENT_VALUES.tsv',sep='\t',index=False)
rows=[]
for g in ['C4','C10','C15','U','All']:
 if g=='All':
  expr=whole.expression;det=whole.detection;share=pd.Series(np.ones(16))
 else:
  q=core[core.group==g];eligible=q[(q.cells>=20)&((q.total_macrophages-q.cells)>=20)];expr=eligible.SPP1_log2CPM;det=eligible.SPP1_detected_fraction;share=a[a.group==g].SPP1_RNA_share
 rows.append(dict(group=g,n_expression=len(expr),expression=expr.median(),detection=100*det.median(),n_contribution=len(share),contribution=100*share.median()))
s=pd.DataFrame(rows);s.to_csv(d/'tables/CANDIDATE_VS_WHOLE_SUMMARY.tsv',sep='\t',index=False);a.to_csv(d/'tables/CANDIDATE_PATIENT_VALUES.tsv',sep='\t',index=False);b.to_csv(d/'tables/WHOLE_MACROPHAGE_PATIENT_VALUES.tsv',sep='\t',index=False)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'svg.fonttype':'none'})
fig,axs=plt.subplots(1,3,figsize=(12,4.8),gridspec_kw={'wspace':.36})
colors=['#993F80','#146A98','#397348','#B9852D','#89929B'];names=['C4','C10','C15','C4 + C10 + C15','All macrophages']
for k,(ax,key,title,xlab,lim) in enumerate(zip(axs,['expression','detection','contribution'],['SPP1 expression','SPP1 detection','SPP1 RNA contribution'],['Median log2(CPM + 1)','Median positive cells (%)','Median share of total SPP1 RNA (%)'],[14,120,145])):
 vals=s[key].to_numpy();ax.barh(np.arange(5),vals,color=colors,height=.58)
 ax.set_yticks(np.arange(5));ax.set_yticklabels(names if k==0 else ['']*5);ax.invert_yaxis();ax.set_xlim(0,lim);ax.set_title(chr(65+k)+'  '+title,loc='left',fontweight='bold',pad=13);ax.set_xlabel(xlab,fontsize=8.5)
 for i,v in enumerate(vals):
  n=int(s.iloc[i]['n_contribution' if k==2 else 'n_expression']);t=f'{v:.2f}' if k==0 else f'{v:.1f}%';ax.text(v+lim*.018,i,t+f' | n={n}',va='center',fontsize=8)
 for edge in ['top','right','left']:ax.spines[edge].set_visible(False)
 ax.tick_params(axis='y',length=0);ax.grid(axis='x',color='#ededed',linewidth=.6);ax.set_axisbelow(True)
fig.suptitle('SPP1 across candidate states, their union and all macrophages',fontsize=13,fontweight='bold',y=.96)
fig.subplots_adjust(left=.135,right=.97,top=.80,bottom=.29)
fig.text(.135,.15,'A-B: full discovery cohort. C: the original 16-patient matched cohort. Each patient has equal weight.',fontsize=9)
fig.text(.135,.10,'A-B: states/union require >=20 target and >=20 other macrophages per patient; whole population requires >=20 macrophages.',fontsize=8)
fig.text(.135,.055,'C: all 16 patients, including zero contribution for absent states; denominator = all macrophage SPP1 RNA in that patient.',fontsize=8)
fig.text(.135,.012,'n = eligible patients, not cells. C10 (n=2) is descriptive only. Groups overlap; medians are not additive.',fontsize=8,color='#555555')
for ext in ['pdf','svg','png']:fig.savefig(d/'figures'/('Candidate_states_vs_whole_macrophages_v6.'+ext),dpi=600 if ext=='png' else 150)
fig.savefig(d/'qc/comparison.jpg',dpi=150)
print(s.to_string(index=False))
