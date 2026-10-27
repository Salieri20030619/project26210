from pathlib import Path
import numpy as np,pandas as pd,json,hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1')
d=r/'phase201_three_candidates_C4_main_report_20260918_v4'
h=pd.read_csv(r/'phase192_sections28_212_figures_20260917_v1/tables/MARKER_PATIENT_EQUAL_LOG2CPM.tsv',sep='\t',index_col=0)
markers=['SPP1','GPNMB','CCL18','CYP27A1','USP53','ALOX5AP','MKI67','TOP2A']
full=h.loc[markers].to_numpy()
z=(full-full.mean(1,keepdims=True))/full.std(1,keepdims=True)
cols=['C4','C10','C15'];zz=z[:,[list(h.columns).index(c) for c in cols]]
m=pd.read_csv(r/'phase192_sections28_212_figures_20260917_v1/inputs/MACRO_DISPLAY_METADATA.tsv.gz',sep='\t',index_col=0,dtype={'candidate':str})
counts=pd.crosstab(m.candidate,m.study).loc[['4','10','15']]
qc=pd.read_csv(r/'part2_rebuild_v1/05_patient_study_validation/run_v4_macrophage_only/CLUSTER_COMPOSITION_QC.tsv',sep='\t').set_index('cluster').loc[[4,10,15]]
assert counts.sum(1).tolist()==[2674,1315,513]
assert qc.patients_paired_marker_test.tolist()==[26,2,9]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
cm=LinearSegmentedColormap.from_list('marker',['#759DBF','#FAF9F5','#993F80'])
fig=plt.figure(figsize=(8.2,3.65))
ax=fig.add_axes([.12,.22,.23,.68])
im=ax.imshow(zz,cmap=cm,vmin=-2.5,vmax=2.5,aspect='auto',interpolation='nearest')
ax.set_xticks(range(3),cols,fontsize=9)
ax.set_yticks(range(8),markers,fontsize=9,fontstyle='italic')
ax.tick_params(length=0)
for y in [0.5,3.5,5.5]:ax.axhline(y,color='white',lw=1.5)
for spine in ax.spines.values():spine.set_visible(False)
ca=fig.add_axes([.37,.34,.014,.39]);cb=fig.colorbar(im,cax=ca,ticks=[-2,0,2]);cb.ax.tick_params(labelsize=7,length=2);ca.set_title('Row z',fontsize=8,pad=6)
fig.text(.04,.96,'A',fontsize=12,weight='bold')
fig.text(.085,.96,'Representative markers',fontsize=10,weight='bold')
ax2=fig.add_axes([.49,.43,.35,.42])
palette=['#6A92B8','#54A3A0','#D9B46B','#967BB6','#C88080']
frac=counts.div(counts.sum(1),axis=0)*100
left=np.zeros(3)
for st,col in zip(frac.columns,palette):
 vals=frac[st].to_numpy()
 ax2.barh(np.arange(3),vals,left=left,height=.56,color=col,edgecolor='white',linewidth=.4,label=st)
 left+=vals
ax2.set_yticks(range(3),cols,fontsize=9);ax2.invert_yaxis();ax2.set_xlim(0,100)
ax2.set_xticks([0,50,100]);ax2.set_xlabel('Cell share by study (%)',fontsize=8)
ax2.tick_params(axis='y',length=0)
ax2.set_ylim(2.55,-.65)
for i,n in enumerate([26,2,9]):ax2.text(104,i,f'n = {n}',fontsize=10,weight='bold' if i==0 else 'normal',va='center',clip_on=False,color='#993F80' if i==0 else '#333333')
ax2.text(104,-.65,'Paired\npatients',fontsize=8,va='bottom',clip_on=False)
ax2.legend(loc='upper left',bbox_to_anchor=(-.04,-.43),ncol=2,frameon=False,fontsize=7.2,handlelength=1.2,columnspacing=1.4)
fig.text(.45,.96,'B',fontsize=12,weight='bold')
fig.text(.49,.96,'Study distribution and patient support',fontsize=10,weight='bold')
fig.text(.12,.09,'Colors scaled across the original 12 states',fontsize=7,color='#555555')
fig.text(.12,.045,'Fixed cluster identities and patient-equal expression',fontsize=7,color='#555555')
for ext in ['pdf','svg','png']:
 fig.savefig(d/'figures'/f'Figure8_three_SPP1_candidates.{ext}',dpi=600)
plt.close(fig)
h.loc[markers,cols].to_csv(d/'qc/MARKER_SOURCE_VALUES.tsv',sep='\t')
pd.DataFrame(zz,index=markers,columns=cols).to_csv(d/'qc/MARKER_DISPLAY_Z.tsv',sep='\t')
counts.to_csv(d/'qc/STUDY_CELL_COUNTS.tsv',sep='\t')
qc.to_csv(d/'qc/CANDIDATE_PATIENT_SUPPORT.tsv',sep='\t')
(d/'qc/FIGURE_VALIDATION.json').write_text(json.dumps({'new_figures':1,'panels':2,'source_states':12,'displayed_states':cols,'paired_patients':[26,2,9],'cells':[2674,1315,513],'no_new_inferential_tests':True,'row_z_uses_all12_original_states':True},indent=2))
print('CANDIDATE_FIGURE_READY')
