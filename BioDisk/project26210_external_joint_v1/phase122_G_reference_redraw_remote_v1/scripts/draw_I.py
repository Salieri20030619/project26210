from pathlib import Path
import io,sys,json,zipfile
import numpy as np,pandas as pd,anndata as ad,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
r=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1')
a=ad.read_h5ad(r/'02_integration_annotation/cell_annotation_v1/part1_scvi_integrated_cell_annotated_v1.h5ad',backed='r')
c=pd.read_csv(r/'03_malignancy_cnv/malignancy_cnv_summary_v1/tumor_epithelial_malignancy_calls_v1.tsv.gz',sep='\t')
c=c[c.malignancy_status_v1.isin(['malignant_consensus','malignant_infercnv_only','non_malignant_consensus'])].copy()
idx=a.obs_names.get_indexer(c.exact_cell_id);assert (idx>=0).all();o=np.argsort(idx);idx=idx[o];c=c.iloc[o].reset_index(drop=True)
genes=['NDUFA4L2','CA9','MET'];x=a[idx,genes].X.toarray();assert len(c)==12448 and np.isfinite(x).all()
c['plot_group']=np.where(c.malignancy_status_v1=='non_malignant_consensus','Non-malignant','Malignant')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'axes.linewidth':.8})
buf=io.BytesIO();z=zipfile.ZipFile(buf,'w',zipfile.ZIP_DEFLATED);rng=np.random.default_rng(26210);counts={};pcounts={}
for kind in ['violin_points']:
 fig,axs=plt.subplots(1,3,figsize=(9.6,4.6),sharey=True)
 for j,(g,ax) in enumerate(zip(genes,axs)):
  for k,(group,col) in enumerate(zip(['Non-malignant','Malignant'],['#8BBBD0','#73A56B'])):
   sel=c.plot_group==group;v=x[sel,j];counts[group]=len(v)
   if kind=='violin_points' and np.ptp(v)>1e-10:
    vi=ax.violinplot([v],positions=[k],widths=.8,showextrema=False,bw_method=.25);bb=vi['bodies'][0];bb.set_facecolor(col);bb.set_edgecolor('#333333');bb.set_alpha(.85);bb.set_linewidth(1)
   ax.boxplot([v],positions=[k],widths=.16 if kind=='violin_points' else .48,showfliers=False,patch_artist=True,medianprops={'color':'#444444','linewidth':1.3},boxprops={'facecolor':col,'edgecolor':col,'alpha':.55},whiskerprops={'color':col,'linewidth':1},capprops={'color':col,'linewidth':1})
   fr=c.loc[sel,['study','curated_patient_id']].copy();fr['expression']=v;pv=fr.groupby(['study','curated_patient_id']).expression.median().to_numpy();pcounts[group]=len(pv)
   ax.scatter(k+rng.uniform(-.15,.15,len(pv)),pv,s=22,c=col,edgecolors=col,alpha=.8,linewidths=.7,zorder=4)
  ax.set_xticks([0,1],['Non-malignant','Malignant'],rotation=35,ha='right');ax.set_title(g,fontstyle='italic');ax.spines[['top','right']].set_visible(False);ax.set_xlim(-.5,1.5);ax.set_ylim(bottom=-.15)
 axs[0].set_ylabel('Log1p expression');fig.subplots_adjust(left=.08,right=.985,top=.88,bottom=.23,wspace=.2);fig.text(.012,.94,'I',fontsize=19,fontweight='bold');fig.text(.08,.015,'Points: patient-level medians. Distributions: all cells in each CNV-defined group.',fontsize=8,color='#555555')
 for ext in ['png','pdf']:
  bb=io.BytesIO();fig.savefig(bb,format=ext,dpi=300,bbox_inches='tight',pad_inches=.12);z.writestr('I_marker_'+kind+'_v2.'+ext,bb.getvalue())
 plt.close(fig)
z.writestr('I_aggregate_validation.json',json.dumps({'status':'PASS_DATA_CHECKS','genes':genes,'cell_counts':counts,'patient_median_points':pcounts,'no_new_inferential_test':True},indent=2))
z.close();dest=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase122_G_reference_redraw_remote_v1')/'figures'; archive=zipfile.ZipFile(io.BytesIO(buf.getvalue()));
for name in archive.namelist():
 target=(dest.parent/'qc' if name.endswith('.json') else dest)/name.replace('_v2.','_v3.');target.write_bytes(archive.read(name))
print('I_RENDER_COMPLETE')
