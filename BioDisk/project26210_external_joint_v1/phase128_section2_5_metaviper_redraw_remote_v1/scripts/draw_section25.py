from pathlib import Path
import json,hashlib,subprocess,sys,shutil
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,Normalize
from matplotlib.ticker import MaxNLocator
from scipy.stats import rankdata
w=Path(__file__).resolve().parents[1];q=w/'qc';out=w/'figures';root=w.parent
src=Path('/mnt/BioDisk/wangmingyuan/project26210_part1_metaviper_downstream_phase29_20260814_v1')
base=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1')
paths=[]
def read(p,**kwargs):
 p=Path(p);paths.append(p);return pd.read_csv(p,**kwargs)
def tab(p):return read(p,sep='\t')
def hashfile(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
genes=['CASP9','CHMP1A','PRSS8','PXN'];allgenes=genes+['HIF1A']
methods=['Monocle3','Slingshot','scTour'];tcols={s:s+'_pseudotime_patient_v1' for s in methods}
features=allgenes+['TYMP_RNA','RRM2B_RNA']
a=tab(w/'data/SELECTED_METAVIPER_ACTIVITY.tsv').set_index('exact_cell_id')
full=tab(base/'07_pyscenic_v1/input/cell_metadata_v1.tsv.gz').set_index('exact_cell_id')
m=full.loc[a.index].copy()
traj=tab(base/'09_metabolism_trajectory/trajectory_v1/source_data/trajectory_all_methods_cell_scores_v1.tsv.gz').set_index('exact_cell_id').loc[a.index]
rna=tab(src/'m13_expression_v1/m13_tymp_rrm2b_cell_expression_v1.tsv.gz').set_index('exact_cell_id').loc[a.index]
vals=a[allgenes].copy();vals['TYMP_RNA']=rna.TYMP_log1p_cp10k;vals['RRM2B_RNA']=rna.RRM2B_log1p_cp10k
assert len(a)==11921 and a.index.is_unique and len(full)==11950
assert m.curated_patient_id.nunique()==20 and np.array_equal(m.curated_patient_id,a.curated_patient_id)
assert np.array_equal(traj.curated_patient_id,m.curated_patient_id) and np.isfinite(vals).all().all()
old=tab(root/'local_records_migration_20260908_150455/project/phase44_part1_client_figures_final_local_v1/source_remote_v1/metaviper_selected_activity_v1.tsv').set_index('exact_cell_id')
assert set(old.index)==set(a.index)
max_export_error=max(float(np.max(np.abs(a[g]-old.loc[a.index,g+'_activity']))) for g in allgenes)
assert max_export_error<1e-11
points=tab(src/'m11_analysis_v1/m11_mr_effect_activity_points_v1.tsv').set_index('protein').loc[genes]
available=tab(src/'m11_analysis_v1/m11_final_mr_target_inference_status_v1.tsv')
mrpat=tab(src/'m11_analysis_v1/m11_patient_mean_activity_v1.tsv')
mrpc=tab(src/'m11_analysis_v1/m11_patient_cross_protein_correlations_v1.tsv')
mrsummary=tab(src/'m11_analysis_v1/m11_patient_correlation_summary_v1.tsv')
mrglobal=tab(src/'m11_analysis_v1/m11_global_correlation_v1.tsv').iloc[0]
assert np.allclose(np.log(available.odds_ratio),available.mr_beta_log_or,atol=1e-13)
assert set(available.loc[available.activity_row_present,'protein'])==set(genes)
pm=a[genes].groupby(m.curated_patient_id).mean()
for r in mrpat.itertuples():assert abs(pm.loc[r.curated_patient_id,r.protein]-r.mean_activity)<1e-11
assert np.allclose(pm.mean(),points.patient_equal_mean_activity,atol=1e-11)
rglobal=float(np.corrcoef(points.mr_beta_log_or,points.patient_equal_mean_activity)[0,1])
assert abs(rglobal-mrglobal.pearson_r)<1e-12
for r in mrpc.itertuples():
 rr=np.corrcoef(points.mr_beta_log_or,pm.loc[r.curated_patient_id,genes])[0,1]
 assert abs(rr-r.pearson_r)<1e-12
effects=tab(src/'m12_analysis_v1/m12_patient_activity_effects_v1.tsv')
summ=tab(src/'m12_analysis_v1/m12_f1_activity_patient_summary_v1.tsv')
tc=tab(src/'m13_analysis_v1/m13_patient_trajectory_correlations_v1.tsv')
ts=tab(src/'m13_analysis_v1/m13_trajectory_correlation_summary_v1.tsv')
td=tab(src/'m13_analysis_v1/m13_patient_terminal_shifts_v1.tsv')
terminal=tab(src/'m13_analysis_v1/m13_terminal_shift_summary_v1.tsv')
bins=tab(src/'m13_analysis_v1/m13_patient_binned_curve_source_v1.tsv.gz')
bs=tab(src/'m13_analysis_v1/m13_patient_balanced_curve_summary_v1.tsv')
pairs=tab(src/'m13_analysis_v1/m13_pairwise_terminal_coactivation_v1.tsv')
endpoints=['continuous_F1_spearman','F1_dominant_median_difference','F1_Q4_minus_Q1_median_difference']
re=[];groupmed=[];group_ids=[]
for pid,sub in m.groupby('curated_patient_id',sort=True):
 x=sub.F1.to_numpy();labels=sub.dominant_factor_rank7_v1=='F1'
 q25,q75=np.quantile(x,[.25,.75]);lo=x<=q25;hi=x>=q75
 ce=len(sub)>=20 and np.ptp(x)>0;ge=labels.sum()>=10 and (~labels).sum()>=10;qe=q75>q25 and lo.sum()>=10 and hi.sum()>=10
 if ge:group_ids.append(pid)
 for g in allgenes:
  v=a.loc[sub.index,g].values
  corr=np.corrcoef(rankdata(x),rankdata(v))[0,1] if ce and np.ptp(v)>0 else np.nan
  gd=float(np.median(v[labels])-np.median(v[~labels])) if ge else np.nan
  qd=float(np.median(v[hi])-np.median(v[lo])) if qe else np.nan
  for ep,val in zip(endpoints,[corr,gd,qd]):
   r=effects[(effects.curated_patient_id==pid)&(effects.protein==g)&(effects.endpoint==ep)].iloc[0]
   assert np.isclose(val,r.effect,equal_nan=True,atol=1e-11)
   re.append({'patient':pid,'protein':g,'endpoint':ep,'effect':val})
  if ge:groupmed.append({'patient':pid,'protein':g,'F1':float(np.median(v[labels])),'Other':float(np.median(v[~labels])),'F1_cells':int(labels.sum()),'Other_cells':int((~labels).sum())})
assert len(group_ids)==9
pd.DataFrame(re).to_csv(q/'INDEPENDENT_F1_EFFECTS.tsv',sep='\t',index=False)
groupmed=pd.DataFrame(groupmed);groupmed.to_csv(q/'PAIRED_PATIENT_MEDIANS.tsv',sep='\t',index=False)
# Verify trajectory correlations, terminal shifts, raw binned means and support fractions.
methodreg=tab(src/'m13_analysis_v1/m13_method_registry_v1.tsv')
trajectory_rows=0;terminal_rows=0;bin_error=0
for pid,sub in m.groupby('curated_patient_id',sort=True):
 for rr in methodreg.itertuples():
  t=traj.loc[sub.index,rr.column].values
  for g in features:
   v=vals.loc[sub.index,g].values
   rho=np.corrcoef(rankdata(t),rankdata(v))[0,1] if len(sub)>=20 and np.ptp(t)>0 and np.ptp(v)>0 else np.nan
   row=tc[(tc.curated_patient_id==pid)&(tc.feature==g)&(tc.method==rr.method)].iloc[0]
   assert np.isclose(rho,row.effect,equal_nan=True,atol=1e-11);trajectory_rows+=1
   if rr.method not in methods:continue
   q20,q80=np.quantile(t,[.2,.8]);early=t<=q20;late=t>=q80
   delta=float(np.median(v[late])-np.median(v[early])) if q80>q20 and early.sum()>=5 and late.sum()>=5 else np.nan
   row=td[(td.curated_patient_id==pid)&(td.feature==g)&(td.method==rr.method)].iloc[0]
   assert np.isclose(delta,row.effect,equal_nan=True,atol=1e-11);terminal_rows+=1
   ranks=rankdata(t)/len(t);ib=np.clip(np.ceil(ranks*20),1,20).astype(int)
   actual=pd.DataFrame({'bin':ib,'v':v}).groupby('bin').agg(n=('v','size'),mu=('v','mean'))
   orig=bins[(bins.curated_patient_id==pid)&(bins.feature==g)&(bins.method==rr.method)].set_index('bin')
   assert set(actual.index)==set(orig.index)
   assert np.array_equal(actual.n,orig.loc[actual.index,'n_cells'])
   err=float(np.max(np.abs(actual.mu-orig.loc[actual.index,'mean_score'])));bin_error=max(bin_error,err)
assert bin_error<1e-11
for r in pairs.itertuples():
 aa=td[(td.method==r.method)&(td.feature==r.feature_a)].set_index('curated_patient_id').effect
 bb=td[(td.method==r.method)&(td.feature==r.feature_b)].set_index('curated_patient_id').effect
 assert int(((aa>0)&(bb>0)).sum())==r.n_both_positive
def bh(p):
 p=np.asarray(p,float);ok=np.isfinite(p);v=p[ok];ord=np.argsort(v);s=np.minimum.accumulate((v[ord]*len(v)/np.arange(1,len(v)+1))[::-1])[::-1]
 z=np.empty(len(v));z[ord]=np.minimum(s,1);out=np.full(len(p),np.nan);out[ok]=z;return out
for frame,keys in [(summ,['evidence_family','endpoint']),(ts,['evidence_family','method']),(terminal,['evidence_family','method'])]:
 for _,ss in frame.groupby(keys):
  assert np.allclose(bh(ss.wilcoxon_p_two_sided),ss.BH_FDR,equal_nan=True,atol=1e-12)
assert np.allclose(bh(mrsummary.wilcoxon_p_two_sided),mrsummary.BH_FDR,atol=1e-12)
for row in bs.itertuples():
 v=bins[(bins.method==row.method)&(bins.feature==row.feature)&(bins.bin==row.bin)].mean_score
 assert len(v)==row.n_patients and abs(v.mean()-row.patient_equal_mean)<1e-11
# Complete source statistics remain in server QC, including sensitivity methods and unavailable targets.
for folder in ['m11_analysis_v1','m12_analysis_v1','m13_analysis_v1']:
 for f in (src/folder).glob('*.tsv'):shutil.copy2(f,q/f.name)
print('ALL_SOURCE_EFFECT_AND_BIN_VALIDATION_PASS',flush=True)
# Fixed patient-local kernel smoother of actual bin means, weighted by observed bin cell counts.
# Bandwidth 0.075 patient pseudotime-rank units (=1.5 of the 20 bins), fixed for all genes/methods.
grid=np.linspace(.025,.975,120);bw=.075;fits={};fitcurves={};fitrows=[]
for method in methods:
 for g in features:
  prediction=pd.Series(index=m.index,dtype=float);curves=[]
  for pid,sub in m.groupby('curated_patient_id',sort=True):
   b=bins[(bins.method==method)&(bins.feature==g)&(bins.curated_patient_id==pid)].sort_values('bin')
   centers=b.bin_midpoint.values;observed=b.mean_score.values;weights=b.n_cells.values
   def smooth(xx):
    ww=np.exp(-.5*((np.asarray(xx)[:,None]-centers)/bw)**2)*weights
    return (ww@observed)/ww.sum(axis=1)
   curve=smooth(grid);ranks=rankdata(traj.loc[sub.index,tcols[method]])/len(sub)
   vv=smooth(ranks);assert vv.min()>=observed.min()-1e-10 and vv.max()<=observed.max()+1e-10
   prediction.loc[sub.index]=vv;curves.append(curve)
  fits[(method,g)]=prediction.to_numpy()
  mat=np.array(curves);fitcurves[(method,g)]=mat
  for j,xx in enumerate(grid):fitrows.append({'method':method,'feature':g,'pseudotime_rank':xx,'n_patients':len(curves),'patient_equal_mean':mat[:,j].mean(),'q25':np.quantile(mat[:,j],.25),'q75':np.quantile(mat[:,j],.75)})
fit_table=pd.DataFrame(fitrows);fit_table.to_csv(w/'data/FITTED_PATIENT_BALANCED_CURVES.tsv.gz',sep='\t',index=False)
pd.DataFrame({method+'__'+g:v for (method,g),v in fits.items()},index=m.index).to_csv(w/'data/FITTED_CELL_PROJECTIONS.tsv.gz',sep='\t')
# Frozen coordinates; keep all cells for the neutral backdrop.
xp=root/'phase123_section2_2_reference_redraw_remote_v1/data/plot_values_server_only.npz';paths.append(xp);z=np.load(xp)
xyall=pd.DataFrame(z['xy'],index=z['cells']);xy=xyall.loc[m.index].values
assert len(xy)==11921 and len(xyall)==11950
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.labelsize':9,'axes.titlesize':11,'axes.linewidth':.65,'pdf.fonttype':42,'savefig.facecolor':'white','xtick.major.width':.65,'ytick.major.width':.65})
gc=dict(zip(genes,['#B97767','#7299B5','#83A88A','#A68CAE']))
groupc={'Other':'#89A6BC','F1':'#C28C7E'}
div=LinearSegmentedColormap.from_list('NES',['#527F9F','#F0EFEC','#B9635F'])
seq=LinearSegmentedColormap.from_list('joint',['#F2F0ED','#B3CBD2','#547F94'])
def clean(ax):
 ax.spines[['top','right']].set_visible(False);ax.tick_params(length=3)
def save(fig,name):
 for ext in ['png','pdf']:fig.savefig(out/(name+'.'+ext),dpi=360,bbox_inches='tight',pad_inches=.16)
 plt.close(fig)
def scaleaxis(ax):
 ax.set_aspect('equal');ax.axis('off')
 ax.annotate('',xy=(.22,.05),xytext=(.045,.05),xycoords='axes fraction',arrowprops={'arrowstyle':'->','lw':.75,'color':'#4A4A4A'})
 ax.annotate('',xy=(.045,.22),xytext=(.045,.05),xycoords='axes fraction',arrowprops={'arrowstyle':'->','lw':.75,'color':'#4A4A4A'})
 ax.text(.135,.0,'UMAP1',transform=ax.transAxes,ha='center',fontsize=7)
 ax.text(.0,.13,'UMAP2',transform=ax.transAxes,rotation=90,va='center',fontsize=7)
# 1: Four actual protein points; patient sensitivity is conditional on this fixed protein panel.
fig,axs=plt.subplots(1,2,figsize=(11.5,5.2));fig.subplots_adjust(left=.09,right=.97,top=.84,bottom=.22,wspace=.4)
ax=axs[0];x=points.mr_beta_log_or.values;y=points.patient_equal_mean_activity.values
for g in genes:
 r=points.loc[g];ax.scatter(r.mr_beta_log_or,r.patient_equal_mean_activity,s=74,c=gc[g],edgecolor='white',lw=.7,zorder=3)
 offset={'CASP9':(-35,12),'CHMP1A':(-43,-17),'PRSS8':(8,9),'PXN':(10,8)}[g]
 ax.annotate(g,(r.mr_beta_log_or,r.patient_equal_mean_activity),xytext=offset,textcoords='offset points',fontsize=9)
xx=np.linspace(x.min(),x.max(),100);ax.plot(xx,np.polyval(np.polyfit(x,y,1),xx),ls='--',lw=1.1,c='#9E9E9E')
ax.set(xlabel='MR effect (log odds ratio)',ylabel='Patient-equal mean metaVIPER NES',xlim=(-.25,1.08),ylim=(-1.8,.95))
ax.set_title('Four inferable final MR proteins',loc='left')
ax.text(.02,.04,'Pearson r = '+format(rglobal,'.4f')+'\nDescriptive: 4 proteins',transform=ax.transAxes,fontsize=9);clean(ax)
ax=axs[1];pr=mrpc.sort_values('curated_patient_id');v=pr.pearson_r.values
ax.scatter(v,np.linspace(-.12,.12,len(v)),s=30,c='#7299B5',lw=0,alpha=.8)
ax.scatter(np.median(v),0,s=70,c='#333333',marker='D',zorder=4)
ax.axvline(0,c='#BDBDBD',ls='--',lw=.8);ax.set(xlim=(-1,.1),ylim=(-.5,.5),yticks=[],xlabel='Within-patient correlation across the same 4 proteins')
ax.set_title('Patient-level sensitivity: 20/20 negative',loc='left')
ss=mrsummary[mrsummary.metric=='pearson_r'].iloc[0]
ax.text(.04,.84,'Median r = '+format(ss.median_patient_correlation,'.4f')+'\nP = '+format(ss.wilcoxon_p_two_sided,'.3g')+'; FDR = '+format(ss.BH_FDR,'.3g'),transform=ax.transAxes,va='top',fontsize=9,linespacing=1.6);clean(ax)
fig.suptitle('The measured MR-effect / inferred-activity relationship is not positive',fontsize=14,x=.09,ha='left',y=.98)
fig.text(.09,.075,'Activity is transcriptome-inferred, not measured protein abundance. Patient repetition does not increase the number of MR proteins.',fontsize=8)
save(fig,'01_A_MR_effect_metaVIPER_association_v1')
# 2: Violin distributions restricted to the same nine eligible patients; connected patient medians.
paired=m[m.curated_patient_id.isin(group_ids)]
fig,axs=plt.subplots(1,4,figsize=(12.3,5.1));fig.subplots_adjust(left=.075,right=.98,top=.81,bottom=.22,wspace=.4)
for ax,g in zip(axs,genes):
 arms=[a.loc[paired.index[paired.dominant_factor_rank7_v1!='F1'],g].values,a.loc[paired.index[paired.dominant_factor_rank7_v1=='F1'],g].values]
 viol=ax.violinplot(arms,positions=[0,1],widths=.7,showextrema=False)
 for body,c in zip(viol['bodies'],groupc.values()):body.set_facecolor(c);body.set_edgecolor(c);body.set_alpha(.55);body.set_linewidth(.7)
 pat=groupmed[groupmed.protein==g].sort_values('patient')
 for i,row in enumerate(pat.itertuples()):
  shift=(i-4)*.018;ax.plot([shift,1+shift],[row.Other,row.F1],color='#777777',lw=.7,alpha=.6)
  ax.scatter([shift,1+shift],[row.Other,row.F1],s=19,c=['#507690','#A66056'],edgecolor='white',lw=.4,zorder=3)
 ss=summ[(summ.protein==g)&(summ.endpoint==endpoints[1])].iloc[0]
 ax.set(xticks=[0,1],xticklabels=['Other programs','F1-dominant'],ylabel='metaVIPER NES',xlim=(-.55,1.55))
 ax.set_title(g,loc='left',color=gc[g]);clean(ax)
 ax.text(.02,1.015,'Median delta '+format(ss.median_patient_effect,'.3f')+' | '+str(ss.n_positive)+'/9 positive',transform=ax.transAxes,fontsize=8)
fig.suptitle('F1-dominant activity contrasts are exploratory in 9 paired patients',fontsize=14,x=.075,ha='left',y=.98)
fig.text(.075,.075,'Violin: cells from the paired subset only. Connected dots: patient medians. Each arm has >=10 cells per patient.',fontsize=8)
save(fig,'02_B_F1_dominant_activity_violin_v1')
# 3: Patient-level evidence for each endpoint, with separate effect scales.
fig,axs=plt.subplots(1,3,figsize=(14.0,6.2));fig.subplots_adjust(left=.075,right=.98,top=.83,bottom=.20,wspace=.40)
titles=['Continuous F1 association','F1-dominant vs other programs','Patient-relative upper vs lower quartile']
for ax,ep,title in zip(axs,endpoints,titles):
 sub=summ[(summ.evidence_family=='final_MR_inferable')&(summ.endpoint==ep)].set_index('protein').loc[genes]
 detail=effects[(effects.protein.isin(genes))&(effects.endpoint==ep)]
 lim=1 if ep==endpoints[0] else float(np.nanmax(np.abs(detail.effect)))*1.08
 for k,g in enumerate(genes):
  dd=detail[(detail.protein==g)&detail.effect.notna()].sort_values('curated_patient_id');vv=dd.effect.values;r=sub.loc[g]
  ax.scatter(vv,np.full(len(vv),k)+np.linspace(-.11,.11,len(vv)),s=15,c=gc[g],alpha=.65,lw=0)
  ax.scatter(r.median_patient_effect,k,s=42,c='#333333',marker='D',edgecolor='white',lw=.4,zorder=3)
  if r.n_eligible_patients>=10:lab='Effect '+format(r.median_patient_effect,'.3f')+' | P '+format(r.wilcoxon_p_two_sided,'.3g')+'\nFDR '+format(r.BH_FDR,'.3g')+' | n='+str(int(r.n_eligible_patients))
  else:lab='Delta '+format(r.median_patient_effect,'.3f')+' | n='+str(int(r.n_eligible_patients))+'\n'+str(int(r.n_positive))+' positive, '+str(int(r.n_negative))+' negative'
  ax.text(-lim*.95,k+.20,lab,fontsize=7.5,va='top',linespacing=1.35)
 ax.axvline(0,c='#BBBBBB',ls='--',lw=.7)
 ax.set(xlim=(-lim,lim),ylim=(3.65,-.5),yticks=range(4),yticklabels=genes,xlabel='Spearman rho' if ep==endpoints[0] else 'Patient median activity difference (NES)')
 ax.set_title(title,loc='left',fontsize=10);ax.xaxis.set_major_locator(MaxNLocator(5));clean(ax)
fig.suptitle('CASP9, CHMP1A and PRSS8 show F1-associated inferred activity',x=.075,ha='left',y=.98,fontsize=14)
fig.text(.075,.075,'Patient is the statistical unit. FDR is adjusted within the four final MR targets for each endpoint; the n=9 panel is exploratory.',fontsize=8)
save(fig,'03_C_F1_patient_endpoint_evidence_v1')
# 4: Original activity in frozen UMAP.
def umaps(values,title,filename,limits=None):
 fig,axs=plt.subplots(2,2,figsize=(8.5,8.7));fig.subplots_adjust(left=.05,right=.98,top=.91,bottom=.06,hspace=.23,wspace=.18)
 for ax,g in zip(axs.flat,genes):
  v=values[g];order=np.argsort(v,kind='stable');lo,hi=limits[g] if limits else (float(np.min(v)),float(np.max(v)))
  ax.scatter(xyall[0],xyall[1],s=2,c='#E5E5E5',lw=0,rasterized=True)
  im=ax.scatter(xy[order,0],xy[order,1],c=v[order],s=3,cmap=div,vmin=lo,vmax=hi,lw=0,rasterized=True)
  ax.set_title(g,loc='left',color=gc[g]);scaleaxis(ax)
  cb=fig.colorbar(im,ax=ax,orientation='horizontal',fraction=.033,pad=.035,shrink=.73);cb.ax.tick_params(labelsize=7,length=2)
  cb.set_label('metaVIPER NES' if limits is None else 'Patient-trajectory fitted metaVIPER NES',fontsize=7.5)
 fig.suptitle(title,x=.05,ha='left',y=.987,fontsize=13)
 save(fig,filename)
umaps({g:a[g].values for g in genes},'Inferred activity of the four MR proteins in malignant cells','04_D_MR_protein_activity_UMAP_v1')
# 5: fixed smooth patient-balanced curves and observed bin means.
def curves(gg,title,filename,rna_mode=False):
 fig,axes=plt.subplots(len(gg),3,figsize=(11.9,2.55*len(gg)+.9),squeeze=False)
 fig.subplots_adjust(left=.09,right=.985,top=.90,bottom=.11,hspace=.42,wspace=.31)
 for ri,g in enumerate(gg):
  for ci,method in enumerate(methods):
   ax=axes[ri,ci];s=fit_table[(fit_table.feature==g)&(fit_table.method==method)]
   colour=gc.get(g,'#AE7764' if g=='TYMP_RNA' else '#7597AE')
   ax.fill_between(s.pseudotime_rank,s.q25,s.q75,color=colour,alpha=.20,lw=0)
   ax.plot(s.pseudotime_rank,s.patient_equal_mean,color=colour,lw=1.6)
   obs=bs[(bs.feature==g)&(bs.method==method)].sort_values('bin')
   ax.scatter(obs.bin_midpoint,obs.patient_equal_mean,s=10,c=colour,alpha=.6,lw=0)
   ax.axhline(0,c='#CCCCCC',lw=.6,ls='--')
   ax.set(xlim=(0,1),xlabel='Within-patient pseudotime rank',ylabel='Log1p CP10K' if rna_mode else 'metaVIPER NES')
   ax.set_title(g.replace('_RNA',' RNA')+' | '+method,loc='left',fontsize=10)
   clean(ax)
 fig.suptitle(title,x=.09,ha='left',fontsize=14,y=.98)
 fig.text(.09,.03,'Line: equal-patient mean of fixed kernel fits (bandwidth 0.075). Band: patient IQR, not confidence interval.\nDots: observed patient-balanced bin means. Same smoothing rule for every feature and method.',fontsize=8,linespacing=1.5)
 save(fig,filename)
curves(genes,'Activity trends differ between the three pseudotime methods','05_E_patient_balanced_activity_trajectories_v1')
# 6-8: fit-derived cell maps. One scale per protein across all three methods.
fitlimits={g:(min(fits[(s,g)].min() for s in methods),max(fits[(s,g)].max() for s in methods)) for g in genes}
for i,method in enumerate(methods,6):
 umaps({g:fits[(method,g)] for g in genes},method+': patient-local fitted activity projected on the frozen UMAP',f'0{i}_F_{method}_fitted_activity_UMAP_v1',fitlimits)
# 9: quantitative terminal shifts and pairwise joint positivity.
fig=plt.figure(figsize=(11.6,7.3));ax=fig.add_axes([.13,.44,.35,.40]);bx=fig.add_axes([.62,.44,.30,.40])
T=terminal[terminal.feature.isin(genes)].pivot(index='feature',columns='method',values='median_patient_effect').loc[genes,methods]
vmax=max(abs(T.values).max(),.05);im=ax.imshow(T,cmap=div,vmin=-vmax,vmax=vmax,aspect='auto')
ax.set(xticks=range(3),xticklabels=methods,yticks=range(4),yticklabels=genes);ax.tick_params(length=0)
for i,g in enumerate(genes):
 for j,method in enumerate(methods):
  r=terminal[(terminal.feature==g)&(terminal.method==method)].iloc[0]
  ax.text(j,i,format(r.median_patient_effect,'.2f')+'\n'+str(r.n_positive)+'/'+str(r.n_eligible_patients)+' positive',ha='center',va='center',fontsize=8,color='white' if abs(r.median_patient_effect)>vmax*.72 else '#333333')
ax.set_title('Terminal minus initial activity',loc='left',pad=12)
prs=[('CASP9','CHMP1A'),('CASP9','PRSS8'),('CHMP1A','PRSS8')]
P=np.zeros((3,3));pairvals=[]
for i,(aa,bb) in enumerate(prs):
 for j,method in enumerate(methods):
  rr=pairs[(pairs.method==method)&(((pairs.feature_a==aa)&(pairs.feature_b==bb))|((pairs.feature_a==bb)&(pairs.feature_b==aa)))].iloc[0]
  P[i,j]=rr.fraction_both_positive;pairvals.append((i,j,rr))
bi=bx.imshow(P,cmap=seq,vmin=0,vmax=1,aspect='auto')
bx.set(xticks=range(3),xticklabels=methods,yticks=range(3),yticklabels=[aa+' / '+bb for aa,bb in prs]);bx.tick_params(length=0)
for i,j,r in pairvals:bx.text(j,i,str(r.n_both_positive)+'/'+str(r.n_patients),ha='center',va='center',fontsize=10,color='white' if r.fraction_both_positive>.65 else '#333333')
bx.set_title('Both proteins increase at the terminal',loc='left',pad=12,fontsize=10)
ca=fig.add_axes([.13,.35,.35,.025]);fig.colorbar(im,cax=ca,orientation='horizontal');ca.set_xlabel('Median patient shift (NES)',fontsize=8)
ca=fig.add_axes([.62,.35,.30,.025]);fig.colorbar(bi,cax=ca,orientation='horizontal');ca.set_xlabel('Fraction of patients with both shifts >0',fontsize=8)
txt=[]
for method in methods:
 sub=terminal[(terminal.feature.isin(genes))&(terminal.method==method)].set_index('feature')
 txt.append(method+': '+ '; '.join(g+' P='+format(sub.loc[g,'wilcoxon_p_two_sided'],'.3g')+', FDR='+format(sub.loc[g,'BH_FDR'],'.3g') for g in genes))
fig.text(.13,.22,'\n'.join(txt),fontsize=7.5,va='top',linespacing=1.7)
fig.suptitle('Terminal coactivation is not reproduced across pseudotime methods',x=.08,ha='left',fontsize=14,y=.98)
fig.text(.13,.055,'Early/late: patient-relative lower/upper 20% of pseudotime; >=5 cells per arm. Joint-positivity counts are descriptive.',fontsize=8)
save(fig,'09_G_terminal_shift_and_coactivation_v1')
# 10: mandatory type separation for TYMP/RRM2B.
curves(['TYMP_RNA','RRM2B_RNA'],'TYMP / RRM2B RNA sensitivity: no protein-activity estimate','10_S_TYMP_RRM2B_RNA_trajectories_v1',True)
v={'status':'PASS','cells':11921,'unscored_background_cells':29,'patients':20,'inferable_MR_proteins':genes,'unavailable_MR_proteins':available.loc[~available.activity_row_present,'protein'].tolist(),'source_export_max_error':max_export_error,'F1_effect_records_recomputed':len(re),'primary_group_patients':9,'primary_group_cells':len(paired),'primary_F1_cells':int((paired.dominant_factor_rank7_v1=='F1').sum()),'primary_other_cells':int((paired.dominant_factor_rank7_v1!='F1').sum()),'trajectory_correlations_recomputed':trajectory_rows,'terminal_shifts_recomputed':terminal_rows,'source_bin_max_error':bin_error,'all_original_BH_families_verified':True,'global_four_protein_r':rglobal,'kernel_bandwidth_rank_units':bw,'kernel_inputs':'Observed patient bin means weighted by bin cell counts','curve_band':'IQR across fitted patient curves, not confidence interval','fitted_maps_scales':fitlimits,'models_rerun':False}
(q/'INDEPENDENT_VALIDATION.json').write_text(json.dumps(v,indent=2))
(q/'SOFTWARE_VERSIONS.json').write_text(json.dumps({'python':sys.version,'executable':sys.executable,'numpy':np.__version__,'pandas':pd.__version__,'matplotlib':matplotlib.__version__},indent=2))
more=[base/'08_metaviper_v1/rds/metaviper_activity_entrez_by_cell_v1.rds',base/'08_metaviper_v1/rds/metaviper_cell_metadata_v1.rds',base/'08_metaviper_v1/tables/metaviper_regulator_id_map_v1.tsv']+[root/('local_records_migration_20260908_150455/references/image'+str(i)+'.png') for i in range(13,20)]
manifest=[{'path':str(p),'bytes':p.stat().st_size,'sha256':hashfile(p)} for p in dict.fromkeys(paths+more+list((w/'scripts').glob('*')))]
(q/'SOURCE_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
checks=[]
for f in sorted(out.glob('*.pdf')):
 info=subprocess.run(['pdfinfo',str(f)],capture_output=True,text=True,check=True).stdout
 assert any(l.startswith('Pages:') and l.split(':')[1].strip()=='1' for l in info.splitlines())
 subprocess.run(['pdftoppm','-scale-to','1700','-singlefile','-png',str(f),str(q/(f.stem+'_pdf'))],capture_output=True,check=True)
 checks.append({'name':f.name,'pages':1,'render':'PASS'})
assert len(checks)==10
(q/'PDF_VALIDATION.json').write_text(json.dumps(checks,indent=2))
print(json.dumps(v),flush=True)
