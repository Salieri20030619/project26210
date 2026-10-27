from pathlib import Path
import json,subprocess
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.ticker import MaxNLocator
from scipy.stats import gaussian_kde,rankdata
w=Path(__file__).resolve().parents[1];q=w/'qc';out=w/'figures_2_5';old=w.parent/'phase128_section2_5_metaviper_redraw_remote_v1'
points=pd.read_csv(old/'qc/m11_mr_effect_activity_points_v1.tsv',sep='\t')
a=pd.read_csv(old/'data/SELECTED_METAVIPER_ACTIVITY.tsv',sep='\t',index_col=0)
m=pd.read_csv('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/07_pyscenic_v1/input/cell_metadata_v1.tsv.gz',sep='\t',index_col=0).loc[a.index]
summ=pd.read_csv(old/'qc/m12_f1_activity_patient_summary_v1.tsv',sep='\t')
genes=['CASP9','CHMP1A','PRSS8','PXN'];assert set(points.protein)==set(genes) and len(a)==11921
cmap=LinearSegmentedColormap.from_list('reference_purple',['#FFFFFF','#EEE9F5','#D4C8E8','#AA96CF','#7D63B1','#54398D','#2C155C'])
levels=np.array([.025,.06,.12,.22,.36,.52,.70,.86,1.000001])
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'savefig.facecolor':'white'})
def kdegrid(x,y,xlim,ylim,weights=None,mult=1.,size=180):
 xy=np.vstack([x,y]);kde=gaussian_kde(xy,weights=weights,bw_method='scott');kde.set_bandwidth(kde.factor*mult)
 gx=np.linspace(*xlim,size);gy=np.linspace(*ylim,size);xx,yy=np.meshgrid(gx,gy)
 density=kde(np.vstack([xx.ravel(),yy.ravel()])).reshape(xx.shape)
 assert np.isfinite(density).all() and density.max()>0
 return gx,gy,density/density.max(),float(kde.factor),float(kde.neff)
def contours(ax,gx,gy,z):
 ax.contourf(gx,gy,z,levels=levels,cmap=cmap,vmin=0,vmax=1,extend='neither')
 ax.contour(gx,gy,z,levels=levels[:-1],colors='#493A63',linewidths=.65,alpha=.75)
 ax.set_axisbelow(True);ax.tick_params(length=3,labelsize=9)
def save(fig,name):
 for ext in ['png','pdf']:fig.savefig(out/(name+'.'+ext),dpi=360)
 plt.close(fig);subprocess.run(['pdftoppm','-r','125','-singlefile','-png',str(out/(name+'.pdf')),str(q/(name+'_pdf'))],check=True)
 print('SAVED',name,flush=True)
x=points.mr_beta_log_or.values;y=points.patient_equal_mean_activity.values;r=float(np.corrcoef(x,y)[0,1]);assert abs(r-(-.648335))<1e-5
pm=a[genes].groupby(a.curated_patient_id).mean().mean()
assert np.allclose(points.set_index('protein').loc[genes,'patient_equal_mean_activity'],pm,atol=1e-12)
probe=gaussian_kde(np.vstack([x,y]),bw_method='scott');pads=3.2*np.sqrt(np.diag(probe.covariance))
xlim=(x.min()-pads[0],x.max()+pads[0]);ylim=(y.min()-pads[1],y.max()+pads[1])
gx,gy,z,bw,neff=kdegrid(x,y,xlim,ylim)
np.savez_compressed(w/'data/MR_DENSITY_GRID.npz',x=gx,y=gy,relative_density=z,observed_x=x,observed_y=y)
fig,ax=plt.subplots(figsize=(7.8,7.2));fig.subplots_adjust(left=.17,right=.79,top=.80,bottom=.20);contours(ax,gx,gy,z)
ax.scatter(x,y,color='#26212B',s=30,edgecolors='white',linewidths=.7,zorder=4)
for rr in points.itertuples():
 offsets={'PXN':(8,10),'CHMP1A':(8,-17),'PRSS8':(-60,12),'CASP9':(8,12)}
 ax.annotate(rr.protein,(rr.mr_beta_log_or,rr.patient_equal_mean_activity),xytext=offsets[rr.protein],textcoords='offset points',fontsize=10,fontweight='bold',zorder=5,bbox={'facecolor':'white','alpha':.8,'edgecolor':'none','pad':1.5})
ax.set_xlim(xlim);ax.set_ylim(ylim);ax.set_xlabel('MR effect (log OR)',fontsize=11);ax.set_ylabel('Patient-equal mean metaVIPER NES',fontsize=11)
ax.text(.04,.97,f'Pearson r = {r:+.3f}\n4 proteins',transform=ax.transAxes,va='top',fontsize=11,bbox={'facecolor':'white','alpha':.8,'edgecolor':'none','pad':3})
cax=fig.add_axes([.84,.34,.025,.32]);cb=fig.colorbar(ScalarMappable(norm=Normalize(0,1),cmap=cmap),cax=cax,ticks=[0,1]);cb.ax.set_yticklabels(['Low','High']);cb.set_label('Relative density',fontsize=10)
fig.text(.055,.958,'MR effects versus inferred protein activity',fontsize=16,fontweight='bold')
fig.text(.055,.905,'Reference-style density contours with the four observed protein points',fontsize=10,color='#53616D')
fig.text(.055,.049,'Each point is one MR protein; activity is averaged equally across 20 patients.\nContours smooth these four points (Gaussian KDE, Scott bandwidth); r is descriptive.\nThis small target set does not support an inferential correlation test.',fontsize=8.5,linespacing=1.5)
save(fig,'01_MR_effect_activity_density_contours_v2')
# Preserve bandwidth sensitivity in QC, not as extra biological evidence.
fig,axs=plt.subplots(1,3,figsize=(12,4))
for ax,mult in zip(axs,[.7,1,1.4]):
 xx,yy,zz,_,_=kdegrid(x,y,xlim,ylim,mult=mult,size=120);contours(ax,xx,yy,zz);ax.scatter(x,y,s=12,c='black');ax.set_title(f'Scott bandwidth x {mult}');ax.set_xlabel('MR log OR')
fig.tight_layout();fig.savefig(q/'MR_KDE_BANDWIDTH_SENSITIVITY.png',dpi=150);plt.close(fig)
# Protein-specific cell distributions, patients contribute equal total KDE weight.
counts=m.curated_patient_id.value_counts();weights=1/m.curated_patient_id.map(counts).to_numpy(float);weights/=weights.sum()
assert np.allclose(pd.Series(weights,index=m.index).groupby(m.curated_patient_id).sum(),1/20)
fig,axs=plt.subplots(4,1,figsize=(7.9,14.5));fig.subplots_adjust(left=.15,right=.69,top=.915,bottom=.12,hspace=.32)
records=[]
for ax,gene in zip(axs,genes):
 x=m.F1.values;y=a[gene].values;span=y.max()-y.min();ylim=(y.min()-.025*span,y.max()+.025*span)
 gx,gy,z,bw1,neff1=kdegrid(x,y,(-.2,1.2),ylim,weights=weights,size=180)
 np.savez_compressed(w/'data'/('F1_'+gene+'_DENSITY_GRID.npz'),x=gx,y=gy,relative_density=z)
 contours(ax,gx,gy,z);ax.set_xlim(-.2,1.2);ax.axvline(0,color='#AAAAAA',lw=.6,ls='--');ax.axvline(1,color='#AAAAAA',lw=.6,ls='--');ax.set_ylim(ylim);ax.set_xticks([0,.5,1]);ax.yaxis.set_major_locator(MaxNLocator(4));ax.set_ylabel('metaVIPER NES');ax.set_title(gene,loc='left',fontsize=12,fontweight='bold')
 mx=np.sum(weights*x);my=np.sum(weights*y);rc=float(np.sum(weights*(x-mx)*(y-my))/np.sqrt(np.sum(weights*(x-mx)**2)*np.sum(weights*(y-my)**2)))
 ax.text(.035,.965,f'Weighted cell r = {rc:+.3f}',transform=ax.transAxes,va='top',fontsize=9,bbox={'facecolor':'white','alpha':.8,'edgecolor':'none','pad':2})
 rr=summ[(summ.protein==gene)&(summ.endpoint=='continuous_F1_spearman')].iloc[0]
 patient_r=[]
 for pid,sub in m.groupby('curated_patient_id'):
  patient_r.append(np.corrcoef(rankdata(sub.F1),rankdata(a.loc[sub.index,gene]))[0,1])
 assert abs(np.median(patient_r)-rr.median_patient_effect)<1e-10
 bounds=ax.get_position();fig.text(.75,bounds.y0+bounds.height*.70,f'Patient-level evidence\nn = {int(rr.n_eligible_patients)}\nMedian rho = {rr.median_patient_effect:+.3f}\nPositive: {int(rr.n_positive)}/{int(rr.n_eligible_patients)}\nP = {rr.wilcoxon_p_two_sided:.3g}\nFDR = {rr.BH_FDR:.3g}',va='top',fontsize=9,linespacing=1.6)
 records.append({'protein':gene,'cells':len(x),'patients':20,'weighted_cell_Pearson_r':rc,'KDE_bandwidth':bw1,'KDE_weight_neff':neff1,'patient_rho':float(rr.median_patient_effect),'patient_P':float(rr.wilcoxon_p_two_sided),'patient_FDR':float(rr.BH_FDR)})
axs[-1].set_xlabel('F1 normalized usage',fontsize=11)
cax=fig.add_axes([.76,.13,.15,.012]);cb=fig.colorbar(ScalarMappable(norm=Normalize(0,1),cmap=cmap),cax=cax,orientation='horizontal',ticks=[0,1]);cb.ax.set_xticklabels(['Low','High']);cb.set_label('Relative density',fontsize=9)
fig.text(.055,.970,'F1 versus inferred protein activity',fontsize=17,fontweight='bold')
fig.text(.055,.943,'Protein-specific density contours; distinct from the MR-effect analysis',fontsize=10,color='#53616D')
fig.text(.055,.043,'11,921 cells from 20 patients | each patient contributes equal total density weight\nGaussian KDE; dashed lines mark observed F1 support [0,1]; outside tails are smoothing only.\nCell r is descriptive. Patient P/FDR retain the original four-MR-protein family.',fontsize=9,linespacing=1.5)
save(fig,'02_F1_protein_activity_density_contours_v2')
pd.DataFrame(records).to_csv(q/'F1_DENSITY_STATISTICS.tsv',sep='\t',index=False)
validation={'status':'PASS','MR_proteins':4,'MR_Pearson_r':r,'MR_activity_patient_count':20,'MR_KDE_bandwidth':bw,'MR_KDE_interpretation':'descriptive smoothing of 4 points; no inferential P','F1_activity_cells':11921,'F1_activity_patients':20,'F1_KDE_patient_equal_weights':True,'F1_original_patient_effects_match':True,'F1_cell_P_values_computed':False,'activity_values_changed':False,'unavailable_MR_proteins_not_imputed':['TYMP','RRM2B','KDELC2/POGLUT3','TMEM25']}
(q/'SECTION25_VALIDATION.json').write_text(json.dumps(validation,indent=2));print(json.dumps(validation),flush=True)
