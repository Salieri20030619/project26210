from pathlib import Path
import json,hashlib,sys,subprocess,shutil
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch,PathPatch
from matplotlib.path import Path as MPath
from scipy.stats import rankdata
from scipy.special import rel_entr
from scipy.spatial.distance import jensenshannon
from scipy.cluster.hierarchy import linkage,dendrogram
w=Path(__file__).resolve().parents[1];root=w.parent;q=w/'qc';out=w/'figures'
base=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/07_pyscenic_v1')
source=Path('/mnt/BioDisk/wangmingyuan/project26210_part1_cnmf_pyscenic_phase26_20260814_v1/analysis_v3')
paths=[]
def read(p,**kwargs):
 p=Path(p);paths.append(p);return pd.read_csv(p,**kwargs)
def tab(p):return read(p,sep='\t')
m=tab(base/'input/cell_metadata_v1.tsv.gz').set_index('exact_cell_id')
methods=['GENIE3','GRNBoost2'];factors=[f'F{i}' for i in range(1,8)];tfs=['ILF2','HOXD3','HNF4A']
A={method:read(base/('intermediate/auc_'+method.lower()+'_all_cells_v1.csv'),index_col=0).loc[m.index] for method in methods}
assert len(m)==11950 and m.index.is_unique and m.curated_patient_id.nunique()==26 and m.study.nunique()==5
assert np.array_equal(m[factors].idxmax(axis=1),m.dominant_factor_rank7_v1)
for method,a in A.items():
 assert a.index.equals(m.index) and a.index.is_unique and np.isfinite(a).all().all() and (a>=0).all().all()
assert A['GENIE3'].shape==(11950,299) and A['GRNBoost2'].shape==(11950,314)
shared=sorted(set(A[methods[0]].columns)&set(A[methods[1]].columns));assert len(shared)==197
assert not any('HIF1A' in a.columns.str.replace(r'\([+-]\)$','',regex=True) for a in A.values())
summary=tab(source/'continuous_regulon_patient_summary_v1.tsv')
group_summary=tab(source/'cnmf_group_regulon_patient_summary_v1.tsv')
cd=tab(source/'patient_continuous_regulon_effects_v1.tsv.gz')
gd=tab(source/'patient_cnmf_group_regulon_effects_v1.tsv.gz')
cons=tab(source/'f1_driver_consensus_v1.tsv')
assert cons.tf.head(3).tolist()==tfs and not cons.high_confidence_f1_driver.any()
assert cons.secondary_support_f1_driver.sum()==2
targets=tab(w/'data/CANDIDATE_REGULON_TARGETS.tsv')
hif=tab(Path('/mnt/BioDisk/wangmingyuan/project26210_part1_hif1a_sensitivity_phase28_20260814_v1/ctx_sensitivity_v2/run_summary_v1.tsv'))
assert len(hif)==8 and not hif.hif1a_positive_present.any()
for p in ['continuous_regulon_patient_summary_v1.tsv','cnmf_group_regulon_patient_summary_v1.tsv','f1_driver_consensus_v1.tsv','prespecified_tf_availability_v1.tsv']:
 shutil.copy2(source/p,q/p)
# Recalculate patient effects from the frozen matrices; no cell-level inference.
checks=[]
for method in methods:
 a=A[method]
 for patient,sub in m.groupby('curated_patient_id',sort=True):
  for tf in tfs:
   reg=tf+'(+)';v=a.loc[sub.index,reg].values;f1=sub.F1.values
   row=cd[(cd.method==method)&(cd.tf==tf)&(cd.curated_patient_id==patient)].iloc[0]
   elig=len(sub)>=20 and np.ptp(f1)>0 and np.ptp(v)>0
   assert bool(row.eligible)==elig
   if elig:
    rho=np.corrcoef(rankdata(f1),rankdata(v))[0,1]
    assert abs(rho-row.spearman_rho)<2e-12
   else:rho=np.nan
   mask=sub.dominant_factor_rank7_v1.values=='F1';ge=mask.sum()>=10 and (~mask).sum()>=10
   gr=gd[(gd.method==method)&(gd.tf==tf)&(gd.factor=='F1')&(gd.curated_patient_id==patient)].iloc[0]
   assert bool(gr.eligible)==ge
   delta=float(np.median(v[mask])-np.median(v[~mask])) if ge else np.nan
   if ge:assert abs(delta-gr.median_auc_difference)<2e-12
   checks.append({'method':method,'tf':tf,'patient':patient,'n_cells':len(sub),'continuous_eligible':elig,'rho':rho,'F1_group_eligible':ge,'delta_AUCell':delta})
pd.DataFrame(checks).to_csv(q/'INDEPENDENT_CANDIDATE_PATIENT_RECHECK.tsv',sep='\t',index=False)
def bh(x):
 x=np.asarray(x,float);ok=np.isfinite(x);res=np.full(len(x),np.nan);v=x[ok];order=np.argsort(v)
 adj=np.minimum.accumulate((v[order]*len(v)/np.arange(1,len(v)+1))[::-1])[::-1]
 r=np.empty(len(v));r[order]=np.clip(adj,0,1);res[ok]=r;return res
for df,keys in [(summary,['method']),(group_summary,['method','factor'])]:
 for key,s in df.groupby(keys):
  assert np.allclose(bh(s.wilcoxon_p_two_sided),s.BH_FDR,equal_nan=True,atol=2e-12)
for tf in tfs:
 for method in methods:
  s=summary[(summary.method==method)&(summary.tf==tf)].iloc[0]
  vals=cd[(cd.method==method)&(cd.tf==tf)&cd.eligible].spearman_rho
  assert len(vals)==s.n_eligible_patients==20 and abs(vals.median()-s.median_patient_effect_internal)<1e-12
  g=group_summary[(group_summary.method==method)&(group_summary.tf==tf)&(group_summary.factor=='F1')].iloc[0]
  dv=gd[(gd.method==method)&(gd.tf==tf)&(gd.factor=='F1')&gd.eligible].median_auc_difference
  assert len(dv)==g.n_eligible_patients==9 and abs(dv.median()-g.median_patient_effect_internal)<1e-12
# Verify historical input hashes and source objects.
original=tab(source/'INPUT_MANIFEST_v1.tsv')
for row in original.itertuples():assert hashlib.sha256(Path(row.path).read_bytes()).hexdigest()==row.sha256
print('SOURCE_AND_PATIENT_EFFECTS_PASS',flush=True)
# Official RSS definition (natural-log JS distance), all 197 shared regulons.
rss={};max_rss_err=0
rsssource=Path('/mnt/BioDisk/ludingyu/project26210/.venv-pyscenic/lib/python3.10/site-packages/pyscenic/rss.py');paths.append(rsssource)
scope={};exec(rsssource.read_text().replace('dtype=np.float','dtype=float'),scope)
for method,a in A.items():
 rss[method]=scope['regulon_specificity_scores'](a[shared],m.dominant_factor_rank7_v1).loc[factors]
 # Independent vectorized KL-to-mixture identity.
 prob=a[shared].values/a[shared].values.sum(axis=0)
 for fac in factors:
  labels=(m.dominant_factor_rank7_v1==fac).to_numpy(float);labels/=labels.sum()
  mixture=(prob+labels[:,None])/2
  re=1-np.sqrt(.5*(rel_entr(prob,mixture).sum(axis=0)+rel_entr(labels[:,None],mixture).sum(axis=0)))
  max_rss_err=max(max_rss_err,float(np.max(np.abs(re-rss[method].loc[fac].values))))
 assert np.isfinite(rss[method]).all().all()
 rss[method].to_csv(q/(method+'_RSS_ALL_197.tsv'),sep='\t')
assert max_rss_err<1e-11
# Fixed selection: top3 averaged within-method RSS ranks per F1-F7, plus 3 historical candidates.
avg_rank=(rss[methods[0]].rank(axis=1,ascending=False,method='average')+rss[methods[1]].rank(axis=1,ascending=False,method='average'))/2
selection=[];selected=[]
for fac in factors:
 top=avg_rank.loc[fac].sort_values(kind='stable').head(3).index.tolist()
 for reg in top:
  selection.append({'factor':fac,'regulon':reg,'criterion':'top3 mean RSS rank across two methods','mean_rank':avg_rank.loc[fac,reg]})
  if reg not in selected:selected.append(reg)
for tf in tfs:
 if tf+'(+)' not in selected:selected.append(tf+'(+)')
pd.DataFrame(selection).to_csv(q/'HEATMAP_SELECTION.tsv',sep='\t',index=False)
# Patient-balanced cNMF-group medians; no invented response/resistance strata.
profiles=[];profile_stats=[]
for (patient,fac),sub in m.groupby(['curated_patient_id','dominant_factor_rank7_v1'],sort=True):
 profile_stats.append({'patient':patient,'study':sub.study.iloc[0],'factor':fac,'cells':len(sub),'included':len(sub)>=10})
 if len(sub)<10:continue
 for method,a in A.items():
  row={'patient':patient,'study':sub.study.iloc[0],'factor':fac,'method':method,'cells':len(sub)}
  row.update(a.loc[sub.index,selected].median().to_dict());profiles.append(row)
pf=pd.DataFrame(profiles);pf.to_csv(w/'data/PATIENT_FACTOR_AUC_MEDIANS.tsv',sep='\t',index=False)
pd.DataFrame(profile_stats).to_csv(q/'ALL_PATIENT_FACTOR_INCLUSION.tsv',sep='\t',index=False)
group_medians={};Z={}
for method in methods:
 gm=pf[pf.method==method].groupby('factor')[selected].median().reindex(factors).T
 assert gm.notna().all().all()
 group_medians[method]=gm
 sd=gm.std(axis=1,ddof=0);Z[method]=gm.sub(gm.mean(axis=1),axis=0).div(sd.where(sd>0,1),axis=0)
 gm.to_csv(q/(method+'_PATIENT_BALANCED_AUC.tsv'),sep='\t')
Zall=np.column_stack([np.column_stack([Z[method][fac].values for method in methods]) for fac in factors])
tree=linkage(Zall,method='average',metric='euclidean',optimal_ordering=True)
leaves=dendrogram(tree,no_plot=True)['leaves'];regs=[selected[i] for i in leaves]
# Labels, colours and publication exports.
fc=dict(zip(factors,['#BA7868','#7295BE','#79A697','#C9AD63','#A58BB4','#8DB6C3','#C78FAD']))
mc={'GENIE3':'#AB7764','GRNBoost2':'#5E91AB'}
div=LinearSegmentedColormap.from_list('regulon_div',['#5E936F','#F7F6F2','#927093'])
activity=LinearSegmentedColormap.from_list('activity',['#E8EAED','#A1B0C4','#345D8C'])
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.labelsize':9,'axes.titlesize':11,'axes.linewidth':.65,'pdf.fonttype':42,'savefig.facecolor':'white','xtick.major.width':.65,'ytick.major.width':.65})
def clean(ax):
 ax.spines[['top','right']].set_visible(False);ax.tick_params(length=3)
def save(fig,name):
 for ext in ['png','pdf']:fig.savefig(out/(name+'.'+ext),dpi=360,bbox_inches='tight',pad_inches=.15)
 plt.close(fig)
# A: both algorithms, all seven program ranks, no selective exclusion of low-scoring factors.
for mi,method in enumerate(methods,1):
 fig,axes=plt.subplots(2,4,figsize=(13.4,7.0));fig.subplots_adjust(left=.065,right=.985,top=.88,bottom=.09,hspace=.43,wspace=.31)
 for j,fac in enumerate(factors):
  ax=axes.flat[j];s=rss[method].loc[fac].sort_values(ascending=False,kind='stable')
  ax.scatter(np.arange(1,198),s.values,s=8,c='#698DAA',linewidths=0,alpha=.75)
  top=s.iloc[:3];ax.scatter([1,2,3],top.values,s=20,c='#B65F60',linewidths=0,zorder=3)
  # Three label rails avoid overlapping the clustered top ranks.
  for k,(reg,val) in enumerate(top.items()):
   ax.annotate(reg,xy=(k+1,val),xytext=(.37,.90-k*.115),textcoords='axes fraction',fontsize=8,
    arrowprops={'arrowstyle':'-','lw':.65,'color':'#A8A8A8'},va='center')
  ax.set_title(fac+'  |  '+format(int((m.dominant_factor_rank7_v1==fac).sum()),',')+' cells',loc='left',color=fc[fac])
  ax.set(xlim=(-4,204),xlabel='Regulon rank',ylabel='Regulon specificity score')
  ax.set_ylim(float(s.min())-.025,float(s.max())+.055);clean(ax)
 axes.flat[7].axis('off')
 axes.flat[7].text(.03,.95,'197 shared regulons\n11,950 malignant cells\n\nRSS = 1 - Jensen-Shannon distance\n\nRanks are descriptive.\nColour highlights the top three\nregulons in each program.',va='top',linespacing=1.65,fontsize=9)
 fig.suptitle(method+' regulon specificity across cNMF programs',x=.06,ha='left',fontsize=15,y=.98)
 save(fig,f'0{mi}_A_RSS_{method}_F1_F7_v1')
# B: clustered heatmap with paired algorithm columns and actual patient counts.
fig=plt.figure(figsize=(9.9,8.8));ax=fig.add_axes([.17,.19,.59,.61])
im=ax.imshow(Zall[leaves],cmap=div,vmin=-2,vmax=2,aspect='auto',interpolation='nearest')
ax.set_yticks(np.arange(len(regs)),regs,fontsize=8);ax.yaxis.tick_right();ax.tick_params(axis='y',length=0,pad=5)
ax.set_xticks(np.arange(14),['G3','GB2']*7,fontsize=8);ax.tick_params(axis='x',length=0)
for k in np.arange(1.5,13,2):ax.axvline(k,color='white',lw=2)
for s in ax.spines.values():s.set_visible(False)
for tick,reg in zip(ax.get_yticklabels(),regs):
 if reg[:-3] in tfs:tick.set_fontweight('bold')
ta=fig.add_axes([.04,.19,.115,.61]);dendrogram(tree,orientation='left',no_labels=True,color_threshold=0,above_threshold_color='#888888',ax=ta);ta.invert_yaxis();ta.axis('off')
ann=fig.add_axes([.17,.815,.59,.022]);ann.imshow(np.repeat(np.arange(7),2)[None,:],cmap=ListedColormap(list(fc.values())),aspect='auto',vmin=0,vmax=6)
ann.set_xticks(np.arange(.5,14,2),factors);ann.xaxis.tick_top();ann.tick_params(length=0);ann.set_yticks([])
nprof=pf[pf.method=='GENIE3'].groupby('factor').patient.nunique().reindex(factors)
for k,fac in enumerate(factors):ax.text(k*2+.5,len(regs)+.85,'n='+str(nprof[fac]),ha='center',fontsize=8)
ax.set_xlabel('cNMF program | G3: GENIE3, GB2: GRNBoost2',labelpad=22)
ca=fig.add_axes([.65,.91,.22,.018]);fig.colorbar(im,cax=ca,orientation='horizontal',ticks=[-2,0,2]);ca.set_title('Within-method row z-score',fontsize=8)
fig.text(.04,.965,'Patient-balanced regulon activity across cNMF programs',fontsize=14,ha='left')
fig.text(.17,.085,'Median AUCell per patient-program (>=10 cells), then median across patients.\nRows standardized across F1-F7 within each algorithm. Bold: the three historical candidates.',fontsize=8,linespacing=1.7)
save(fig,'03_B_cNMF_patient_balanced_regulon_heatmap_v1')
# C: full-statistic continuous candidate effects. Raw points show biological consistency.
fig,axes=plt.subplots(1,2,figsize=(11.0,4.9));fig.subplots_adjust(left=.10,right=.97,top=.79,bottom=.22,wspace=.50)
continuous_display=[]
for ax,method in zip(axes,methods):
 for k,tf in enumerate(tfs):
  rows=cd[(cd.method==method)&(cd.tf==tf)&cd.eligible].sort_values('curated_patient_id')
  vals=rows.spearman_rho.values;s=summary[(summary.method==method)&(summary.tf==tf)].iloc[0]
  jitter=np.linspace(-.11,.11,len(vals))
  ax.scatter(vals,np.full(len(vals),k)+jitter,s=21,c=mc[method],alpha=.6,linewidths=0)
  ax.scatter(np.median(vals),k,s=55,c='#303030',marker='D',edgecolor='white',linewidth=.6,zorder=4)
  lab='rho='+format(s.median_patient_effect_internal,'.3f')+'; P='+format(s.wilcoxon_p_two_sided,'.4g')+'; FDR='+format(s.BH_FDR,'.4g')
  ax.text(-.94,k+.33,lab,fontsize=8,color='#333333')
  continuous_display.append(s.to_dict())
 ax.axvline(0,ls='--',c='#BBBBBB',lw=.8)
 ax.set(xlim=(-1,1),ylim=(2.62,-.5),yticks=range(3),yticklabels=tfs,xlabel='Within-patient F1-AUCell Spearman rho')
 ax.set_title(method+'  |  20 patients',loc='left',color=mc[method]);clean(ax)
fig.suptitle('Three candidate regulons show positive F1 associations in both algorithms',fontsize=14,x=.10,ha='left',y=.98)
fig.text(.10,.075,'One point per eligible patient; black diamond: median. FDR is calculated over all regulons within each algorithm.',fontsize=8)
save(fig,'04_C_candidate_F1_continuous_associations_v1')
pd.DataFrame(continuous_display).to_csv(q/'CANDIDATE_CONTINUOUS_DISPLAY.tsv',sep='\t',index=False)
# D: n=9 exploratory grouped effects, no stars and no p-value-led titles.
fig,axes=plt.subplots(1,2,figsize=(11,4.9));fig.subplots_adjust(left=.10,right=.97,top=.80,bottom=.23,wspace=.50)
group_display=[]
for ax,method in zip(axes,methods):
 chunks=[]
 for tf in tfs:
  rr=gd[(gd.method==method)&(gd.tf==tf)&(gd.factor=='F1')&gd.eligible].sort_values('curated_patient_id')
  chunks.append(rr.median_auc_difference.values)
 bound=max(abs(np.concatenate(chunks)).max()*1.10,.008)
 for k,(tf,vals) in enumerate(zip(tfs,chunks)):
  jitter=np.linspace(-.11,.11,len(vals))
  ax.scatter(vals,np.full(len(vals),k)+jitter,s=29,c=mc[method],alpha=.75,linewidths=0)
  ax.scatter(np.median(vals),k,s=55,c='#303030',marker='D',edgecolor='white',linewidth=.6,zorder=4)
  s=group_summary[(group_summary.method==method)&(group_summary.tf==tf)&(group_summary.factor=='F1')].iloc[0]
  ax.text(-bound*.96,k+.34,'Median delta='+format(np.median(vals),'.4f')+' | '+str(s.n_positive)+' positive, '+str(s.n_negative)+' negative, '+str(s.n_zero)+' ties',fontsize=8)
  group_display.append(s.to_dict())
 ax.axvline(0,ls='--',c='#BBBBBB',lw=.8)
 ax.set(xlim=(-bound,bound),ylim=(2.62,-.5),yticks=range(3),yticklabels=tfs,xlabel='Patient median AUCell: F1-dominant minus other cells')
 ax.set_title(method+'  |  9 paired patients',loc='left',color=mc[method]);ax.ticklabel_format(axis='x',style='sci',scilimits=(-2,2));clean(ax)
fig.suptitle('F1-dominant group contrasts provide exploratory candidate evidence',fontsize=14,x=.10,ha='left',y=.98)
fig.text(.10,.075,'Both arms have >=10 cells per patient. One point per patient; diamond: median. Full P values and FDR are retained in the server tables.',fontsize=8)
save(fig,'05_D_candidate_F1_group_exploratory_v1')
pd.DataFrame(group_display).to_csv(q/'CANDIDATE_GROUP_COMPLETE_STATS.tsv',sep='\t',index=False)
# E: circular, method-labelled true edges, cap 10 targets in each support category.
circles={};edge_rows=[];network_counts=[]
for tf in tfs:
 sub=targets[(targets.tf==tf)&(targets.target!=tf)].copy()
 piv=sub.pivot(index='target',columns='method',values='weight')
 both=piv.notna().all(axis=1);onlyg=piv.GENIE3.notna()&piv.GRNBoost2.isna();onlyb=piv.GRNBoost2.notna()&piv.GENIE3.isna()
 ranks=piv.rank(pct=True)
 categories={'Both methods':piv.index[both],'GENIE3 only':piv.index[onlyg],'GRNBoost2 only':piv.index[onlyb]}
 selected_targets=[];cats={}
 for cat,genes in categories.items():
  ss=ranks.loc[genes].mean(axis=1).sort_values(ascending=False,kind='stable').head(10)
  selected_targets.extend(ss.index)
  for gene in ss.index:cats[gene]=cat
 for gene in selected_targets:
  for method in methods:
   if pd.notna(piv.loc[gene,method]):edge_rows.append({'tf':tf,'target':gene,'method':method,'weight':piv.loc[gene,method],'support':cats[gene]})
 circles[tf]=(selected_targets,cats,piv)
 network_counts.append({'tf':tf,'GENIE3_targets_excluding_self':int(piv.GENIE3.notna().sum()),'GRNBoost2_targets_excluding_self':int(piv.GRNBoost2.notna().sum()),'shared_targets_excluding_self':int(both.sum()),'displayed_targets':len(selected_targets),'TYMP_GENIE3':bool('TYMP' in piv.index and pd.notna(piv.loc['TYMP','GENIE3'])),'TYMP_GRNBoost2':bool('TYMP' in piv.index and pd.notna(piv.loc['TYMP','GRNBoost2']))})
pd.DataFrame(edge_rows).to_csv(q/'CIRCULAR_DISPLAY_EDGES.tsv',sep='\t',index=False)
pd.DataFrame(network_counts).to_csv(q/'NETWORK_COUNTS.tsv',sep='\t',index=False)
catc={'Both methods':'#94758D','GENIE3 only':'#BA8067','GRNBoost2 only':'#6398AA'}
fig,axes=plt.subplots(1,3,figsize=(14.8,6.1));fig.subplots_adjust(left=.01,right=.99,top=.82,bottom=.17,wspace=.12)
for ax,tf,counts in zip(axes,tfs,network_counts):
 genes,cats,piv=circles[tf]
 angles=np.linspace(np.deg2rad(145),np.deg2rad(395),len(genes))
 pos={g:np.array([np.cos(t),np.sin(t)]) for g,t in zip(genes,angles)}
 hubs={'GENIE3':np.array([-.28,.97]),'GRNBoost2':np.array([.28,.97])}
 for g in genes:
  for method in methods:
   if pd.isna(piv.loc[g,method]):continue
   a=hubs[method];b=pos[g]
   path=MPath([a,a*.45,b*.45,b],[MPath.MOVETO,MPath.CURVE4,MPath.CURVE4,MPath.CURVE4])
   ax.add_patch(PathPatch(path,fc='none',ec=mc[method],lw=.75,alpha=.75,zorder=1))
 for g,t in zip(genes,angles):
  xx,yy=pos[g];ax.scatter(xx,yy,s=24,c=catc[cats[g]],edgecolors='white',lw=.4,zorder=3)
  deg=np.rad2deg(t)%360;left=90<deg<270;rotation=deg+180 if left else deg
  ax.text(xx*1.06,yy*1.06,g,rotation=rotation,rotation_mode='anchor',ha='right' if left else 'left',va='center',fontsize=7.5,color='#454545')
 for method,p in hubs.items():
  ax.scatter(*p,s=80,c=mc[method],edgecolor='white',lw=.6,zorder=4)
  ax.text(p[0],p[1]+.12,method,ha='center',fontsize=8,color=mc[method])
 ax.text(0,0,tf,ha='center',va='center',fontsize=14,fontweight='bold',bbox={'facecolor':'white','edgecolor':'none','alpha':.9,'pad':3})
 ax.set(xlim=(-1.43,1.43),ylim=(-1.38,1.5));ax.set_aspect('equal');ax.axis('off')
 ax.set_title(tf+'(+)  |  '+str(counts['displayed_targets'])+' displayed targets',fontsize=11,pad=14)
 ax.text(.5,-.035,'All non-self targets: G3 '+str(counts['GENIE3_targets_excluding_self'])+' | GB2 '+str(counts['GRNBoost2_targets_excluding_self'])+' | shared '+str(counts['shared_targets_excluding_self']),transform=ax.transAxes,ha='center',fontsize=8)
fig.suptitle('Candidate regulon networks retain method-specific target support',fontsize=15,x=.04,ha='left',y=.98)
fig.legend(handles=[Line2D([],[],ls='',marker='o',c=c,label=cat,ms=6) for cat,c in catc.items()],loc='lower center',bbox_to_anchor=(.5,.075),frameon=False,ncol=3,fontsize=9)
fig.text(.5,.025,'Top 10 targets per support category by mean within-method weight percentile. Curves are frozen TF-target predictions.',ha='center',fontsize=8)
save(fig,'06_E_candidate_circular_regulon_networks_v1')
# F: original UMAP, each algorithm has its own AUCell scale.
xyfile=root/'phase123_section2_2_reference_redraw_remote_v1/data/plot_values_server_only.npz';paths.append(xyfile)
z=np.load(xyfile);xy=pd.DataFrame(z['xy'],index=z['cells']).loc[m.index].values
assert len(xy)==11950 and np.isfinite(xy).all()
fig,axes=plt.subplots(2,3,figsize=(11.5,8.8));fig.subplots_adjust(left=.07,right=.98,top=.91,bottom=.07,hspace=.28,wspace=.20)
for ri,method in enumerate(methods):
 for ci,tf in enumerate(tfs):
  ax=axes[ri,ci];v=A[method][tf+'(+)'].values;order=np.argsort(v,kind='stable')
  im=ax.scatter(xy[order,0],xy[order,1],c=v[order],s=2.8,cmap=activity,vmin=0,vmax=v.max(),lw=0,rasterized=True)
  ax.set_title(tf+'(+)  |  '+method,loc='left',fontsize=10);ax.set_aspect('equal');ax.axis('off')
  cb=fig.colorbar(im,ax=ax,orientation='horizontal',fraction=.033,pad=.025,shrink=.75);cb.set_label('AUCell regulon enrichment',fontsize=8);cb.ax.tick_params(labelsize=7,length=2)
fig.suptitle('Candidate regulon enrichment in the frozen malignant-cell UMAP',fontsize=14,x=.07,ha='left',y=.985)
save(fig,'07_F_candidate_AUCell_UMAP_v1')
validation={'status':'PASS','cells':11950,'patients':26,'methods_regulons':{k:len(a.columns) for k,a in A.items()},'shared_regulons':197,'candidate_continuous_n':20,'candidate_F1_group_n':9,'all_input_hashes_match':True,'candidate_patient_rows_recomputed':len(checks),'BH_recomputed_by_original_families':True,'RSS_independent_max_error':max_rss_err,'heatmap_regulons':len(selected),'heatmap_profiles':len(pf)//2,'heatmap_patient_counts':nprof.to_dict(),'HIF1A_available':False,'HIF1A_sensitivity_recovered':0,'network_counts':network_counts,'scientific_model_rerun':False}
(q/'INDEPENDENT_VALIDATION.json').write_text(json.dumps(validation,indent=2))
(q/'SOFTWARE_VERSIONS.json').write_text(json.dumps({'python':sys.version,'executable':sys.executable,'numpy':np.__version__,'pandas':pd.__version__,'matplotlib':matplotlib.__version__},indent=2))
more=[base/('intermediate/regulons_'+method.lower()+'_v1.dat') for method in methods]
manifest=[]
for f in dict.fromkeys(paths+more+[Path(__file__)]):
 h=hashlib.sha256()
 with f.open('rb') as rd:
  for b in iter(lambda:rd.read(1048576),b''):h.update(b)
 manifest.append({'path':str(f),'bytes':f.stat().st_size,'sha256':h.hexdigest()})
(q/'SOURCE_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
pdfchecks=[]
for f in sorted(out.glob('*.pdf')):
 info=subprocess.run(['pdfinfo',str(f)],capture_output=True,text=True,check=True).stdout
 assert any(l.startswith('Pages:') and l.split(':')[1].strip()=='1' for l in info.splitlines())
 subprocess.run(['pdftoppm','-scale-to','1800','-singlefile','-png',str(f),str(q/(f.stem+'_pdf'))],capture_output=True,check=True)
 pdfchecks.append({'name':f.name,'pages':1,'render':'PASS'})
assert len(pdfchecks)==7
(q/'PDF_VALIDATION.json').write_text(json.dumps(pdfchecks,indent=2))
print(json.dumps(validation),flush=True)
