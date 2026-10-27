from pathlib import Path
import json,hashlib,sys,subprocess,platform
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,TwoSlopeNorm,LogNorm
from matplotlib.lines import Line2D
from scipy.stats import rankdata
from PIL import Image,ImageOps,ImageDraw
w=Path(__file__).resolve().parents[1];root=w.parent;q=w/'qc';out=w/'figures'
src=Path('/mnt/BioDisk/wangmingyuan/project26210_part1_scdrs_phase25_20260814_v1')
paths=[]
def read(p):
 paths.append(p);return pd.read_csv(p,sep='\t')
def js(p):
 paths.append(p);return json.loads(p.read_text())
d=read(src/'analysis_v1/scdrs_factor_tymp_cell_source_v1.tsv.gz').set_index('exact_cell_id')
pt=read(src/'analysis_v1/patient_within_association_v1.tsv')
ov=js(src/'analysis_v1/overall_patient_association_v1.json')
desc=js(src/'analysis_v1/cell_level_descriptive_v1.json')
scores=read(src/'compute_v2/scores/RCC_GWAS_Catalog_loci.score.gz').set_index('cell') if False else pd.read_csv(src/'compute_v2/scores/RCC_GWAS_Catalog_loci.score.gz',sep='\t',index_col=0)
paths.append(src/'compute_v2/scores/RCC_GWAS_Catalog_loci.score.gz')
usage=read(Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/05_cnmf/full_ranks4_8_v1/consensus_usage_rank7_v1.tsv.gz')).set_index('exact_cell_id')
xyfile=root/'phase123_section2_2_reference_redraw_remote_v1/data/plot_values_server_only.npz';paths.append(xyfile);coords=np.load(xyfile)
assert len(d)==11950 and d.index.is_unique and set(d.index)==set(usage.index)==set(scores.index)==set(coords['cells'])
assert np.isfinite(d[['F1','scdrs_norm_score','UMAP1','UMAP2']]).all().all()
assert np.allclose(d.F1,usage.loc[d.index,'F1'],atol=1e-12)
assert np.allclose(d.scdrs_norm_score,scores.loc[d.index,'norm_score'],atol=1e-8)
xy=pd.DataFrame(coords['xy'],index=coords['cells']).loc[d.index].values
assert np.allclose(xy,d[['UMAP1','UMAP2']],atol=2e-6)
re=[]
for row in pt.itertuples():
 a=d[d.curated_patient_id==row.curated_patient_id]
 eligible=len(a)>=20 and a.F1.nunique()>1 and a.scdrs_norm_score.nunique()>1
 assert len(a)==row.n_cells and eligible==bool(row.eligible_primary)
 rho=float(np.corrcoef(rankdata(a.F1),rankdata(a.scdrs_norm_score))[0,1]) if eligible else None
 if eligible:assert abs(rho-row.spearman_rho_f1_vs_scdrs)<1e-12
 re.append({'patient':row.curated_patient_id,'study':row.study,'cells':len(a),'eligible':eligible,'rho':rho})
rhos=np.array([r['rho'] for r in re if r['eligible']]);assert len(rhos)==20
rng=np.random.default_rng(26210)
boot=np.median(rhos[rng.integers(0,20,size=(10000,20))],axis=1)
ci=np.quantile(boot,[.025,.975])
perm=np.median(rng.choice(np.array([-1.,1.]),size=(9999,20))*rhos,axis=1)
p=(1+(np.abs(perm)>=abs(np.median(rhos))).sum())/10000
assert abs(p-ov['signflip_p_two_sided'])<1e-12
assert np.allclose(ci,[ov['bootstrap_ci95_low'],ov['bootstrap_ci95_high']],atol=1e-12)
pearson=float(np.corrcoef(d.F1,d.scdrs_norm_score)[0,1]);assert abs(pearson-desc['pooled_cell_pearson_r_descriptive_only'])<1e-12
assert (d.scdrs_qval_bh_descriptive<.2).sum()==0
pd.DataFrame(re).to_csv(q/'INDEPENDENT_PATIENT_RECHECK.tsv',sep='\t',index=False)
d.to_csv(w/'data/PLOT_CELL_SOURCE.tsv.gz',sep='\t')
pt.to_csv(q/'ALL_PATIENT_STATISTICS.tsv',sep='\t',index=False)
print('IDENTITY_AND_STATISTICS_PASS',flush=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.labelsize':9,'axes.titlesize':11,'axes.linewidth':.65,'pdf.fonttype':42,'savefig.facecolor':'white','xtick.major.width':.65,'ytick.major.width':.65})
studies=['GSE156632','GSE159115','GSE178481','GSE207493','GSE242299']
sd=dict(zip(studies,['#7296BC','#E3AD66','#76A893','#9F86B4','#C47F85']))
div=LinearSegmentedColormap.from_list('scdrs',['#456F99','#A7BCD0','#F1F0EE','#DAA08A','#A43D40'])
density=LinearSegmentedColormap.from_list('density',['#E5EEF3','#9ABACE','#528BA9','#205878'])
low=float(d.scdrs_norm_score.min());high=float(d.scdrs_norm_score.max())
norm=TwoSlopeNorm(vmin=low,vcenter=0,vmax=high)
v=d.scdrs_norm_score.to_numpy();x=d.F1.to_numpy()
def clean(ax):
 ax.spines[['top','right']].set_visible(False);ax.tick_params(length=3)
def umap(fig,ax):
 order=np.argsort(np.abs(v),kind='stable')
 im=ax.scatter(xy[order,0],xy[order,1],c=v[order],cmap=div,norm=norm,s=3.3,linewidths=0,rasterized=True)
 ax.set_aspect('equal');ax.axis('off');ax.set_title('RCC GWAS-locus scDRS score',loc='left',pad=12)
 ax.annotate('',xy=(.23,.06),xytext=(.05,.06),xycoords='axes fraction',arrowprops={'arrowstyle':'->','lw':.9,'color':'#343434'})
 ax.annotate('',xy=(.05,.24),xytext=(.05,.06),xycoords='axes fraction',arrowprops={'arrowstyle':'->','lw':.9,'color':'#343434'})
 ax.text(.15,.00,'UMAP1',transform=ax.transAxes,fontsize=8,ha='center')
 ax.text(.005,.15,'UMAP2',transform=ax.transAxes,fontsize=8,rotation=90,va='center')
 cb=fig.colorbar(im,ax=ax,orientation='horizontal',fraction=.038,pad=.055,shrink=.73)
 cb.set_label('scDRS normalized score',fontsize=8);cb.ax.tick_params(labelsize=8,length=2)
 ax.text(.02,1.01,'11,950 malignant cells | 26 patients',transform=ax.transAxes,fontsize=8,color='#666666')
def scatter(fig,ax):
 # Every observation is shown; density is empirical occupancy, not synthetic or smoothed scores.
 counts,ex,ey=np.histogram2d(x,v,bins=(65,65))
 ix=np.clip(np.searchsorted(ex,x,side='right')-1,0,64);iy=np.clip(np.searchsorted(ey,v,side='right')-1,0,64)
 den=counts[ix,iy];order=np.argsort(den,kind='stable')
 im=ax.scatter(x[order],v[order],c=den[order],cmap=density,norm=LogNorm(vmin=1,vmax=den.max()),s=4.1,linewidths=0,rasterized=True)
 coef=np.polyfit(x,v,1);xx=np.array([x.min(),x.max()]);ax.plot(xx,np.polyval(coef,xx),color='#B56F60',lw=1.4,zorder=5)
 ax.set(xlabel='Factor-TYMP score (F1 normalized usage)',ylabel='scDRS normalized score')
 ax.set_title('F1 and scDRS across malignant cells',loc='left',pad=12)
 ax.text(.97,.96,'Pearson r = '+format(pearson,'.4f')+'\nCell-level descriptive association',transform=ax.transAxes,ha='right',va='top',fontsize=8.5,bbox={'facecolor':'white','edgecolor':'none','alpha':.92,'pad':4})
 ax.margins(x=.025,y=.045);clean(ax)
 cb=fig.colorbar(im,ax=ax,orientation='horizontal',fraction=.038,pad=.19,shrink=.73)
 cb.set_label('Cells per local bin (log color scale)',fontsize=8);cb.ax.tick_params(labelsize=8,length=2)
def patient(fig,ax,compact=False):
 z=pt[pt.eligible_primary].copy().sort_values(['study','curated_patient_id'],kind='stable')
 yy=np.arange(len(z));vals=z.spearman_rho_f1_vs_scdrs.to_numpy()
 for i,row in enumerate(z.itertuples()):
  ax.plot([0,row.spearman_rho_f1_vs_scdrs],[i,i],color=sd[row.study],lw=1.2,alpha=.55)
  ax.scatter(row.spearman_rho_f1_vs_scdrs,i,s=30,c=sd[row.study],edgecolor='white',linewidth=.5,zorder=3)
 labels=[row.curated_patient_id.split('__',1)[-1] for row in z.itertuples()]
 ax.set_yticks(yy,labels,fontsize=7.5)
 ax.axvline(0,c='#BBBBBB',lw=.7,ls='--')
 ax.errorbar(np.median(rhos),len(z)+1,xerr=[[np.median(rhos)-ci[0]],[ci[1]-np.median(rhos)]],fmt='D',c='#343434',capsize=3,ms=5,lw=1.4)
 ax.text(-.51,len(z)+1,'Median',ha='right',va='center',fontsize=8)
 ax.set(xlim=(-.51,.51),ylim=(len(z)+2,-1.8),xlabel='Within-patient Spearman rho',xticks=[-.4,-.2,0,.2,.4])
 ax.set_title('Patient-level association',loc='left',pad=12);clean(ax)
 ax.text(.0,1.01,'20 eligible patients | 8 positive, 12 negative',transform=ax.transAxes,fontsize=8,color='#666666')
 ax.text(.0,-.15,'Median rho = -0.0490  |  P = 0.4833\n95% patient bootstrap interval: -0.1292 to 0.0636',transform=ax.transAxes,fontsize=8,va='top',linespacing=1.7)
 handles=[Line2D([],[],marker='o',ls='',c=sd[s],label=s,ms=5) for s in studies if s in set(z.study)]
 fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.015 if compact else .02),frameon=False,ncol=len(handles),fontsize=7.5,handletextpad=.3,columnspacing=1.2)
def save(fig,name):
 for ext in ['png','pdf']:fig.savefig(out/(name+'.'+ext),dpi=360,bbox_inches='tight',pad_inches=.16)
 plt.close(fig)
fig,axs=plt.subplots(1,3,figsize=(14.4,5.7),gridspec_kw={'width_ratios':[1.05,1.1,1]})
fig.subplots_adjust(left=.035,right=.982,top=.84,bottom=.26,wspace=.40)
umap(fig,axs[0]);scatter(fig,axs[1]);patient(fig,axs[2],True)
for lab,ax in zip('ABC',axs):ax.text(-.12,1.14,lab,transform=ax.transAxes,fontsize=17,fontweight='bold')
fig.suptitle('F1 usage is not consistently associated with RCC GWAS-locus scDRS scores',fontsize=14,x=.035,ha='left',y=.98)
save(fig,'00_section2_3_scDRS_overview_v1')
fig,ax=plt.subplots(figsize=(5.7,5.9));fig.subplots_adjust(left=.07,right=.96,top=.90,bottom=.10);umap(fig,ax);save(fig,'01_A_scDRS_malignant_UMAP_v1')
fig,ax=plt.subplots(figsize=(6.2,5.8));fig.subplots_adjust(left=.14,right=.96,top=.90,bottom=.20);scatter(fig,ax);save(fig,'02_B_F1_scDRS_density_scatter_v1')
fig,ax=plt.subplots(figsize=(6.2,7.3));fig.subplots_adjust(left=.23,right=.95,top=.91,bottom=.24);patient(fig,ax);save(fig,'03_C_patient_level_association_v1')
versions={'python':sys.version,'executable':sys.executable,'numpy':np.__version__,'pandas':pd.__version__,'matplotlib':matplotlib.__version__}
(q/'SOFTWARE_VERSIONS.json').write_text(json.dumps(versions,indent=2))
audit={'status':'PASS','cells':len(d),'patients_total':int(d.curated_patient_id.nunique()),'studies':int(d.study.nunique()),'eligible_patients':20,'patient_effect_recalculation':'rank-Pearson independent recomputation','median_rho':float(np.median(rhos)),'ci95':ci.tolist(),'signflip_p':float(p),'pooled_pearson_descriptive':pearson,'cell_p_value_calculated':False,'multiple_testing':'single prespecified primary endpoint; no FDR required','cell_fdr_below_0_20':0,'UMAP_color_limits':[low,high],'UMAP_clipping':False,'new_scDRS_or_cNMF_fit':False,'density_bins':[65,65]}
(q/'INDEPENDENT_VALIDATION.json').write_text(json.dumps(audit,indent=2))
more=[src/'local_bundle_v1/PHASE25_SCDRS_ANALYSIS_REGISTRY_v1.md',src/'compute_v2/run_metadata_v1.json',src/'independent_validation_v2/status.json',root/'local_records_migration_20260908_150455/references/image11.png']
manifest=[]
for f in dict.fromkeys(paths+more):
 h=hashlib.sha256()
 with f.open('rb') as rd:
  for chunk in iter(lambda:rd.read(1048576),b''):h.update(chunk)
 manifest.append({'path':str(f),'bytes':f.stat().st_size,'sha256':h.hexdigest()})
(q/'SOURCE_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
pdfchecks=[]
for f in sorted(out.glob('*.pdf')):
 info=subprocess.run(['pdfinfo',str(f)],text=True,capture_output=True,check=True).stdout
 assert any(l.startswith('Pages:') and l.split(':')[1].strip()=='1' for l in info.splitlines())
 subprocess.run(['pdftoppm','-scale-to','1500','-singlefile','-png',str(f),str(q/(f.stem+'_pdf'))],check=True,capture_output=True)
 pdfchecks.append({'file':f.name,'pages':1,'render':'PASS'})
(q/'PDF_VALIDATION.json').write_text(json.dumps(pdfchecks,indent=2))
print(json.dumps(audit),flush=True)
