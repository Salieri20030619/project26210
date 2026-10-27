from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,to_rgba
from matplotlib.patches import FancyBboxPatch,Rectangle
from matplotlib.collections import LineCollection
from scipy.ndimage import gaussian_filter
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1')
d=r/'phase192_sections28_212_figures_20260917_v1';q=r/'part2_rebuild_v1'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':10,'axes.titleweight':'bold','axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
purple='#993F80';teal='#408E91';gold='#BD9753'
cm=LinearSegmentedColormap.from_list('expression',['#759DBF','#FAF9F5',purple])
seq=LinearSegmentedColormap.from_list('density',['#F6F3F7','#BC99BC','#72316C'])
reg=json.loads((d/'inputs/FEATURE_REGISTRY.json').read_text())
m=pd.read_csv(d/'inputs/MACRO_DISPLAY_METADATA.tsv.gz',sep='\t',index_col=0,dtype={'candidate':str})
pal=['#538B9E','#8BA672',purple,'#AA7998','#BD9753','#9B929F','#50978D','#C88773','#7892BE','#B3A877','#77768E','#B684AC']
def save(fig,name):
 for ext in ['pdf','svg','png']:fig.savefig(d/'figures'/f'{name}.{ext}',dpi=600)
 fig.savefig(d/'qc'/f'{name}.jpg',dpi=170)
 plt.close(fig);print('SAVED',name,flush=True)
def title(fig,k,s):fig.suptitle(k+'  |  '+s,x=.035,y=.985,ha='left',fontsize=12,fontweight='bold')
def zrow(a):
 a=np.asarray(a,float);assert np.isfinite(a).all()
 return (a-a.mean(1,keepdims=True))/np.maximum(a.std(1,keepdims=True),1e-8)
def density(x,y,weights=None,bins=180,bounds=None):
 if bounds is None:
  xr=np.ptp(x)*.05;yr=np.ptp(y)*.05;bounds=[[min(x)-xr,max(x)+xr],[min(y)-yr,max(y)+yr]]
 h,xe,ye=np.histogram2d(x,y,bins=bins,range=bounds,weights=weights)
 return gaussian_filter(h.T,2),[(xe[0]+xe[1])/2,(xe[-2]+xe[-1])/2,(ye[0]+ye[1])/2,(ye[-2]+ye[-1])/2],bounds
x=m.UMAP1.values;y=m.UMAP2.values
base,ext,bounds=density(x,y)
def outline(ax,xx=x,yy=y,mask=None,b=bounds,color=purple):
 mask=(m.candidate=='4').values if mask is None else mask
 if mask.sum()<5:return
 h,e,_=density(xx[mask],yy[mask],bounds=b)
 ax.contour(h,levels=[h.max()*.12],extent=e,colors=[color],linewidths=.8,zorder=4)
def clean(ax):ax.set_xticks([]);ax.set_yticks([]);ax.set_aspect('equal');ax.spines[['top','right','left','bottom']].set_visible(False)
# 01 actual biochemical relations only
fig,ax=plt.subplots(figsize=(7.48,5.0));ax.set(xlim=(0,10),ylim=(0,7));ax.axis('off')
title(fig,'2.8','A biochemical route to MTA')
ax.text(5,6.6,'TYMP-associated metabolic features',ha='center',weight='bold',color=purple,fontsize=11)
ax.text(5,6.22,'Cysteine / methionine     ·     Arginine / proline     ·     ODC1 / SRM',ha='center',fontsize=8.5)
def box(x,y,t,c='#EDF3F2',w=1.65):
 ax.add_patch(FancyBboxPatch((x-w/2,y-.28),w,.56,boxstyle='round,pad=.07,rounding_size=.08',fc=c,ec='none'))
 ax.text(x,y,t,ha='center',va='center',fontsize=8.5)
def arr(a,b,t='',offset=(0,0)):
 ax.annotate('',xy=b,xytext=a,arrowprops=dict(arrowstyle='-|>',lw=1.25,color='#596568'))
 if t:ax.text((a[0]+b[0])/2+offset[0],(a[1]+b[1])/2+offset[1],t,ha='center',va='center',fontsize=8,color=teal)
ax.text(.15,5.45,'SAM supply',weight='bold',color=teal)
box(1.4,4.65,'Methionine');box(4,4.65,'SAM');box(7,4.65,'dcSAM')
arr((2.3,4.65),(3.1,4.65),'MAT2A',(0,.32));arr((4.9,4.65),(6.1,4.65),'AMD1',(0,.32))
ax.text(.15,3.6,'Polyamine synthesis',weight='bold',color=teal)
for xx,tt in [(1.2,'Ornithine'),(3.65,'Putrescine'),(6.1,'Spermidine'),(8.6,'Spermine')]:box(xx,2.95,tt,w=1.65)
arr((2.1,2.95),(2.77,2.95),'ODC1',(0,.32))
arr((4.53,2.95),(5.22,2.95),'SRM',(0,.32));arr((6.99,2.95),(7.72,2.95),'SMS',(0,.32))
ax.text(4.88,3.58,'+ dcSAM',ha='center',fontsize=7,color='#777777');ax.text(7.36,3.58,'+ dcSAM',ha='center',fontsize=7,color='#777777')
box(6.1,1.52,"5'-MTA",c='#ECDDE9',w=2.1)
arr((4.88,2.75),(5.63,1.84));arr((7.36,2.75),(6.57,1.84))
ax.text(6.1,2.17,'reaction product',ha='center',fontsize=7,color=purple)
box(2.3,1.52,'MTR-1-P + adenine',w=2.65)
arr((4.97,1.52),(3.7,1.52),'MTAP',(0,.3))
ax.plot([.93,.05,.05,1.4],[1.52,1.52,4.1,4.1],ls='--',lw=1,color='#999999')
ax.annotate('',xy=(1.4,4.28),xytext=(1.4,4.1),arrowprops=dict(arrowstyle='->',ls='--',lw=1,color='#999999'))
ax.text(.18,1.95,'Methionine\nsalvage',fontsize=7,color='#777777')
ax.text(5,.45,'Solid arrows: biochemical reactions     Dashed arrow: multistep salvage',ha='center',fontsize=7,color='#666666')
fig.subplots_adjust(top=.9,bottom=.03,left=.025,right=.99);save(fig,'01_MTA_metabolic_map')
# 02
fig=plt.figure(figsize=(7.48,5.2));title(fig,'2.9','A distinct SPP1-enriched macrophage state')
gs=fig.add_gridspec(2,2,height_ratios=[3.2,1],top=.87,bottom=.06,wspace=.12)
for j in [0,1]:
 ax=fig.add_subplot(gs[0,j]);ax.contourf(base,levels=[base.max()*.006,base.max()],extent=ext,colors=['#E7E8E8'])
 if j==0:
  for c,col in zip(reg['clusters'],pal):
   mask=(m.candidate==c).values;h,e,_=density(x[mask],y[mask],bounds=bounds)
   ax.contourf(h,levels=[h.max()*.16,h.max()*.48,h.max()+1e-5],extent=e,colors=[to_rgba(col,.40),to_rgba(col,.85)])
   # Label at density mode instead of outlier-prone mean.
   iy,ix=np.unravel_index(h.argmax(),h.shape)
   ax.text(e[0]+ix/(h.shape[1]-1)*(e[1]-e[0]),e[2]+iy/(h.shape[0]-1)*(e[3]-e[2]),'C'+c,fontsize=7,weight='bold',ha='center',bbox=dict(fc='white',ec='none',alpha=.8,pad=.6))
 else:
  mask=(m.candidate=='4').values;h,e,_=density(x[mask],y[mask],bounds=bounds)
  ax.contourf(h,levels=np.linspace(h.max()*.03,h.max()+1e-5,8),extent=e,cmap=seq);outline(ax)
 ax.set_title(['A   Frozen macrophage states','B   C4 localization'][j],loc='left',fontsize=9);clean(ax)
ax=fig.add_subplot(gs[1,:]);ax.axis('off')
for k,(c,col) in enumerate(zip(reg['clusters'],pal)):
 xx=(k%3)/3;yy=.94-(k//3)*.25
 ax.add_patch(Rectangle((xx,yy-.07),.017,.08,transform=ax.transAxes,fc=col,ec='none'))
 ax.text(xx+.025,yy-.03,'C'+c+'  '+reg['labels'][c],transform=ax.transAxes,va='center',fontsize=6.8)
save(fig,'02_macrophage_density_atlas')
# heatmap helper
def heat(df,name,section,heading,figheight,groups=None):
 zz=zrow(df.values);fig=plt.figure(figsize=(7.48,figheight));title(fig,section,heading)
 left=.27 if len(df)<15 else .18
 ax=fig.add_axes([left,.12,.64,.76]);im=ax.imshow(zz,cmap=cm,vmin=-2.5,vmax=2.5,aspect='auto',interpolation='nearest')
 ax.set_xticks(range(len(df.columns)),df.columns,fontsize=8)
 ax.set_yticks(range(len(df)),df.index,fontsize=7)
 if len(df)>15:
  for t in ax.get_yticklabels():t.set_fontstyle('italic')
 ax.tick_params(length=0);c4=list(df.columns).index('C4')
 ax.add_patch(Rectangle((c4-.5,-.5),1,len(df),fill=False,ec=purple,lw=1.6,clip_on=False))
 if groups:
  pos=0
  for i,(g,genes) in enumerate(groups.items()):
   ax.axhline(pos-.5,color='white',lw=1.2)
   pos+=len(genes)
 ax.set_xlabel('Frozen macrophage candidate',labelpad=9)
 ca=fig.add_axes([.865 if len(df)>15 else .915,.31,.015,.25]);fig.colorbar(im,cax=ca,ticks=[-2,0,2]);ca.set_title('Row z',fontsize=7,pad=6)
 save(fig,name)
hm=pd.read_csv(d/'tables/MARKER_PATIENT_EQUAL_LOG2CPM.tsv',sep='\t',index_col=0)
heat(hm,'03_macrophage_marker_heatmap','2.9','Marker programs resolve C4 identity',7.0,reg['marker_groups'])
# 04
bulk=pd.read_csv(d/'tables/TISSUE_MTA_ORDERED_PATIENTS.tsv',sep='\t',index_col=0)
features=['MET:5-methylthioadenosine (MTA)','SPP1','LST1','AIF1','CSF1R','CD68','FCER1G','SMS','MTAP']
fig=plt.figure(figsize=(7.48,4.8));title(fig,'2.10','Matched MTA and tissue expression landscapes')
left=.13;width=.82;start=0
for i,(co,dat) in enumerate(bulk.groupby('cohort',sort=False)):
 w=width*len(dat)/255;ax=fig.add_axes([left+width*start/255,.23,w-.009,.59])
 v=dat[features].values.T.copy();v[0]=np.log2(v[0]+1);im=ax.imshow(zrow(v),cmap=cm,vmin=-2.5,vmax=2.5,aspect='auto',interpolation='nearest')
 ax.set_title(co.replace('_',' ')+'\nn = '+str(len(dat)),fontsize=8)
 ax.set_yticks(range(len(features)),['MTA']+features[1:] if i==0 else ['']*len(features));ax.set_xticks([]);ax.tick_params(length=0)
 ax.add_patch(Rectangle((-.5,.5),len(dat),1,fill=False,ec=purple,lw=1.1))
 ax.annotate('',xy=(.97,-.035),xytext=(.03,-.035),xycoords='axes fraction',arrowprops=dict(arrowstyle='->',lw=.8,color='#666666'))
 start+=len(dat)
fig.text(.54,.15,'Patients ordered by measured MTA within each cohort',ha='center',fontsize=8)
ca=fig.add_axes([.34,.08,.36,.02]);fig.colorbar(im,cax=ca,orientation='horizontal',label='Within-cohort row z-score',ticks=[-2,0,2])
save(fig,'04_MTA_SPP1_tissue_landscape')
pm=pd.read_csv(d/'tables/PROGRAM_PATIENT_EQUAL_SCORES.tsv',sep='\t',index_col=0)
pm.index=[s.replace(' processing and presentation',' presentation').replace('Interferon Gamma Response','IFN-gamma response').replace('Cholesterol Homeostasis','Cholesterol homeostasis') for s in pm.index]
heat(pm,'05_macrophage_program_heatmap','2.11','Lipid and lysosomal programs characterize C4',4.7)
# 06 maps actual scores, no quantile distortions
fig,axs=plt.subplots(2,2,figsize=(7.48,6.25));title(fig,'2.11','Functional expression across macrophage states')
fields=['expr:SPP1','score:Cholesterol metabolism','score:Lysosome','score:Antigen processing and presentation']
for ax,f,letter in zip(axs.ravel(),fields,'ABCD'):
 vals=m[f].values;h,e,_=density(x,y,weights=vals,bounds=bounds)
 val=np.divide(h,base,out=np.zeros_like(h),where=base>0)
 val=np.ma.masked_where(base<base.max()*.008,val)
 im=ax.imshow(val,extent=e,origin='lower',cmap=seq,interpolation='bilinear',vmin=0,vmax=np.percentile(val.compressed(),99))
 outline(ax);clean(ax)
 lab=f.split(':')[1].replace('Antigen processing and presentation','Antigen presentation')
 ax.set_title(letter+'   '+lab,loc='left',fontsize=9)
 cb=fig.colorbar(im,ax=ax,fraction=.032,pad=.015,shrink=.7);cb.ax.tick_params(labelsize=6)
 cb.ax.set_title('Mean\nlog1p',fontsize=6,pad=5)
fig.subplots_adjust(top=.9,left=.025,right=.94,bottom=.05,wspace=.15,hspace=.18)
save(fig,'06_C4_expression_landscapes')
# 07 original graph and actual partition membership validation
pt=pd.read_csv(q/'07_trajectory/monocle3_v3/ALL_CELL_PSEUDOTIME.tsv.gz',sep='\t',index_col=0,dtype={'partition':str,'candidate':str},low_memory=False)
roots=pd.read_csv(q/'07_trajectory/monocle3_v3/ROOTS.tsv',sep='\t',dtype={'partition':str})
fig,axs=plt.subplots(2,3,figsize=(7.48,5.4));title(fig,'2.12','Original partitioned myeloid trajectories')
audit=[]
for ax,study in zip(axs.ravel(),sorted(pt.study.unique())):
 p=d/'inputs/trajectory/primary'/study
 cc=pd.read_csv(p/'cells.tsv',sep='\t',index_col=0);nn=pd.read_csv(p/'nodes.tsv',sep='\t',index_col=0);ee=pd.read_csv(p/'edges.tsv',sep='\t');cl=pd.read_csv(p/'closest.tsv',sep='\t',index_col=0)
 pp=pt[(pt['mode']=='primary')&(pt.study==study)].loc[cc.index]
 edgeidx=nn.index.get_indexer(ee['from']);edgej=nn.index.get_indexer(ee['to']);assert min(edgeidx.min(),edgej.min())>=0
 graph=coo_matrix((np.ones(len(ee)),(edgeidx,edgej)),shape=(len(nn),len(nn)))
 ncomp,comp=connected_components(graph,directed=False)
 mapping=pd.DataFrame({'component':comp[nn.index.get_indexer(cl.loc[cc.index,'node'])],'partition':pp.partition.values})
 assert mapping.groupby('component').partition.nunique().max()==1
 xx=cc.x.values;yy=cc.y.values;hh,ex,bo=density(xx,yy)
 ax.contourf(hh,levels=[hh.max()*.006,hh.max()+1e-5],extent=ex,colors=['#E0E2E3'])
 ok=np.isfinite(pp.pt_primary.values);vv=pp.pt_primary.values.copy()
 if ok.any():
  # Normalization is within original graph partition, never across disconnected graphs.
  for part in pp.partition.unique():
   ii=ok&(pp.partition.values==part)
   if ii.any():vv[ii]=(vv[ii]-vv[ii].min())/max(np.ptp(vv[ii]),1e-9)
  den,_,_=density(xx[ok],yy[ok],bounds=bo);num,_,_=density(xx[ok],yy[ok],weights=vv[ok],bounds=bo)
  v=np.ma.masked_where(den<max(den.max()*.018,1e-9),num/np.maximum(den,1e-9))
  im=ax.imshow(v,extent=ex,origin='lower',cmap='viridis',vmin=0,vmax=1,zorder=2)
 seg=np.stack([nn.loc[ee['from'],['x','y']].values,nn.loc[ee['to'],['x','y']].values],axis=1)
 ax.add_collection(LineCollection(seg,colors='#343838',linewidths=.45,alpha=.85,zorder=3))
 outline(ax,xx,yy,(pp.candidate=='4').values,bo)
 rr=roots[(roots['mode']=='primary')&(roots.study==study)&(roots.root_version=='primary')]
 for rc in rr.root_cell:
  if rc in cc.index:ax.plot(cc.loc[rc,'x'],cc.loc[rc,'y'],marker='*',ms=5,mfc='white',mec='black',mew=.6)
 clean(ax);ax.set_title(study,loc='left',fontsize=9)
 audit.append({'study':study,'cells':len(cc),'nodes':len(nn),'edges':len(ee),'graph_components':ncomp,'cross_partition_components':0})
ax=axs.ravel()[-1];ax.axis('off')
ax.text(.03,.87,'Original graph edges',fontsize=8,transform=ax.transAxes)
ax.text(.03,.67,'Purple contour: C4',color=purple,fontsize=8,transform=ax.transAxes)
ax.text(.03,.47,'Star: monocyte root',fontsize=8,transform=ax.transAxes)
ax.text(.03,.27,'Gray: no finite pseudotime',color='#777777',fontsize=8,transform=ax.transAxes)
fig.subplots_adjust(left=.035,right=.965,bottom=.12,top=.88,hspace=.25,wspace=.13)
ca=fig.add_axes([.34,.072,.32,.018]);fig.colorbar(im,cax=ca,orientation='horizontal',label='Within-partition normalized pseudotime',ticks=[0,.5,1])
pd.DataFrame(audit).to_csv(d/'qc/ORIGINAL_GRAPH_VALIDATION.tsv',sep='\t',index=False)
save(fig,'07_partitioned_trajectory_maps')
# 08 independent partition then patient equal smoothing; joint row scaling.
dy=np.load(d/'inputs/DYNAMIC_DISPLAY.npz')
a=dy['primary'];b=dy['sensitivity'];joint=zrow(np.c_[a,b]);fig=plt.figure(figsize=(7.48,7.1));title(fig,'2.12','Expression dynamics along rooted myeloid states')
for j,(z,label) in enumerate([(joint[:,:40],'Primary | 14 patients'),(joint[:,40:],'Sensitivity | 24 patients')]):
 ax=fig.add_axes([.19+j*.37,.12,.33,.70]);im=ax.imshow(z,aspect='auto',cmap=cm,vmin=-2.5,vmax=2.5,interpolation='nearest')
 ax.set_yticks(range(len(dy['genes'])),dy['genes'] if j==0 else ['']*len(dy['genes']),fontsize=7)
 for t in ax.get_yticklabels():t.set_fontstyle('italic')
 ax.set_xticks([0,19.5,39],['Early','Intermediate','Late']);ax.tick_params(length=0);ax.set_xlabel('Within-partition pseudotime rank')
 pos=0
 for k,(group,genes) in enumerate(reg['dynamic_groups'].items()):
  ax.axhline(pos-.5,color='white',lw=1)
  ax.add_patch(Rectangle((-.95,pos-.5),.3,len(genes),fc=pal[k],ec='none',clip_on=False));pos+=len(genes)
 ax.set_title(label,fontsize=9,pad=34)
 cax=fig.add_axes([.19+j*.37,.835,.33,.035]);ci=cax.imshow(dy['composition_primary' if j==0 else 'composition_sensitivity'],aspect='auto',cmap='Greys',vmin=0,vmax=1)
 cax.set_xticks([]);cax.set_yticks([0,1],['C4','Mono'] if j==0 else ['',''],fontsize=6);cax.tick_params(length=0)
ca=fig.add_axes([.925,.36,.016,.25]);fig.colorbar(im,cax=ca,label='Joint row z-score',ticks=[-2,0,2])
fig.text(.56,.91,'Composition strips: patient-mean fraction (white 0; black 1)',ha='center',fontsize=6.5,color='#666666')
save(fig,'08_pseudotime_expression_heatmap')
(d/'qc/RENDER_STATUS.json').write_text(json.dumps({'figures':8,'formats':['pdf','svg','png'],'png_dpi':600,'new_tests':False,'raw_scatter_or_boxplots':False,'graph_checks':audit},indent=2))
