from pathlib import Path
import json,hashlib,sys,platform,subprocess
import numpy as np,pandas as pd,h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize,TwoSlopeNorm,LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from PIL import Image
from scipy.spatial import cKDTree
P=Path(__file__).resolve().parents[1];R=P.parent
A=R/'phase234_part4_replan_execution_20260921_v1'
B=R/'phase250_independent_sixpatient_mainline_20260924_v1'
C=R/'phase248_independent_cases_mainline_20260924_v1'
D=R/'phase246_spatial_marker_structure_20260924_v1'
sources=[]
def read(f):
 f=Path(f);sources.append(f);return pd.read_csv(f,sep='\t')
su=read(A/'tables/A1_COHORT_SUMMARY.tsv')
su=su[(su.layer=='broad')&(su.scale=='W0')&(su.model=='A0')&(su.qualification=='Q0')].copy()
reg=read(A/'tables/INPUT_SECTION_REGISTER.tsv')
pat=read(B/'tables/PATIENT_EFFECTS.tsv')
ext=read(C/'tables/EXTERNAL_CASE_MAINLINE_EFFECTS.tsv')
paired=read(B/'tables/PAIRED_TUMOR_ADJACENT_DESCRIPTIVE.tsv')
dist=read(D/'tables/MARKER_STRUCTURE_PATIENT_EFFECTS.tsv')
distsu=read(D/'tables/MARKER_STRUCTURE_COHORT_SUMMARY_WITH_EFFECTIVE_N.tsv')
data={};imgs={};scales={};lims={}
for sample in ['GLMF1','GLMF2']:
 d=read(A/f'inputs/{sample}_INPUTS.tsv.gz');data[sample]=d
 row=reg[reg['sample']==sample].iloc[0]
 f=Path(row.input_path);sources.append(f)
 with h5py.File(f,'r') as h:
  grp=h['uns/spatial'][sample];imgs[sample]=grp['images/hires'][()]
  sf=grp['scalefactors/tissue_hires_scalef'][()]
  diam=grp['scalefactors/spot_diameter_fullres'][()]*sf
 # Whole section crop, one square viewport shared by every channel.
 xmin,xmax=d.x_hires_px.min(),d.x_hires_px.max()
 ymin,ymax=d.y_hires_px.min(),d.y_hires_px.max()
 width=max(xmax-xmin,ymax-ymin)+130
 cx=(xmin+xmax)/2;cy=(ymin+ymax)/2
 lims[sample]=(cx-width/2,cx+width/2,cy-width/2,cy+width/2)
 scales[sample]=dict(um_per_hires_px=float(row.pixel_to_um/sf),spot_diameter_hires=float(diam))
 assert imgs[sample].shape[2]==3
 assert np.isfinite(d[['x_hires_px','y_hires_px','TYMP_RNA','polyamine_mta_context','V1_SPP1_q05']]).all().all()
 d[['spot_id','x_hires_px','y_hires_px','TYMP_RNA','polyamine_mta_context','V1_SPP1_q05']].to_csv(P/f'plot_data/{sample}_mapped_spots.tsv',sep='\t',index=False)
 # Coordinates independently agree with stored physical units.
 xy=d[['x_hires_px','y_hires_px']].to_numpy();uv=d[['physical_x_um','physical_y_um']].to_numpy()
 distance,ix=cKDTree(xy).query(xy,k=2)
 ratio=np.median(np.linalg.norm(uv-uv[ix[:,1]],axis=1)/distance[:,1])
 scales[sample]['nearest_neighbor_scale_check']=float(ratio)
 assert abs(ratio/scales[sample]['um_per_hires_px']-1)<.01
cohorts=['GSE175540','TJRCC'];colors={'GSE175540':'#307F9B','TJRCC':'#CE7659'}
labels={'GSE175540':'GSE175540','TJRCC':'TJ-RCC'}
targetsE=['V1_SPP1_q05','V2_SPP1_fraction','V3_SPP1_macro_share','V4_SPP1_program']
targetsF=['invasion_emt','immunosuppression_core','tcell_dysfunction','cytotoxicity']
targetsG=['endothelial_q05','pericyte_q05','endothelial_fraction','pericyte_fraction']
e=su[(su.family=='MTA_context_SPP1')&(su.target.isin(targetsE))&(su.dataset.isin(cohorts))].copy()
f=su[(su.family=='SPP1_function')&(su.exposure=='V4_SPP1_program')&(su.target.isin(targetsF))&(su.dataset.isin(cohorts))].copy()
g=su[(su.family=='SPP1_vascular')&(su.exposure=='V1_SPP1_q05')&(su.target.isin(targetsG))&(su.dataset.isin(cohorts))].copy()
assert len(e)==len(f)==len(g)==8
pairs=[('TYMP_RNA','polyamine_mta_context'),('TYMP_RNA','spp1_state_no_spp1'),('TYMP_RNA','SPP1_RNA'),('polyamine_mta_context','spp1_state_no_spp1'),('polyamine_mta_context','SPP1_RNA')]
pt=pat[(pat.tissue=='tumor')&(pat.role=='mainline')&(pat.scale=='block4x4')&(pat.model=='A0')].copy();pt['dataset']='E-MTAB-12767'
ex=ext[(ext.role=='mainline')&(ext.scale=='block4x4')&(ext.model=='A0')].copy()
rows=sorted(pt.patient.unique().tolist())+ex[['patient','sample']].drop_duplicates().sort_values('sample').patient.tolist()
assert len(rows)==9,len(rows)
allh=pd.concat([pt,ex],ignore_index=True);hm=np.full((9,5),np.nan)
for i,patient in enumerate(rows):
 for j,(x,y) in enumerate(pairs):
  rr=allh[(allh.patient==patient)&(allh.exposure==x)&(allh.target==y)]
  assert len(rr)==1,(patient,x,y,len(rr))
  hm[i,j]=rr.effect.iloc[0]
assert (hm[:6,3]>0).all() and (hm[6:,4]>0).all()
ep='SPP1_high_minus_low_distance_to_CD8lineage'
jdata=dist[(dist.endpoint==ep)&(dist.model=='fixed_20percent')].copy()
js=distsu[(distsu.endpoint==ep)&(distsu.model=='fixed_20percent')].copy()
assert int(js[js.dataset=='GSE175540'].n.iloc[0])==16
for key,tab in [('E_context_SPP1',e),('F_functional_programs',f),('G_vascular',g),('H_independent_patients',allh),('I_paired_tissues',paired),('J_CD8_patients',jdata),('J_CD8_summary',js)]:
 tab.to_csv(P/f'plot_data/{key}.tsv',sep='\t',index=False)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'axes.titlesize':7.5,'axes.labelsize':7,'xtick.labelsize':6.5,'ytick.labelsize':6.5,'axes.linewidth':.55,'xtick.major.width':.5,'ytick.major.width':.5,'xtick.major.size':2.5,'ytick.major.size':2.5,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
fig=plt.figure(figsize=(210/25.4,297/25.4),facecolor='white')
def axmm(x,y,w,h):return fig.add_axes([x/210,y/297,w/210,h/297])
def textmm(x,y,s,**kw):return fig.text(x/210,y/297,s,**kw)
def heading(letter,title,x,y):
 textmm(x,y,letter,fontweight='bold',fontsize=10,va='center')
 textmm(x+5,y,title,fontsize=7.7,fontweight='bold',va='center')
def clean(ax):
 ax.spines[['top','right']].set_visible(False)
 ax.tick_params(pad=2)
def mark(row):
 return '*' if row.FDR<.05 else ('†' if row.P<.05 else '')
# A-D: complete tissue fields, same crop and scale per section.
heads=[('A','Tissue morphology'),('B','TYMP expression'),('C','Metabolic context'),('D','SPP1 macrophages')]
xx=[8,58,108,158];maps=[]
tympcm=LinearSegmentedColormap.from_list('tymp',['#f4eff6','#c8afd8','#8452a1','#3d155f'])
macm=LinearSegmentedColormap.from_list('mac',['#f0f6f4','#98c6b9','#39927e','#005443'])
cmaps=[None,tympcm,'RdBu_r',macm]
norms=[None,Normalize(0,2),Normalize(-1.25,1.25),Normalize(0,3)]
variables=[None,'TYMP_RNA','polyamine_mta_context','V1_SPP1_q05']
for col,(letter,title) in enumerate(heads):
 heading(letter,title,xx[col],286)
 for row,sample in enumerate(['GLMF1','GLMF2']):
  ay=237-row*45
  ax=axmm(xx[col]+1,ay,42,42);maps.append(ax)
  x0,x1,y0,y1=lims[sample]
  if col==0:ax.imshow(imgs[sample],interpolation='bilinear')
  else:
   # very faint tissue background makes spot colors and geometry legible
   ax.imshow(imgs[sample],alpha=.16,interpolation='bilinear')
   d=data[sample]
   patches=[Circle((x,y),scales[sample]['spot_diameter_hires']/2) for x,y in zip(d.x_hires_px,d.y_hires_px)]
   coll=PatchCollection(patches,cmap=cmaps[col],norm=norms[col],edgecolor='none',rasterized=True)
   coll.set_array(d[variables[col]].to_numpy());ax.add_collection(coll)
  ax.set_xlim(x0,x1);ax.set_ylim(y1,y0);ax.set_aspect('equal');ax.axis('off')
  if col==0:
   ax.text(.01,1.018,sample,transform=ax.transAxes,fontsize=7,fontweight='bold')
   length=500/scales[sample]['um_per_hires_px'];sx=x0+(x1-x0)*.055;sy=y1-(y1-y0)*.06
   ax.plot([sx,sx+length],[sy,sy],lw=1.7,color='#111111')
   ax.text(sx,sy-(y1-y0)*.028,'500 μm',fontsize=6,va='bottom',color='#111111',bbox=dict(facecolor='white',edgecolor='none',alpha=.75,pad=.4))
 if col:
  cax=axmm(xx[col]+8,186,29,2)
  cb=fig.colorbar(plt.cm.ScalarMappable(norm=norms[col],cmap=cmaps[col]),cax=cax,orientation='horizontal',extend='both' if col==2 else 'max')
  cb.outline.set_linewidth(.4);cb.ax.tick_params(labelsize=6,length=1.8,pad=1)
  cb.set_ticks([0,1,2] if col==1 else ([-1,0,1] if col==2 else [0,1.5,3]))
  textmm(xx[col]+22.5,179.5,['','log1p(CP10k)','Mean gene z-score','Estimated abundance'][col],ha='center',fontsize=6.3)
# E median effect and CI
heading('E','Metabolic context and SPP1',8,170)
legend=[Line2D([0],[0],marker='o',color=colors[c],lw=1,ms=3,label=labels[c]+' (n='+str(23 if c=='GSE175540' else 12)+')') for c in cohorts]
fig.legend(handles=legend,loc='upper left',bbox_to_anchor=(8/210,166/297),frameon=False,fontsize=6.2,borderaxespad=0,handlelength=1.2,labelspacing=.25)
ae=axmm(34,122,31,34)
for k,targ in enumerate(targetsE):
 for ic,c in enumerate(cohorts):
  rr=e[(e.dataset==c)&(e.target==targ)].iloc[0];yy=k+(ic-.5)*.27
  ae.plot([rr.CI_low,rr.CI_high],[yy,yy],color=colors[c],lw=1.1,zorder=2)
  ae.plot(rr['median'],yy,'o',ms=3.3,color=colors[c],zorder=3)
  ae.text(.315,yy,mark(rr),va='center',ha='left',fontsize=7,color=colors[c])
ae.axvline(0,color='#999999',ls='--',lw=.6,zorder=0)
ae.set_yticks(range(4),['Abundance','All-cell fraction','Macrophage share','State score'])
ae.set_ylim(3.55,-.65);ae.set_xlim(-.16,.35);ae.set_xticks([-.1,0,.1,.2,.3]);ae.set_xlabel('Patient median ρ (95% CI)',labelpad=4);clean(ae)
# F all four prespecified functional endpoints
heading('F','SPP1 state and function',76,170)
af=axmm(108,122,24,34)
fmat=np.array([[f[(f.dataset==c)&(f.target==t)].iloc[0]['median'] for c in cohorts] for t in targetsF])
cf=af.imshow(fmat,cmap='RdBu_r',vmin=-.2,vmax=.2,aspect='auto')
af.set_yticks(range(4),['Invasion / EMT','Myeloid suppression','T-cell dysfunction','Cytotoxicity'])
af.set_xticks([0,1],['GSE175540','TJ-RCC'],rotation=35,ha='right')
af.tick_params(length=0,pad=3)
for i,targ in enumerate(targetsF):
 for k,c in enumerate(cohorts):
  rr=f[(f.dataset==c)&(f.target==targ)].iloc[0]
  af.text(k,i,f'{rr["median"]:.2f}'+mark(rr),ha='center',va='center',fontsize=6.8,color='white' if abs(rr['median'])>.14 else '#202020')
for sp in af.spines.values():sp.set_visible(False)
cax=axmm(111,111,20,1.8);cb=fig.colorbar(cf,cax=cax,orientation='horizontal',ticks=[-.2,0,.2]);cb.outline.set_linewidth(.4);cb.ax.tick_params(labelsize=6,length=1.5,pad=1)
textmm(108,105.8,'Patient median ρ',fontsize=6.2)
# G matched absolute versus compositional vascular readouts
heading('G','Vascular co-distribution',141,170)
ag=axmm(171,122,29,34)
for k,targ in enumerate(targetsG):
 for ic,c in enumerate(cohorts):
  rr=g[(g.dataset==c)&(g.target==targ)].iloc[0];yy=k+(ic-.5)*.27
  ag.plot([rr.CI_low,rr.CI_high],[yy,yy],color=colors[c],lw=1.1)
  ag.plot(rr['median'],yy,'o',ms=3.3,color=colors[c])
  ag.text(.55,yy,mark(rr),va='center',ha='left',fontsize=7,color=colors[c])
ag.axvline(0,color='#999999',ls='--',lw=.6)
ag.set_yticks(range(4),['Endothelium: amount','Pericytes: amount','Endothelium: fraction','Pericytes: fraction'])
ag.set_ylim(3.55,-.65);ag.set_xlim(-.42,.6);ag.set_xticks([-.4,0,.4]);ag.set_xlabel('Patient median ρ (95% CI)',labelpad=4);clean(ag)
textmm(143,160,'SPP1 macrophage abundance vs vascular readout',fontsize=6)
# H independent case matrix
heading('H','Independent spatial cases',8,99)
ah=axmm(25,42,53,46)
cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#e4e4e4')
ih=ah.imshow(np.ma.masked_invalid(hm),vmin=-.8,vmax=.8,cmap=cmap,aspect='auto')
rowlabels=[x.split(':')[-1] for x in rows[:6]]+[ex[ex.patient==x]['sample'].iloc[0] for x in rows[6:]]
ah.set_yticks(range(9),rowlabels);ah.set_xticks(range(5),['TYMP\n~ context','TYMP\n~ state','TYMP\n~ SPP1','Context\n~ state','Context\n~ SPP1'],fontsize=6.1)
ah.tick_params(length=0,pad=3)
ah.axhline(5.5,color='white',lw=2)
for i in range(9):
 for k in range(5):
  value=hm[i,k]
  ah.text(k,i,'NA' if not np.isfinite(value) else f'{value:.2f}',ha='center',va='center',fontsize=6.2,color='white' if np.isfinite(value) and abs(value)>.52 else '#222222')
for sp in ah.spines.values():sp.set_visible(False)
textmm(9,72,'E-MTAB-12767',rotation=90,fontsize=6.1,ha='center',va='center')
textmm(9,49,'GSE294006',rotation=90,fontsize=6.1,ha='center',va='center')
textmm(25,91,'Unadjusted 4 × 4 regions; one row per patient',fontsize=6.2)
cax=axmm(40,28,25,1.8);cb=fig.colorbar(ih,cax=cax,orientation='horizontal',ticks=[-.8,0,.8]);cb.outline.set_linewidth(.4);cb.ax.tick_params(labelsize=6,length=1.5,pad=1)
textmm(26,20,'Spearman ρ     Gray: not estimable',fontsize=6.1)
# I matched adjacent and tumor, observed signals and no inferential P for 2 pairs
heading('I','Tumor–adjacent pairs',86,99)
pcolors={'EMTAB12767:HG1':'#9468A6','EMTAB12767:HG3':'#D69549'}
textmm(90,92,'HG1',color=pcolors['EMTAB12767:HG1'],fontsize=6.7,fontweight='bold')
textmm(105,92,'HG3',color=pcolors['EMTAB12767:HG3'],fontsize=6.7,fontweight='bold')
ivars=['TYMP_RNA','SPP1_RNA','polyamine_mta_context','spp1_state_no_spp1']
ititles=['TYMP','SPP1','Metabolic context','SPP1 state score']
for k,(var,title) in enumerate(zip(ivars,ititles)):
 x=91+(k%2)*28;y=62-(k//2)*34
 ai=axmm(x,y,22,22)
 for patient,col in pcolors.items():
  rr=paired[(paired.patient==patient)&(paired.variable==var)].iloc[0]
  ai.plot([0,1],[rr.adjacent,rr.tumor],'-o',color=col,lw=.9,ms=3)
 ai.set_xticks([0,1],['Adj.','Tumor'],fontsize=6)
 ai.set_xlim(-.3,1.3);ai.set_title(title,fontsize=6.5,pad=3)
 ai.locator_params(axis='y',nbins=3);ai.tick_params(axis='y',labelsize=6);clean(ai)
textmm(116,18,'Mean log1p(CP10k); n = 2 pairs',fontsize=6.2,ha='center')
# J supplemental regional distance, all three sources including nonreplication
heading('J','CD8-region distance',150,99)
aj=axmm(164,37,37,46)
jco=['GSE175540','TJRCC','GLMF_GSE210041'];jcolors=[colors['GSE175540'],colors['TJRCC'],'#888888']
rng=np.random.default_rng(262254)
for k,c in enumerate(jco):
 vals=jdata[jdata.dataset==c].effect.dropna().to_numpy();assert len(vals)>0
 aj.scatter(k+rng.uniform(-.15,.15,len(vals)),vals,s=13,color=jcolors[k],alpha=.88,edgecolors='white',linewidths=.25,zorder=3)
 med=np.median(vals);aj.plot([k-.25,k+.25],[med,med],color='#111111',lw=1.3,zorder=4)
 aj.text(k,1.02,'n = '+str(len(vals)),ha='center',va='bottom',transform=aj.get_xaxis_transform(),fontsize=6.3)
aj.axhline(0,color='#999999',ls='--',lw=.6)
aj.set_xticks(range(3),['GSE175540','TJ-RCC','GLMF'],rotation=35,ha='right',fontsize=6.3)
aj.set_ylabel('High − low SPP1 distance (μm)',fontsize=6.5,labelpad=3)
aj.set_xlim(-.5,2.5);clean(aj)
textmm(153,91,'GSE175540: P = 0.035; q = 0.053',fontsize=6.1)
textmm(150,21,'Bars: median; positive = farther from CD8 regions',fontsize=6)
# Small shared symbol key only; explanatory caption is separate.
textmm(8,8,'* FDR < 0.05     † P < 0.05 only     State score excludes SPP1     All n values denote independent patients',fontsize=6.4,color='#333333')
assert e.CI_low.min()>=ae.get_xlim()[0] and e.CI_high.max()<=ae.get_xlim()[1]
assert g.CI_low.min()>=ag.get_xlim()[0] and g.CI_high.max()<=ag.get_xlim()[1]
base=P/'delivery/26210_Spatial_A4_10panels_v1'
fig.savefig(str(base)+'.pdf',dpi=450)
fig.savefig(str(base)+'.svg',dpi=450)
fig.savefig(str(base)+'.png',dpi=300)
fig.savefig(str(base)+'.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'})
fig.canvas.draw()
# All text must remain inside the canvas and at least 6 pt.
renderer=fig.canvas.get_renderer();bad=[];small=[]
for obj in fig.findobj(matplotlib.text.Text):
 if not obj.get_visible() or not obj.get_text():continue
 if obj.get_fontsize()<6:small.append(obj.get_text())
 bb=obj.get_window_extent(renderer)
 if bb.width>0 and (bb.x0<-.5 or bb.y0<-.5 or bb.x1>fig.bbox.width+.5 or bb.y1>fig.bbox.height+.5):bad.append(obj.get_text())
assert not small,small
assert not bad,bad
manifest=[]
for fpath in dict.fromkeys(sources):
 stat=fpath.stat()
 manifest.append(dict(path=str(fpath),size=stat.st_size,mtime_ns=stat.st_mtime_ns,sha256=hashlib.sha256(fpath.read_bytes()).hexdigest() if stat.st_size<100_000_000 else 'metadata_only_large_readonly_h5ad'))
(P/'qc/INPUT_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
(P/'qc/BUILD_CHECK.json').write_text(json.dumps(dict(page_mm=[210,297],panels=10,letters=list('ABCDEFGHIJ'),font_min_pt=6,new_inference=False,new_model=False,rows_H=rows,case_scaling=scales,plot_limits=dict(TYMP=[0,2],context=[-1.25,1.25],SPP1abundance=[0,3]),statistics_source='Frozen original BH families retained; no re-correction for chosen panels',versions=dict(python=platform.python_version(),matplotlib=matplotlib.__version__,numpy=np.__version__,pandas=pd.__version__),text_bounds_pass=True),ensure_ascii=False,indent=2))
print('FIGURE_BUILD_PASS',flush=True)
