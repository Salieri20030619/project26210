from pathlib import Path
import json,hashlib,re
import numpy as np,pandas as pd
from scipy.spatial import Delaunay
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection,LineCollection
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle
P=Path(__file__).resolve().parents[1];R=P.parent
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'})
samples=['GLMF1','GLMF2'];case={};manifest=[];display=[];plot_rows=[]
def read(path,**kw):
 path=Path(path);manifest.append(dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),bytes=path.stat().st_size))
 return pd.read_csv(path,sep='\t',**kw)
register=read(R/'phase234_part4_replan_execution_20260921_v1/tables/INPUT_SECTION_REGISTER.tsv')
paths=['Pyrimidine metabolism','Cysteine and methionine metabolism','Arginine and proline metabolism','Glycerophospholipid metabolism']
tfs=['HNF4A','CEBPD','JUNB','STAT3']
for name in samples:
 d=read(R/f'phase234_part4_replan_execution_20260921_v1/inputs/{name}_INPUTS.tsv.gz').set_index('spot_id')
 row=register[register['sample']==name].iloc[0]
 ab=read(row.abundance_path,index_col=0);ab.columns=ab.columns.str.replace('q05cell_abundance_w_sf_','',regex=False);ab=ab.loc[d.index]
 for col in ['Malignant','SPP1-state Macrophage','Endothelial','T cell']:d[col]=ab[col].to_numpy()
 rel='s02_gsva_pilot_glmf1_v1' if name=='GLMF1' else 's02_gsva_sections_v1/GLMF2'
 m=read(R/f'phase37_part4_s02_s04_spatial_metabolism_lr_v1/{rel}/s02_gsva_spot_scores_v1.tsv')
 mat=m[m.source=='CUSTOM_PRIMARY'].pivot(index='spot_id',columns='pathway',values='gsva_ssgsea_score')
 assert set(mat.index)==set(d.index)
 for col in paths:d[col]=mat.loc[d.index,col].to_numpy()
 a=read(P/f'plot_data/{name}_spatial_ATAC_projection.tsv',index_col=0)
 assert set(a.index)==set(d.index)
 for col in tfs+['represented_fraction','projection_eligible']:d[col]=a.loc[d.index,col].to_numpy()
 xy=d[['x_hires_px','y_hires_px']].to_numpy();uv=d[['physical_x_um','physical_y_um']].to_numpy()
 coef=np.linalg.lstsq(np.c_[uv,np.ones(len(uv))],xy,rcond=None)[0]
 assert np.max(np.linalg.norm(np.c_[uv,np.ones(len(uv))]@coef-xy,axis=1))<2
 theta=np.linspace(0,2*np.pi,17)[:-1]
 offset=np.c_[np.cos(theta),np.sin(theta)]*27.5@coef[:2,:]
 vertices=xy[:,None,:]+offset[None,:,:]
 tri=Delaunay(uv).simplices
 edges_len=np.linalg.norm(uv[tri]-uv[tri[:,[1,2,0]]],axis=2)
 tri=tri[(edges_len<151).all(axis=1)]
 edges=np.sort(np.concatenate([tri[:,[0,1]],tri[:,[1,2]],tri[:,[2,0]]]),axis=1)
 uni,cnt=np.unique(edges,axis=0,return_counts=True);boundary=xy[uni[cnt==1]]
 mins=xy.min(axis=0);maxs=xy.max(axis=0);center=(mins+maxs)/2
 width=max(maxs-mins)+150
 case[name]=dict(d=d,xy=xy,vertices=vertices,bound=[center[0]-width/2,center[0]+width/2,center[1]-width/2,center[1]+width/2],boundary=boundary,scale500=np.mean(np.linalg.svd(coef[:2,:],compute_uv=False))*500)
groups=[
 ('A','Cell abundance','cell2location posterior q05',['Malignant','SPP1-state Macrophage','Endothelial','T cell'],['Malignant','SPP1-state macrophages','Endothelial','T cells'],'magma','q05 abundance'),
 ('B','Metabolic pathway scores','Frozen signatures / ssGSEA',paths,['Pyrimidine','Cysteine / methionine','Arginine / proline','Glycerophospholipid'],'magma','ssGSEA score'),
 ('C','ATAC footprint activity','RCC100 reference / chr1 pilot',tfs,tfs,'viridis','HINT activity (PS + TC)')
]
# Native A4 page, three bands. Every map uses the same orientation and crop per case.
fig=plt.figure(figsize=(210/25.4,297/25.4),facecolor='white')
def tx(x,y,t,**kw):return fig.text(x/210,y/297,t,**kw)
def axmm(x,y,w,h):return fig.add_axes([x/210,y/297,w/210,h/297])
tx(9,290,'Supplementary figure | Spatial composition and regulatory activity',fontsize=11,fontweight='bold',va='center')
tx(9,284,'GLMF1 / Patient1 and GLMF2 / Patient7 | external ATAC reference: RCC100',fontsize=8,va='center',color='#444444')
xs=[16,63,110,157]
for g,(letter,title,method,cols,labels,cmap,unit) in enumerate(groups):
 top=275-g*89
 tx(7,top,letter,fontsize=12,fontweight='bold',va='center');tx(16,top,title,fontsize=9.5,fontweight='bold',va='center')
 tx(203,top,method,fontsize=6.6,ha='right',va='center',color='#444444')
 for j,(col,label) in enumerate(zip(cols,labels)):
  combined=np.concatenate([case[n]['d'][col].to_numpy() for n in samples]);combined=combined[np.isfinite(combined)];assert len(combined)>0
  lo,hi=np.quantile(combined,[.02,.98])
  if g==0:lo=0.
  assert hi>lo
  nrm=Normalize(lo,hi,clip=True)
  cmap=plt.get_cmap(cmap).copy();cmap.set_bad('#626262')
  tx(xs[j]+22,top-7,label,fontsize=7.2,ha='center',va='center',fontweight='bold')
  for i,name in enumerate(samples):
   y=top-43-i*34
   o=case[name];v=o['d'][col].to_numpy();ax=axmm(xs[j],y,44,33);ax.set_facecolor('black')
   pc=PolyCollection(o['vertices'],array=np.ma.masked_invalid(v),cmap=cmap,norm=nrm,edgecolors='none',rasterized=True);ax.add_collection(pc)
   ax.add_collection(LineCollection(o['boundary'],colors='#d7d7d7',linewidths=.35,linestyles=(0,(2,2)),alpha=.85))
   x0,x1,y0,y1=o['bound'];ax.set_xlim(x0,x1);ax.set_ylim(y1,y0);ax.set_aspect('equal',adjustable='box');ax.set_xticks([]);ax.set_yticks([])
   for spine in ax.spines.values():spine.set_visible(False)
   # Width-consistent observed tissue axes; label and scale kept inside lower black margin.
   sx=x0+.06*(x1-x0);sy=y1-.055*(y1-y0)
   ax.plot([sx,sx+o['scale500']],[sy,sy],color='white',lw=.8)
   ax.text(sx,sy-.025*(y1-y0),'500 μm',color='white',fontsize=6,va='bottom')
   if j==0:tx(10,y+16.5,name,rotation=90,fontsize=7,fontweight='bold',va='center',ha='center')
   display.append(dict(group=letter,feature=col,sample=name,spots=len(v),vmin=float(lo),vmax=float(hi),clip_low=int((v<lo).sum()),clip_high=int((v>hi).sum()),same_scale_across_cases=True,masked=int(np.isnan(v).sum())))
   plot_rows.append(o['d'][[col]].rename(columns={col:'value'}).assign(sample=name,feature=col).reset_index())
  cbax=axmm(xs[j]+7,top-80,30,1.5)
  cb=fig.colorbar(plt.cm.ScalarMappable(norm=nrm,cmap=cmap),cax=cbax,orientation='horizontal',extend='both' if g else 'max')
  cb.set_ticks([lo,hi]);cb.ax.set_xticklabels([f'{lo:.3g}',f'{hi:.3g}']);cb.ax.tick_params(labelsize=6,length=1.5,pad=1);cb.outline.set_linewidth(.4)
 tx(203,top-85,unit+'; common scale per column',fontsize=6.1,ha='right',color='#444444')
tx(9,6,'Gray, reference coverage <80%. Dashed lines, sampled tissue outline. Scale bars, 500 μm.',fontsize=6.3,color='#333333')
fig.canvas.draw();renderer=fig.canvas.get_renderer()
bad=[]
for obj in fig.findobj(matplotlib.text.Text):
 if not obj.get_visible() or not obj.get_text():continue
 box=obj.get_window_extent(renderer)
 if box.width and (box.x0<-.5 or box.y0<-.5 or box.x1>fig.bbox.width+.5 or box.y1>fig.bbox.height+.5):bad.append(obj.get_text())
assert not bad,bad
base=P/'delivery/26210_spatial_reference_ATAC_chr1_pilot_v1'
fig.savefig(str(base)+'.pdf',dpi=450);fig.savefig(str(base)+'.png',dpi=300);fig.savefig(str(base)+'.svg',dpi=450)
pd.concat(plot_rows).to_csv(P/'plot_data/ALL_DISPLAYED_SPOT_VALUES.tsv.gz',sep='\t',index=False)
(P/'qc/DISPLAY_METADATA.json').write_text(json.dumps(display,indent=2))
(P/'qc/FIGURE_INPUT_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
qc=json.loads((P/'qc/ATAC_PROJECTION_VALIDATION.json').read_text());assert qc['status']=='PASS'
(P/'qc/BUILD_VALIDATION.json').write_text(json.dumps(dict(status='PASS_VISUAL_PENDING',maps=24,patients=2,page_mm=[210,297],spots_per_case={n:len(case[n]['d']) for n in samples},atac_data_in_figure=True,new_inferential_tests=0,font_min_pt=6,source_inputs_unchanged=True),indent=2))
print('FIGURE_READY',base,flush=True)
