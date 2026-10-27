from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd,h5py
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection,PathCollection
from matplotlib.colors import Normalize,ListedColormap,LinearSegmentedColormap
from matplotlib.patches import Rectangle
P=Path(__file__).resolve().parents[1];R=P.parent
A=R/'phase234_part4_replan_execution_20260921_v1';B=R/'phase250_independent_sixpatient_mainline_20260924_v1'
SRC=R/'phase249_independent_spatial_sources_20260924_v1/downloads'
manifest=[]
def record(f):
 f=Path(f);st=f.stat();manifest.append(dict(path=str(f),bytes=st.st_size,mtime_ns=st.st_mtime_ns,sha256=hashlib.sha256(f.read_bytes()).hexdigest() if st.st_size<100_000_000 else None))
def read(f):record(f);return pd.read_csv(f,sep='\t')
reg=read(A/'tables/INPUT_SECTION_REGISTER.tsv');identity=read(B/'tables/CASE_IDENTITY_FROZEN.tsv');qc=read(B/'tables/INPUT_QC.tsv')
selected=qc[qc.tissue=='tumor'].sort_values(['patient','QC_spots','sample'],ascending=[True,False,True]).drop_duplicates('patient').sort_values('patient')
assert len(selected)==6
selection=selected[['patient','sample','QC_spots']].copy();selection['rule']='largest QC spot count per donor; ties lexical sample'
selection.to_csv(P/'plot_data/INDEPENDENT_MAP_SELECTION.tsv',sep='\t',index=False)
plan=dict(panels=list('ABCDEFGHIJ'),page_mm=[210,297],primary='observed tissue maps',scatter_plots=0,statistical_heatmaps=0,new_tests=0,new_models=0,display_regions='Fixed within-section 80th-percentile cutoffs; observed tiles only; no interpolation or smoothing',independent_selection=selection.to_dict('records'),source_figure='Phase254 layout rejected; new tissue-led design; A H&E retained',report_change=False)
(P/'qc/PLAN_BEFORE_DRAWING.json').write_text(json.dumps(plan,ensure_ascii=False,indent=2))
angles=np.deg2rad(np.arange(30,390,60));offset=np.c_[np.cos(angles),np.sin(angles)]*(100/np.sqrt(3))
tissue={}
composition_names=['Malignant','Myeloid','Lymphoid','Vascular','Epithelial','Mast']
composition_colors=['#E9AB42','#527CAF','#A58BBE','#5DAA94','#A8A8A8','#DE9BAD']
groups=[['Malignant'],['C1QC-state Macrophage','SPP1-state Macrophage','Other Macrophage','Monocyte','Dendritic cell'],['T cell','B cell','NK cell'],['Endothelial','SMC/Pericyte'],['Epithelial-Normal'],['Mast cell']]
def package(name,d,uv,xy,img):
 coef=np.linalg.lstsq(np.c_[uv,np.ones(len(uv))],xy,rcond=None)[0]
 err=np.linalg.norm(np.c_[uv,np.ones(len(uv))]@coef-xy,axis=1);assert np.max(err)<2,(name,np.max(err))
 verts=xy[:,None,:]+offset[None,:,:]@coef[:2,:]
 scale500=float(np.mean(np.linalg.svd(coef[:2,:],compute_uv=False))*500)
 extent=np.ptp(xy,axis=0);width=max(extent)+150;c=(xy.min(0)+xy.max(0))/2
 bound=[c[0]-width/2,c[0]+width/2,c[1]-width/2,c[1]+width/2]
 tissue[name]=dict(data=d,uv=uv,xy=xy,vertices=verts,image=img,bound=bound,scale500=scale500,max_registration_error_px=float(err.max()))
for name in ['GLMF1','GLMF2']:
 d=read(A/f'inputs/{name}_INPUTS.tsv.gz');row=reg[reg['sample']==name].iloc[0];record(row.input_path)
 with h5py.File(row.input_path,'r') as f:img=f[f'uns/spatial/{name}/images/hires'][()]
 record(row.abundance_path);ab=pd.read_csv(row.abundance_path,sep='\t',index_col=0)
 ab.columns=ab.columns.str.replace('q05cell_abundance_w_sf_','',regex=False);ab=ab.loc[d.spot_id]
 cl=np.column_stack([ab[g].sum(axis=1).to_numpy() for g in groups]);d['dominant_component']=np.argmax(cl,axis=1)
 assert np.allclose(cl.sum(axis=1),ab.sum(axis=1))
 package(name,d,d[['physical_x_um','physical_y_um']].to_numpy(),d[['x_hires_px','y_hires_px']].to_numpy(),img)
 for target,out in [('spp1_state_no_spp1','context_state_region'),('endothelial_q05','spp1_endothelial_region')]:
  x='polyamine_mta_context' if out=='context_state_region' else 'V1_SPP1_q05';qx=float(d[x].quantile(.8));qy=float(d[target].quantile(.8))
  d[out]=(d[x].ge(qx).astype(int)+2*d[target].ge(qy).astype(int)).to_numpy()
  tissue[name][out+'_thresholds']={x:qx,target:qy}
 d[['spot_id','x_hires_px','y_hires_px','dominant_component','TYMP_RNA','polyamine_mta_context','V1_SPP1_q05','spp1_state_no_spp1','context_state_region','immunosuppression_core','invasion_emt','spp1_endothelial_region']].to_csv(P/f'plot_data/{name}_tissue_layers.tsv',sep='\t',index=False)
names=selected['sample'].tolist()+['pHG_1.1','pHG_3.1']
for name in names:
 d=read(B/f'inputs/{name}_MAINLINE_INPUTS.tsv.gz')
 posfile=SRC/f'{name}_tissue_positions_list.csv';record(posfile)
 pos=pd.read_csv(posfile,header=None,names=['barcode','in_tissue','array_row','array_col','px_row','px_col'])
 if pos.barcode.iloc[0]=='barcode':pos=pos.iloc[1:]
 pos=pos.set_index('barcode').loc[d.barcode].apply(pd.to_numeric)
 sfpath=SRC/f'{name}_scalefactors_json.json';record(sfpath);sf=json.loads(sfpath.read_text())
 ipath=SRC/f'{name}_tissue_hires_image.png';record(ipath);img=np.array(Image.open(ipath).convert('RGB'))
 xy=pos[['px_col','px_row']].to_numpy()*sf['tissue_hires_scalef'];uv=d[['x_um','y_um']].to_numpy()
 package(name,d,uv,xy,img)
 qx=float(d.polyamine_mta_context.quantile(.8));qy=float(d.spp1_state_no_spp1.quantile(.8))
 d['context_state_region']=(d.polyamine_mta_context.ge(qx).astype(int)+2*d.spp1_state_no_spp1.ge(qy).astype(int)).to_numpy()
 tissue[name]['context_state_region_thresholds']={'polyamine_mta_context':qx,'spp1_state_no_spp1':qy}
 d[['barcode','x_um','y_um','TYMP_RNA','SPP1_RNA','polyamine_mta_context','spp1_state_no_spp1','context_state_region']].to_csv(P/f'plot_data/{name}_tissue_layers.tsv',sep='\t',index=False)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'axes.labelsize':7,'xtick.labelsize':6.2,'ytick.labelsize':6.2,'axes.linewidth':.45,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
fig=plt.figure(figsize=(210/25.4,297/25.4));allmaps=[]
def axmm(x,y,w,h):return fig.add_axes([x/210,y/297,w/210,h/297])
def tx(x,y,s,**kw):return fig.text(x/210,y/297,s,**kw)
def head(letter,title,x,y,size=7.5):
 tx(x,y,letter,fontsize=10,fontweight='bold',va='center');tx(x+4.5,y,title,fontsize=size,fontweight='bold',va='center')
def maps(name,x,y,w,h=None,values=None,cmap=None,norm=None,scale=True):
 h=h or w;obj=tissue[name];d=obj['data'];ax=axmm(x,y,w,h)
 ax.imshow(obj['image'],alpha=1 if values is None else .12,interpolation='bilinear')
 if values is not None:
  v=d[values].to_numpy() if isinstance(values,str) else np.asarray(values);assert len(v)==len(d) and np.isfinite(v).all()
  pc=PolyCollection(obj['vertices'],array=v,cmap=cmap,norm=norm,edgecolors='none',antialiased=False,rasterized=True,alpha=.96);ax.add_collection(pc)
 x0,x1,y0,y1=obj['bound'];ax.set_xlim(x0,x1);ax.set_ylim(y1,y0);ax.set_aspect('equal');ax.axis('off')
 if scale:
  sx=x0+(x1-x0)*.055;sy=y1-(y1-y0)*.055
  ax.plot([sx,sx+obj['scale500']],[sy,sy],lw=1.4,color='#151515')
  ax.text(sx,sy-(y1-y0)*.033,'500 μm',fontsize=6,va='bottom',color='#151515',bbox=dict(facecolor='white',edgecolor='none',alpha=.88,pad=.3))
 allmaps.append(dict(sample=name,feature='H&E' if values is None else str(values),axis=ax));return ax
def cbar(x,y,w,cmap,norm,ticks,label):
 ax=axmm(x,y,w,1.8);cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=ax,orientation='horizontal',extend='both' if norm.vmin<0 else 'max')
 cb.outline.set_linewidth(.4);cb.set_ticks(ticks);cb.ax.tick_params(length=1.5,pad=1,labelsize=6);tx(x+w/2,y-5,label,ha='center',fontsize=6.2)
def swatch(x,y,color,label,size=6.1):
 fig.add_artist(Rectangle((x/210,(y-.6)/297),1.8/210,1.8/297,transform=fig.transFigure,facecolor=color,edgecolor='none'));tx(x+2.6,y,label,va='center',fontsize=size)
tympcm=LinearSegmentedColormap.from_list('tymp',['#F1EAF4','#CEACDD','#8750A2','#421A61'])
spp1cm=LinearSegmentedColormap.from_list('mac',['#EDF5EE','#A1CDAE','#499A78','#07513D'])
regioncolors=['#E5E3DF','#D59A42','#4C9F9A','#9A406E'];regioncm=ListedColormap(regioncolors)
vascularcolors=['#E5E3DF','#327EAA','#D87562','#764D95'];vascularcm=ListedColormap(vascularcolors)
xs=[8,47,86,125,164];titles=['Histology','Tissue composition','TYMP expression','Metabolic context','SPP1 macrophages']
for i,title in enumerate(titles):head('ABCDE'[i],title,xs[i],287,size=7.2)
for row,name in enumerate(['GLMF1','GLMF2']):
 y=242-row*39;tx(8,y+36.5,name,fontsize=7,fontweight='bold')
 maps(name,xs[0],y,35,values=None)
 maps(name,xs[1],y,35,values='dominant_component',cmap=ListedColormap(composition_colors),norm=Normalize(-.5,5.5),scale=False)
 maps(name,xs[2],y,35,values='TYMP_RNA',cmap=tympcm,norm=Normalize(0,2),scale=False)
 maps(name,xs[3],y,35,values='polyamine_mta_context',cmap='RdBu_r',norm=Normalize(-1.25,1.25),scale=False)
 maps(name,xs[4],y,35,values='V1_SPP1_q05',cmap=spp1cm,norm=Normalize(0,3),scale=False)
cbar(90,196,26,tympcm,Normalize(0,2),[0,1,2],'log1p(CP10k)')
cbar(129,196,26,'RdBu_r',Normalize(-1.25,1.25),[-1,0,1],'Mean gene z-score')
cbar(168,196,26,spp1cm,Normalize(0,3),[0,1.5,3],'Estimated abundance')
tx(8,196,'B  Dominant mapped component',fontsize=6.4)
for i,(name,color) in enumerate(zip(composition_names,composition_colors)):swatch(8+(i%3)*23,191-(i//3)*4.5,color,name)
head('F','Context and SPP1-state regions',8,175,size=7.2)
head('G','Spatial functional programs',76,175,size=7.2)
head('H','SPP1 and endothelium',145,175,size=7.2)
for row,name in enumerate(['GLMF1','GLMF2']):
 y=135-row*39
 tx(12,y+34,name,fontsize=6.7,fontweight='bold');maps(name,12,y,32,values='context_state_region',cmap=regioncm,norm=Normalize(-.5,3.5),scale=True)
 tx(76,y+33.8,name,fontsize=6.4,fontweight='bold')
 maps(name,76,y,28,values='immunosuppression_core',cmap='RdBu_r',norm=Normalize(-1.25,1.25),scale=False)
 maps(name,108,y,28,values='invasion_emt',cmap='RdBu_r',norm=Normalize(-1.25,1.25),scale=False)
 tx(146,y+34,name,fontsize=6.7,fontweight='bold');maps(name,146,y,32,values='spp1_endothelial_region',cmap=vascularcm,norm=Normalize(-.5,3.5),scale=True)
tx(90,165.5,'Immune regulation',ha='center',fontsize=6.1);tx(122,165.5,'Invasion / EMT',ha='center',fontsize=6.1)
for i,(label,col) in enumerate(zip(['Other','Context only','State only','Both high'],regioncolors)):swatch(8+(i%2)*27,89-(i//2)*4.5,col,label)
tx(8,79,'High: within-section top 20%',fontsize=6)
cbar(94,90,24,'RdBu_r',Normalize(-1.25,1.25),[-1,0,1],'Mean gene z-score')
for i,(label,col) in enumerate(zip(['Other','SPP1 only','Endothelium only','Both high'],vascularcolors)):swatch(143+(i%2)*29,89-(i//2)*4.5,col,label)
tx(143,79,'High: within-section top 20%',fontsize=6)
head('I','Independent tissue maps',8,73,size=7.5);tx(8,68.8,'E-MTAB-12767; region colors as in F',fontsize=6.2)
for k,row in enumerate(selected.itertuples()):
 x=8+(k%3)*37;y=36.7-(k//3)*32.8;name=row.sample;patient=row.patient.split(':')[-1]
 tx(x,y+29,patient+'  '+name,fontsize=6.2,fontweight='bold');maps(name,x,y,28,values='context_state_region',cmap=regioncm,norm=Normalize(-.5,3.5),scale=True)
head('J','Paired SPP1 expression',127,73,size=7.5);tx(133,68.8,'Adjacent',fontsize=6.5);tx(169,68.8,'Tumor',fontsize=6.5)
tumors=dict(zip(selected.patient,selected['sample']));pairrows=[('HG1','pHG_1.1',tumors['EMTAB12767:HG1']),('HG3','pHG_3.1',tumors['EMTAB12767:HG3'])]
for k,(patient,adj,tumor) in enumerate(pairrows):
 y=36.7-k*32.8;tx(126,y+14,patient,fontsize=6.5,fontweight='bold',rotation=90,va='center')
 maps(adj,132,y,28,values='SPP1_RNA',cmap=tympcm,norm=Normalize(0,5),scale=True);maps(tumor,165,y,28,values='SPP1_RNA',cmap=tympcm,norm=Normalize(0,5),scale=True)
cx=axmm(197,16,1.8,36);cb=fig.colorbar(plt.cm.ScalarMappable(cmap=tympcm,norm=Normalize(0,5)),cax=cx,orientation='vertical',extend='max',ticks=[0,2.5,5])
cb.outline.set_linewidth(.4);cb.ax.tick_params(labelsize=6,length=1.5,pad=1);cb.set_label('log1p(CP10k)',fontsize=6,labelpad=2)
assert len(allmaps)==28,len(allmaps)
assert not any(isinstance(x,PathCollection) for ax in fig.axes for x in ax.collections)
fig.canvas.draw();renderer=fig.canvas.get_renderer();bad=[];small=[]
for obj in fig.findobj(matplotlib.text.Text):
 if not obj.get_visible() or not obj.get_text():continue
 if obj.get_fontsize()<6:small.append(obj.get_text())
 bb=obj.get_window_extent(renderer)
 if bb.width>0 and (bb.x0<-.5 or bb.y0<-.5 or bb.x1>fig.bbox.width+.5 or bb.y1>fig.bbox.height+.5):bad.append(obj.get_text())
assert not small,small
assert not bad,bad
base=P/'delivery/26210_Spatial_A4_tissue_10panels_v2'
fig.savefig(str(base)+'.pdf',dpi=450);fig.savefig(str(base)+'.svg',dpi=450);fig.savefig(str(base)+'.png',dpi=300);fig.savefig(str(base)+'.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'})
(P/'qc/INPUT_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
(P/'qc/DISPLAY_METADATA.json').write_text(json.dumps({name:{k:v for k,v in obj.items() if k in ['scale500','max_registration_error_px','bound'] or k.endswith('_thresholds')} for name,obj in tissue.items()},ensure_ascii=False,indent=2))
(P/'qc/BUILD_CHECK.json').write_text(json.dumps(dict(status='PASS',page_mm=[210,297],panels=10,tissue_maps=len(allmaps),observed_hex_tile_maps=26,original_HE_maps=2,scatter_plots=0,statistical_heatmaps=0,model_fitting=0,statistical_tests=0,minimum_font=6,text_bounds_pass=True,independent_patients=6,paired_patients=2,independent_sample_selection=selection.to_dict('records'),composition_counts={name:dict(zip(composition_names,np.bincount(tissue[name]['data'].dominant_component,minlength=6).tolist())) for name in ['GLMF1','GLMF2']}),ensure_ascii=False,indent=2))
print('BUILD_PASS spatial_maps=28 panels=10 scatter=0 statistical_heatmap=0')
