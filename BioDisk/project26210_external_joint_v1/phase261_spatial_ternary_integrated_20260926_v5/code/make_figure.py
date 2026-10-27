from pathlib import Path
import json,hashlib,subprocess
import numpy as np,pandas as pd,h5py
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection,PathCollection,LineCollection
from matplotlib.colors import Normalize,ListedColormap,LinearSegmentedColormap
from matplotlib.patches import Rectangle
P=Path(__file__).resolve().parents[1];R=P.parent
A=R/'phase234_part4_replan_execution_20260921_v1';B=R/'phase250_independent_sixpatient_mainline_20260924_v1'
SRC=R/'phase249_independent_spatial_sources_20260924_v1/downloads'
manifest=[]
def gate(label):
 r=subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True)
 with (P/'logs/memory_gate.log').open('a') as f:f.write(label+'\t'+r.stdout+r.stderr)
 if r.returncode:raise RuntimeError(r.stdout+r.stderr)

def record(f):
 f=Path(f);st=f.stat();manifest.append(dict(path=str(f),bytes=st.st_size,mtime_ns=st.st_mtime_ns,sha256=hashlib.sha256(f.read_bytes()).hexdigest() if st.st_size<100_000_000 else None))
def read(f):record(f);return pd.read_csv(f,sep='\t')
reg=read(A/'tables/INPUT_SECTION_REGISTER.tsv');identity=read(B/'tables/CASE_IDENTITY_FROZEN.tsv');qc=read(B/'tables/INPUT_QC.tsv')
selected=qc[qc.tissue=='tumor'].sort_values(['patient','QC_spots','sample'],ascending=[True,False,True]).drop_duplicates('patient').sort_values('patient')
assert len(selected)==6
selection=selected[['patient','sample','QC_spots']].copy();selection['rule']='largest QC spot count per donor; ties lexical sample'
selection.to_csv(P/'plot_data/INDEPENDENT_MAP_SELECTION.tsv',sep='\t',index=False)
plan=dict(panels=list('ABCDEFGHI'),page_mm=[210,297],primary='observed tissue maps plus functional GSEA',scatter_plots=2,statistical_heatmaps=0,new_tests=0,new_models=0,display_regions='Fixed within-section 80th-percentile cutoffs; observed tiles only; no interpolation or smoothing',independent_selection=selection.to_dict('records'),source_figure='Phase259 tissue and GSEA data; Phase260 two separate ternary cases; no overlay',report_change=False)
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
 gate('load_'+name)
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
 gate('load_'+name)
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

# New layout: omit separate raw H&E and group composition/TYMP.
gate('gsea_display_data')
GSEA=R/'phase256_spatial_region_function_20260926_v1'
allresults=read(GSEA/'tables/GSEA_FULL_RESULTS.tsv.gz')
gsea_curves={};gsea_rows=[]
cohorts=[('GSE175540',23),('TJRCC',12)]
pathways=[('Lysosome','Lysosome','#247F9E'),('Cholesterol metabolism','Chol. met.','#AB597A')]
for cohort,npat in cohorts:
 rankpath=GSEA/f'gsea/{cohort}__absolute__exclude_definition_genes.rnk'
 setpath=GSEA/f'gsea/{cohort}__absolute__exclude_definition_genes__KEGG.gmt'
 ranked=read(rankpath);record(setpath)
 genes=ranked.gene.to_numpy();metric=ranked['rank'].to_numpy()
 sets={f[0]:f[2:] for f in [l.split('\t') for l in setpath.read_text().splitlines()]}
 assert len(set(genes))==len(genes) and np.all(np.diff(metric)<=1e-12)
 for term,short,color in pathways:
  selected_gene=np.isin(genes,sets[term]);N=len(metric);Nh=int(selected_gene.sum())
  weights=np.abs(metric)*selected_gene;weights/=weights.sum()
  walk=np.cumsum(weights-(~selected_gene)/(N-Nh))
  hit=np.flatnonzero(selected_gene);peak_index=int(np.argmax(np.abs(walk)));ES=float(walk[peak_index])
  row=allresults[(allresults.dataset==cohort)&(allresults.model=='absolute')&(allresults.version=='exclude_definition_genes')&(allresults.library=='KEGG')&(allresults.term==term)]
  assert len(row)==1;row=row.iloc[0]
  assert row.patient_n==npat and abs(ES-row.ES)<1e-8 and row.NES>0 and row.FDR_GSEA<.05
  assert abs(walk[-1])<1e-10 and Nh==row['size']
  xx=(np.arange(N)+1)/N
  gsea_curves[(cohort,term)]=dict(x=xx,walk=walk,hits=(hit+1)/N,ES=ES,NES=float(row.NES),FDR=float(row.FDR_GSEA),P=float(row.P_GSEA),genes=N,set_size=Nh,color=color)
  gsea_rows.append(dict(cohort=cohort,patient_n=npat,pathway=term,model='absolute_SPP1_q05_top20_vs_rest',version='exclude_definition_genes',ES=ES,NES=float(row.NES),P_GSEA=float(row.P_GSEA),FDR_GSEA=float(row.FDR_GSEA),ranked_genes=N,set_size=Nh,rank_source=str(rankpath),set_source=str(setpath)))
  pd.DataFrame({'gene':genes,'rank':metric,'rank_fraction':xx,'gene_set_hit':selected_gene,'running_ES':walk}).to_csv(P/f'plot_data/GSEA_{cohort}_{short.replace(" ","_")}.tsv.gz',sep='\t',index=False)
pd.DataFrame(gsea_rows).to_csv(P/'plot_data/GSEA_FOUR_CURVES_STATISTICS.tsv',sep='\t',index=False)
for suffix in ['pdf','svg','png','tiff']:
 record(R/f'phase255_spatial_tissue_figure_20260925_v2/delivery/26210_Spatial_A4_tissue_10panels_v2.{suffix}')
 record(R/f'phase258_spatial_GSEA_figure_20260926_v3/delivery/26210_Spatial_A4_GSEA_v3.{suffix}')


head('A','Tissue context',6,290)
head('B','Metabolism',72,290,size=7)
head('C','SPP1-state',105,290,size=7)
head('D','Context / state',138,290,size=6.8)
head('E','SPP1 / endoth.',171,290,size=6.8)
tx(19.5,284,'Composition',ha='center',fontsize=6.5)
tx(52.5,284,'TYMP RNA',ha='center',fontsize=6.5)
for i,name in enumerate(['GLMF1','GLMF2']):
 y=249-i*34;tx(6,y+29,name,fontsize=6.6,fontweight='bold')
 for x,feat,cm,no in [(6,'dominant_component',ListedColormap(composition_colors),Normalize(-.5,5.5)),(39,'TYMP_RNA',tympcm,Normalize(0,2)),(72,'polyamine_mta_context','RdBu_r',Normalize(-1.25,1.25)),(105,'V1_SPP1_q05',spp1cm,Normalize(0,3)),(138,'context_state_region',regioncm,Normalize(-.5,3.5)),(171,'spp1_endothelial_region',vascularcm,Normalize(-.5,3.5))]:
  maps(name,x,y,27,values=feat,cmap=cm,norm=no,scale=x==6)
for i,(name,color) in enumerate(zip(composition_names,composition_colors)):
 swatch(6+(i%2)*23,208-(i//2)*4.5,color,name,size=6)
cbar(43,208,20,tympcm,Normalize(0,2),[0,1,2],'log1p(CP10k)')
cbar(76,208,20,'RdBu_r',Normalize(-1.25,1.25),[-1,0,1],'Mean gene z-score')
cbar(109,208,20,spp1cm,Normalize(0,3),[0,1.5,3],'q05 abundance')
for i,(label,col) in enumerate(zip(['Other','Context','State','Both'],regioncolors)):
 swatch(137+(i%2)*17,208-(i//2)*5,col,label,size=6)
for i,(label,col) in enumerate(zip(['Other','SPP1','Endoth.','Both'],vascularcolors)):
 swatch(171+(i%2)*17,208-(i//2)*5,col,label,size=6)
head('F','Regional GSEA',6,193)
curve_axes=[];stem_collections=[]
ylim_high=np.ceil(max(obj['walk'].max() for obj in gsea_curves.values())/.1)*.1
ylim_low=min(-.025,np.floor(min(obj['walk'].min() for obj in gsea_curves.values())/.1)*.1)
gsea_cmap=LinearSegmentedColormap.from_list('ES_blue_gray_purple',['#426B98','#C8CBD0','#AA5EAB'])
gsea_norm=Normalize(0,ylim_high);stem_rows=[]
for r,(cohort,npat) in enumerate(cohorts):
 tx(56+r*98,186.5,('TJ-RCC' if cohort=='TJRCC' else cohort)+f' (n = {npat})',ha='center',fontsize=6.8,fontweight='bold')
 for c,(term,short,color) in enumerate(pathways):
  obj=gsea_curves[(cohort,term)];x=12+(r*2+c)*49
  tx(x+19.5,181.5,short,ha='center',fontsize=6.6,fontweight='bold')
  ax=axmm(x,151,39,23)
  curve_x=np.r_[0,obj['x']];curve_y=np.r_[0,obj['walk']]
  points=np.column_stack([curve_x,curve_y]);segments=np.stack([points[:-1],points[1:]],axis=1)
  lc=LineCollection(segments,cmap=gsea_cmap,norm=gsea_norm,linewidths=1.25,zorder=3)
  lc.set_array((curve_y[:-1]+curve_y[1:])/2);ax.add_collection(lc)
  hit_idx=np.rint(obj['hits']*obj['genes']).astype(int)-1;heights=obj['walk'][hit_idx]
  stems=ax.vlines(obj['hits'],0,heights,colors=gsea_cmap(gsea_norm(heights)),lw=.5,alpha=.92,zorder=2)
  ax.axhline(0,color='#333333',lw=.7,ls=(0,(4,3)),zorder=4)
  ax.set_xlim(0,1);ax.set_ylim(ylim_low,ylim_high*1.04)
  ax.set_yticks([0,ylim_high/2,ylim_high]);ax.set_yticklabels([f'{v:g}' for v in [0,ylim_high/2,ylim_high]])
  ax.tick_params(axis='y',labelsize=6,length=2,pad=1.4)
  ax.set_xticks([0,.5,1]);ax.set_xticklabels(['0','0.5','1'])
  ax.tick_params(axis='x',labelsize=6,length=1.5,pad=1)
  for side in ['left','right','top','bottom']:
   ax.spines[side].set_visible(True);ax.spines[side].set_linewidth(.5);ax.spines[side].set_color('#414141')
  if c==0:ax.set_ylabel('ES',fontsize=6.3,labelpad=1.5)
  ax.set_title(f"NES {obj['NES']:.2f}   q {obj['FDR']:.1e}",loc='left',fontsize=6,pad=2)
  for k,hv in zip(hit_idx,heights):
   stem_rows.append(dict(cohort=cohort,pathway=term,rank_index_zero_based=int(k),rank_fraction=float(obj['x'][k]),stem_bottom=0.,stem_top=float(hv),color_scalar_ES=float(hv)))
  curve_axes.append(ax);stem_collections.append(stems)
pd.DataFrame(stem_rows).to_csv(P/'plot_data/GSEA_COLORED_STEMS.tsv',sep='\t',index=False)
tx(105,141.5,'Ranked genes (fraction)',ha='center',fontsize=6.2)
cbax=axmm(173,141.5,24,1.3)
es_cb=fig.colorbar(plt.cm.ScalarMappable(norm=gsea_norm,cmap=gsea_cmap),cax=cbax,orientation='horizontal',ticks=[0,ylim_high/2,ylim_high])
es_cb.outline.set_linewidth(.4);es_cb.ax.tick_params(labelsize=6,length=1,pad=1)
tx(168,141.5,'ES',fontsize=6,va='center')
head('G','Relative spatial distances from malignant-rich spots',6,134)
ternary_axes=[];ternary_counts={}
T=R/'phase260_spatial_ternary_preview_20260926_v1'
h=np.sqrt(3)/2
def xy(a,b,c):return [.5*a+c,h*a]
for i,(sample,color) in enumerate([('GLMF1','#4B94B5'),('GLMF2','#D47C6C')]):
 src=T/f'plot_data/{sample}_TERNARY_POINTS.tsv';dd=read(src)
 (P/'plot_data'/src.name).write_bytes(src.read_bytes())
 ternary_counts[sample]=len(dd)
 tx(52+i*105,127.5,f'{sample} | {len(dd):,} spots',ha='center',fontsize=6.8,fontweight='bold',color=color)
 ax=axmm(6+i*105,70,92,54)
 for t in [.2,.4,.6,.8]:
  for p,q in [(xy(t,1-t,0),xy(t,0,1-t)),(xy(1-t,t,0),xy(0,t,1-t)),(xy(1-t,0,t),xy(0,1-t,t))]:
   ax.plot([p[0],q[0]],[p[1],q[1]],color='#DBDEE2',lw=.45,zorder=0)
 ax.plot([.5,0,1,.5],[h,0,0,h],color='#44484E',lw=.85,zorder=2)
 ax.text(.91,.48,'SPP1-state',ha='center',va='center',rotation=-60,fontsize=6.5,fontweight='bold')
 ax.text(.09,.48,'Endoth.',ha='center',va='center',rotation=60,fontsize=6.5,fontweight='bold')
 ax.text(.5,-.10,'T cells',ha='center',va='top',fontsize=6.5,fontweight='bold')
 for val in [0,20,40,60,80,100]:
  t=val/100
  ax.text(t,-.017,str(val),ha='center',va='top',fontsize=6)
  ax.text(.5*t-.025,h*t+.008,str(100-val),ha='right',va='center',fontsize=6)
  ax.text(.5+.5*t+.025,h*(1-t)+.008,str(100-val),ha='left',va='center',fontsize=6)
 ax.scatter(dd.ternary_x,dd.ternary_y,s=3.5,c=color,alpha=.35,lw=0,rasterized=True,zorder=1)
 ax.set_xlim(-.16,1.16);ax.set_ylim(-.16,1.03);ax.set_aspect('equal');ax.axis('off')
 ternary_axes.append(ax)
tx(105,67,'Distance shares (%); larger share = relatively farther',ha='center',fontsize=6.2)
head('H','Independent tissue maps',6,63)
tx(6,59.5,'E-MTAB-12767; region colors as in D',fontsize=6)
for k,row in enumerate(selected.itertuples()):
 x=6+(k%3)*34;y=32-(k//3)*27;name=row.sample;patient=row.patient.split(':')[-1]
 tx(x,y+23.5,patient+'  '+name,fontsize=6,fontweight='bold')
 maps(name,x,y,22,values='context_state_region',cmap=regioncm,norm=Normalize(-.5,3.5),scale=True)
head('I','Paired SPP1 expression',122,63)
tx(130,59.5,'Adjacent',fontsize=6.3);tx(163,59.5,'Tumor',fontsize=6.3)
tumors=dict(zip(selected.patient,selected['sample']))
pairrows=[('HG1','pHG_1.1',tumors['EMTAB12767:HG1']),('HG3','pHG_3.1',tumors['EMTAB12767:HG3'])]
for k,(patient,adj,tumor) in enumerate(pairrows):
 y=32-k*27;tx(120,y+11,patient,fontsize=6.3,fontweight='bold',rotation=90,va='center')
 maps(adj,126,y,22,values='SPP1_RNA',cmap=tympcm,norm=Normalize(0,5),scale=True)
 maps(tumor,158,y,22,values='SPP1_RNA',cmap=tympcm,norm=Normalize(0,5),scale=True)
cx=axmm(188,12,1.8,34)
cb=fig.colorbar(plt.cm.ScalarMappable(cmap=tympcm,norm=Normalize(0,5)),cax=cx,orientation='vertical',extend='max',ticks=[0,2.5,5])
cb.outline.set_linewidth(.4);cb.ax.tick_params(labelsize=6,length=1.5,pad=1);cb.set_label('log1p(CP10k)',fontsize=6,labelpad=2)
assert len(allmaps)==22 and len(curve_axes)==4 and len(ternary_axes)==2
assert ternary_counts=={'GLMF1':1282,'GLMF2':3113}
assert sum(isinstance(c,PathCollection) for ax in fig.axes for c in ax.collections)==2
fig.canvas.draw();renderer=fig.canvas.get_renderer();bad=[];small=[]
for obj in fig.findobj(matplotlib.text.Text):
 if not obj.get_visible() or not obj.get_text():continue
 if obj.get_fontsize()<6:small.append(obj.get_text())
 bb=obj.get_window_extent(renderer)
 if bb.width>0 and (bb.x0<-.5 or bb.y0<-.5 or bb.x1>fig.bbox.width+.5 or bb.y1>fig.bbox.height+.5):bad.append(obj.get_text())
assert not small,small
assert not bad,bad
gate('figure_export')
base=P/'delivery/26210_Spatial_A4_GSEA_Ternary_v5'
for ext,dpi in [('pdf',450),('svg',450),('png',300),('tiff',600)]:
 fig.savefig(str(base)+'.'+ext,dpi=dpi,**({'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {}))
(P/'qc/INPUT_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
(P/'qc/DISPLAY_METADATA.json').write_text(json.dumps({name:{k:v for k,v in obj.items() if k in ['scale500','max_registration_error_px','bound'] or k.endswith('_thresholds')} for name,obj in tissue.items()},ensure_ascii=False,indent=2))
(P/'qc/BUILD_CHECK.json').write_text(json.dumps(dict(status='PASS',page_mm=[210,297],panels=9,tissue_maps=22,gsea_curves=4,ternary_plots=2,ternary_spots=ternary_counts,overlay=False,colored_gene_stems=True,bottom_rug_strips=False,minimum_font=6,text_bounds_pass=True,new_tests=0,new_models=0),indent=2))
print('BUILD_PASS: A4, 9 panels, 22 tissue maps, 4 GSEA curves, 2 individual ternary plots',flush=True)
