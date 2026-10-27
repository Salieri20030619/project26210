from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,Normalize
from matplotlib.ticker import FixedLocator,FixedFormatter,NullLocator
from PIL import Image
P=Path(__file__).resolve().parents[1];R=P.parent
read=lambda f:pd.read_csv(P/'data'/f,sep='\t')
blue='#237C9C';red='#B75463';teal='#2E857A';purple='#8A70B0';ink='#22343D';gray='#69757C';light='#DBE3E6'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':6.5,'axes.labelsize':6.5,'xtick.labelsize':6.3,'ytick.labelsize':6.3,'axes.linewidth':.6,'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
fig=plt.figure(figsize=(210/25.4,297/25.4))
def txt(x,y,t,size=6.5,**kw):return fig.text(x/210,1-y/297,t,fontsize=size,color=kw.pop('color',ink),va=kw.pop('va','center'),**kw)
def axmm(x,y,w,h):return fig.add_axes([x/210,1-(y+h)/297,w/210,h/297])
def line(x1,y1,x2,y2,col=light,lw=.6):fig.add_artist(plt.Line2D([x1/210,x2/210],[1-y1/297,1-y2/297],transform=fig.transFigure,color=col,lw=lw))
def panel(letter,x,y,title,sub=None):
 txt(x,y,letter,10,fontweight='bold');txt(x+6,y,title,8.3,fontweight='bold')
 if sub:txt(x+6,y+4,sub,6.3,color=gray)
def style(a):
 a.spines[['top','right']].set_visible(False);a.spines[['left','bottom']].set_color(gray);a.tick_params(length=2,width=.5,pad=1.5,color=gray)
txt(9,7,'TYMP–SPP1 immune programs, spatial niches and clinical outcomes',11.2,fontweight='bold')
line(9,11,201,11)
panel('A',9,17,'Cellular organization of the macrophage–stromal niche')
a=axmm(9,23,119,40);a.imshow(Image.open(P/'assets/niche_illustration.png'));a.axis('off')
txt(14,24,'Tumor nest',6.8,color=red);line(30,25.5,36,34,gray,.5)
txt(73,24,'Fibroblasts / ECM',6.8,color=teal);line(86,26,87,32,gray,.5)
txt(80,62,'CD8 T cells',6.8,color=blue);line(95,60.5,112,52,gray,.5)
txt(16,64,'TYMP/SPP1-associated macrophages',6.8,color=purple);line(54,61.5,46,49,gray,.5)
txt(165,25,'Macrophage RNA programs',7.2,fontweight='bold',ha='center')
txt(140,33,'TYMP',8,fontweight='bold',color=blue);txt(176,33,'SPP1',8,fontweight='bold',color=red)
for yy,l,r in [(39,'IFN-γ  +','IFN-γ  −'),(44,'MHC-I  +','MHC-II  −'),(49,'Phagocytosis  +','')]:
 txt(140,yy,l,6.7);txt(176,yy,r,6.7)
line(140,53,198,53);txt(169,57,'Shared glycolysis association (+)',6.7,color=teal,ha='center')
txt(169,63,'Conceptual schematic; + / − indicate associations',6.3,color=gray,ha='center')
# B: frozen per-patient values, same 16 rows and patient ordering as v1.
panel('B',9,73,'Distinct immune programs across patients','Obradovic | 11 patients | same patient order in each block')
targets=['Interferon Gamma Response','Interferon Alpha Response','GO:0002474__antigen processing and presentation of peptide antigen via MHC class I','GO:0002495__antigen processing and presentation of peptide antigen via MHC class II','GO:0006909__phagocytosis','Complement','Inflammatory Response','TNF-alpha Signaling via NF-kB','IL-6/JAK/STAT3 Signaling','IL-2/STAT5 Signaling','TGF-beta Signaling','Allograft Rejection','GO:0048247__lymphocyte chemotaxis','GO:0050868__negative regulation of T cell activation','GO:0050870__positive regulation of T cell activation','Glycolysis']
labels=['IFN-γ response','IFN-α response','MHC-I antigen presentation','MHC-II antigen presentation','Phagocytosis','Complement','Inflammatory response','TNF-α / NF-κB signaling','IL-6 / JAK / STAT3 signaling','IL-2 / STAT5 signaling','TGF-β signaling','Allograft rejection','Lymphocyte chemotaxis','T-cell activation: negative reg.','T-cell activation: positive reg.','Glycolysis']
patients=['OBRADOVIC:Patient'+str(x) for x in [1,2,3,4,5,6,7,8,'A','B','C']]
eff=read('PANEL_B_PATIENT_VALUES.tsv');sm=read('PANEL_B_STATISTICS.tsv');cm=LinearSegmentedColormap.from_list('corr',['#3C78A1','#A5C6D9','#F8F7F3','#E0A99F','#B34D57']);norm=Normalize(-.4,.4)
for gene,x,col in [('TYMP',66,blue),('SPP1',118,red)]:
 arr=eff[eff.source==gene].pivot(index='target',columns='patient',values='rho').loc[targets,patients].values
 m=axmm(x,88,48,40);m.imshow(arr,cmap=cm,norm=norm,aspect='auto',interpolation='nearest');m.set(xticks=[],yticks=[]);m.set_xticks(np.arange(-.5,11,1),minor=True);m.set_yticks(np.arange(-.5,16,1),minor=True);m.grid(which='minor',color='white',lw=.5);m.tick_params(which='minor',bottom=False,left=False)
 for sp in m.spines.values():sp.set_visible(False)
 txt(x+24,84,gene,7.2,color=col,fontweight='bold',ha='center');line(x,86,x+48,86,col,.9)
 for j,pat in enumerate(patients):txt(x+(j+.5)*48/11,130.6,pat.split('Patient')[1],6.3,ha='center')
 for k,t in enumerate(targets):
  r=sm[(sm.source==gene)&(sm.target==t)].iloc[0];assert np.isclose(np.median(arr[k]),r.median_rho)
  txt(178 if gene=='TYMP' else 193,88+(k+.5)*2.5,'0.00' if abs(r.median_rho)<.005 else f'{r.median_rho:+.2f}',6.3,ha='center',fontweight='bold' if r.FDR<.05 else 'normal',color=ink if r.FDR<.05 else gray)
for i,l in enumerate(labels):txt(63,88+(i+.5)*2.5,l,6.3,ha='right')
txt(186,79,'Median ρ',6.3,ha='center');txt(178,84,'TYMP',6.3,ha='center',color=blue);txt(193,84,'SPP1',6.3,ha='center',color=red)
cax=axmm(68,134,35,1.8);cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cm),cax=cax,orientation='horizontal',ticks=[-.4,0,.4]);cb.outline.set_visible(False);cb.ax.tick_params(labelsize=6.3,length=1.5,pad=.7)
txt(9,135,'Within-patient partial ρ',6.3,color=gray);txt(201,134.5,'Bold median: FDR < 0.05',6.3,ha='right',color=gray)
# C: real spatial coordinates, no smoothing; same full section and scale in all maps.
panel('C',9,145,'Spatial programs in the same tissue','GSE175540 | GSM5924033 | whole-section view')
z=read('SPATIAL_MAP_VALUES.tsv');lims=[z.x_um.min()-100,z.x_um.max()+100,z.y_um.min()-100,z.y_um.max()+100]
quant=[]
for j,(key,label) in enumerate([('SPP1_MAC25_disjoint','SPP1 program'),('FAP_CAF25_disjoint','FAP program'),('ECM_org_heldout','ECM program')]):
 x=9+j*33.5;v=z[key];lo,hi=np.quantile(v,[.02,.98]);quant.append({'program':key,'color_min':lo,'color_max':hi})
 aa=axmm(x,156,32,26);aa.scatter(z.x_um,z.y_um,c=v,s=1.7,cmap='magma',vmin=lo,vmax=hi,edgecolors='none',rasterized=True);aa.set(xlim=lims[:2],ylim=lims[2:]);aa.invert_yaxis();aa.set_aspect('equal');aa.axis('off');txt(x+16,153,label,6.3,ha='center')
 aa.plot([lims[0]+150,lims[0]+1150],[lims[3]-180]*2,color='white',lw=1.4)
 if j==0:aa.text(lims[0]+650,lims[3]-260,'1 mm',fontsize=6.3,ha='center',color='white',bbox=dict(fc='#222222',ec='none',alpha=.6,pad=.6))
ca=axmm(48,184,25,1.6);cl=fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,1),cmap='magma'),cax=ca,orientation='horizontal',ticks=[]);cl.outline.set_visible(False);txt(44,184.8,'Low',6.3,ha='right');txt(76,184.8,'High',6.3)
# D: all qualified patients; fixed original patient order, not sorted by effect.
panel('D',116,145,'CD8 signal at the interface','GSE175540 | 400 μm | 20 patients')
bp=read('BOUNDARY_PATIENT_EFFECTS.tsv');bs=read('BOUNDARY_COHORT_SUMMARY.tsv')
criteria=lambda d:(d.dataset=='GSE175540')&(d.width_um==400)&(d.endpoint=='band_primary_readouts__CD8_RNA')&(d.model=='depth_macro_adjusted')
bd=bp[criteria(bp)].dropna(subset=['effect']).sort_values('patient');st=bs[criteria(bs)].iloc[0];assert len(bd)==20
bd.to_csv(P/'data/PANEL_D_VALUES.tsv',sep='\t',index=False)
d=axmm(127,157,72,24);style(d);d.bar(np.arange(20),bd.effect,color=[blue if v<0 else red for v in bd.effect],width=.76);d.axhline(0,color=gray,lw=.65);d.set(xlim=(-.8,19.8),xticks=[],ylabel='CD8 RNA: inside − outside')
# Patient-direction label placed below the axes.
txt(163,184,'19/20 patients: lower inside',6.3,ha='center',color=gray);txt(163,188,f'P = {st.P:.2g}; FDR = {st.FDR:.2g}',6.3,ha='center')
# E: paired patient estimates show both stromal relationships, including negatives.
panel('E',9,197,'SPP1–stromal program coupling','GSE175540 | 23 patients | depth adjusted')
pa=read('SPATIAL_PATIENT_EFFECTS_ALL.tsv');ss=read('SPATIAL_COHORT_SUMMARY_ALL.tsv');keys=['FAP_CAF25_disjoint__SPP1_MAC25_disjoint','SPP1_MAC25_disjoint__ECM_org_heldout']
q=pa[(pa.dataset=='GSE175540')&(pa.stratum=='overall')&(pa.model=='depth_adjusted')&pa.endpoint.isin(keys)].pivot(index='patient',columns='endpoint',values='effect')[keys];assert q.shape==(23,2)
q.to_csv(P/'data/PANEL_E_PATIENT_PAIRS.tsv',sep='\t')
e=axmm(22,208,74,21);style(e);e.axhline(0,color=gray,ls='--',lw=.5)
for i,row in enumerate(q.values):
 offset=(i-11)*.006;e.plot([0+offset,1+offset],row,color='#ADB8BE',lw=.55,alpha=.65,zorder=1);e.scatter([offset,1+offset],row,c=[purple,teal],s=8,edgecolors='white',lw=.15,zorder=3)
for j,k in enumerate(keys):
 med=q[k].median();e.plot([j-.16,j+.16],[med,med],color=ink,lw=1.5,zorder=4)
e.set(xlim=(-.28,1.28),ylim=(-.16,.8),xticks=[0,1],xticklabels=['FAP ↔ SPP1','SPP1 ↔ ECM'],ylabel='Patient-level partial ρ',yticks=[0,.4,.8]);e.tick_params(axis='x',pad=2)
se=ss[(ss.dataset=='GSE175540')&(ss.stratum=='overall')&(ss.model=='depth_adjusted')&ss.endpoint.isin(keys)].set_index('endpoint')
txt(59,237,'FDR: 4.77 × 10⁻⁷ / 4.88 × 10⁻⁴',6.3,ha='center')
# F: experimental endpoints, individual animal/pool units preserved.
panel('F',116,197,'Published Spp1 intervention','CRC mouse model | WT vs Spp1 knockout')
ko=read('KO_BIOLOGICAL_OR_DECLARED_UNIT_SUMMARY.tsv');rng=np.random.default_rng(329)
for x,pn,gps,title,ylab in [(124,'Ext Fig10H',['wt','Spp1-/-'],'CAF abundance','Count / tissue weight'),(170,'Ext Fig10C',['WT','KO'],'CD8 proliferation','Ki67⁺ CD8 (%)')]:
 f=axmm(x,212,29,17);style(f);txt(x+14.5,207,title,6.3,ha='center')
 zz=ko[ko.panel==pn].copy();assert len(zz)==(13 if pn=='Ext Fig10H' else 10)
 for j,g in enumerate(gps):
  vals=zz[zz.group==g]['mean'].to_numpy();f.scatter(j+np.linspace(-.13,.13,len(vals)),vals,s=12,color=gray if j==0 else teal,edgecolors='white',lw=.25,zorder=3);f.plot([j-.19,j+.19],[vals.mean()]*2,color=ink,lw=1)
 f.set(xlim=(-.45,1.45),xticks=[0,1],xticklabels=['WT','KO'],ylabel=ylab);f.yaxis.label.set_size(6.3)
 if pn=='Ext Fig10H':f.set(ylim=(0,4200),yticks=[0,2000,4000]);f.ticklabel_format(axis='y',style='sci',scilimits=(3,3),useMathText=True)
 else:f.set(ylim=(0,65),yticks=[0,30,60])
 txt(x+14.5,237,'n = 6 / 7' if pn=='Ext Fig10H' else 'n = 5 / 5',6.3,ha='center',color=gray)
# G-H: original clinical curves with all risk-table values shown at selected ticks.
panel('G',9,244,'SPP1 TAM score and survival','CheckMate nivolumab | 181 patients')
panel('H',116,244,'C4 state and PFS over time','CheckMate nivolumab | 164 patients')
km=read('KM_POINTS.tsv');kr=read('KM_RISK_TABLE.tsv');ks=read('KM_STATISTICS.tsv').iloc[0];g=axmm(22,255,78,26);style(g)
for group,col,n in [('Low',blue,97),('High',red,84)]:
 zz=km[km.group==group].sort_values('time');g.step(zz.time,zz.survival*100,where='post',color=col,lw=1,label=f'{group} ({n})');g.fill_between(zz.time,zz.lower*100,zz.upper*100,step='post',color=col,alpha=.1,lw=0);cc=zz[zz.n_censor>0];g.plot(cc.time,cc.survival*100,'+',ms=2,mew=.4,color=col)
g.set(xlim=(0,84),ylim=(0,103),xticks=[0,24,48,72,84],yticks=[0,50,100],ylabel='Overall survival (%)');g.legend(loc='upper right',frameon=False,fontsize=6.3,handlelength=1,labelspacing=.1)
g.text(.025,.07,'P = 0.0148; FDR = 0.0593',transform=g.transAxes,fontsize=6.3)
txt(61,287,'Months',6.3,ha='center')
for y,group,col in [(290,'Low',blue),(293,'High',red)]:
 txt(18,y,group,6.3,ha='right',color=col)
 for rr in kr[(kr.group==group)&kr.time.isin([0,24,48,72,84])].itertuples():txt(22+rr.time/84*78,y,str(rr.n_risk),6.3,ha='center')
cv=read('C4_TIME_HR_CURVE.tsv');rr=read('C4_TIME_RISK_TABLE.tsv');h=axmm(128,255,71,26);style(h);h.set_yscale('log');h.axvspan(12,24,color='#F0F2F3',zorder=0);h.fill_between(cv.month,cv.CI_low,cv.CI_high,color=teal,alpha=.18,lw=0);h.plot(cv.month,cv.HR,color=teal,lw=1.25);h.axhline(1,color=gray,ls='--',lw=.55);h.set(xlim=(0,24),ylim=(.25,3.2),xticks=[0,6,12,18,24],ylabel='HR per SD of C4 score');h.yaxis.set_major_locator(FixedLocator([.25,.5,1,2]));h.yaxis.set_major_formatter(FixedFormatter(['0.25','0.5','1','2']));h.yaxis.set_minor_locator(NullLocator());h.text(.98,.95,'Joint P = 0.000717\nFDR = 0.001435',ha='right',va='top',transform=h.transAxes,fontsize=6.3,linespacing=1.2)
txt(164,287,'Months',6.3,ha='center');txt(120,291.5,'At risk',6.3,ha='right',color=gray)
for row in rr.itertuples():txt(128+row.month/24*71,291.5,str(row.n_risk),6.3,ha='center')
name='26210_immune_spatial_clinical_A4_v2';fig.savefig(P/'delivery'/f'{name}.pdf',dpi=300);fig.savefig(P/'delivery'/f'{name}.png',dpi=300);fig.savefig(P/'figures'/f'{name}.svg');fig.savefig(P/'qc/preview.png',dpi=160)
fig.canvas.draw();render=fig.canvas.get_renderer();outside=[]
for t in fig.findobj(matplotlib.text.Text):
 if t.get_visible() and t.get_text():
  bb=t.get_window_extent(render)
  if bb.x0<0 or bb.y0<0 or bb.x1>fig.bbox.x1 or bb.y1>fig.bbox.y1:outside.append(t.get_text())
(P/'qc/BUILD.json').write_text(json.dumps({'panels':8,'size_mm':[210,297],'spatial_spots':len(z),'spatial_color_limits':quant,'patient_heatmap_cells':352,'outside_text':outside,'new_models':False,'imagegen':'built-in; conceptual schematic only; all measured panels plotted from frozen tables'},indent=2))
print('BUILT',outside)
