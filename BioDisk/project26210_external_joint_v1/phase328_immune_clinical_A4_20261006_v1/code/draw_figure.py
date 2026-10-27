from pathlib import Path
import json, hashlib, shutil
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse,Circle,PathPatch,FancyBboxPatch
from matplotlib.path import Path as MPath
from matplotlib.colors import LinearSegmentedColormap,Normalize
from matplotlib.ticker import FixedLocator,FixedFormatter,NullLocator
P=Path(__file__).resolve().parents[1]; ROOT=P.parent
I=ROOT/'phase314_TYMP_participation_immune_exploration_20261004_v1/tables'
J=ROOT/'phase321_TYMP_SPP1_conditional_gene_robustness_20261005_v1/tables'
C=ROOT/'phase273_clinical_suppfig_wording_20260928_v1/data'
sources=[I/'ALL_PATIENT_EFFECTS.tsv',I/'ALL_EFFECT_SUMMARIES.tsv',J/'FORMAL_SUMMARIES.tsv']+[C/f for f in ['KM_POINTS.tsv','KM_RISK_TABLE.tsv','KM_STATISTICS.tsv','LITERATURE_INTERACTIONS.tsv','C4_TIME_HR_CURVE.tsv','C4_TIME_RISK_TABLE.tsv','C4_TIME_MODEL.tsv']]
identities=[]
for f in sources:
 identities.append({'path':str(f),'sha256':hashlib.sha256(f.read_bytes()).hexdigest()});shutil.copy2(f,P/'data'/f.name)
read=lambda name:pd.read_csv(P/'data'/name,sep='\t')
eff=read('ALL_PATIENT_EFFECTS.tsv');summ=read('ALL_EFFECT_SUMMARIES.tsv');formal=read('FORMAL_SUMMARIES.tsv')
# Prespecified 15-program immune view, with previously reported glycolysis appended.
targets=['Interferon Gamma Response','Interferon Alpha Response','GO:0002474__antigen processing and presentation of peptide antigen via MHC class I','GO:0002495__antigen processing and presentation of peptide antigen via MHC class II','GO:0006909__phagocytosis','Complement','Inflammatory Response','TNF-alpha Signaling via NF-kB','IL-6/JAK/STAT3 Signaling','IL-2/STAT5 Signaling','TGF-beta Signaling','Allograft Rejection','GO:0048247__lymphocyte chemotaxis','GO:0050868__negative regulation of T cell activation','GO:0050870__positive regulation of T cell activation','Glycolysis']
labels=['IFN-γ response','IFN-α response','MHC-I antigen presentation','MHC-II antigen presentation','Phagocytosis','Complement','Inflammatory response','TNF-α / NF-κB signaling','IL-6 / JAK / STAT3 signaling','IL-2 / STAT5 signaling','TGF-β signaling','Allograft rejection','Lymphocyte chemotaxis','T-cell activation: negative reg.','T-cell activation: positive reg.','Glycolysis']
patients=['OBRADOVIC:Patient'+str(x) for x in [1,2,3,4,5,6,7,8,'A','B','C']]
sel=eff[(eff.branch=='immune')&(eff.condition=='L2_broad')&(eff.model=='core')&eff.target.isin(targets)].copy()
assert len(sel)==16*22 and sel.rho.notna().all() and not sel.duplicated(['source','target','patient']).any()
sel.to_csv(P/'data/PANEL_B_PATIENT_VALUES.tsv',sep='\t',index=False)
sm=summ[(summ.branch=='immune')&(summ.condition=='L2_broad')&(summ.model=='core')&summ.target.isin(targets)].copy()
sm.to_csv(P/'data/PANEL_B_STATISTICS.tsv',sep='\t',index=False)
blue='#237C9C';red='#B75463';teal='#2E857A';ink='#22343D';gray='#69757C';light='#DBE3E6'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7.2,'axes.labelsize':7.2,'xtick.labelsize':6.8,'ytick.labelsize':6.8,'axes.linewidth':.65,'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
fig=plt.figure(figsize=(210/25.4,297/25.4))
def txt(x,y,t,size=7.2,**kw): return fig.text(x/210,1-y/297,t,fontsize=size,color=kw.pop('color',ink),va=kw.pop('va','center'),**kw)
def axmm(x,y,w,h):return fig.add_axes([x/210,1-(y+h)/297,w/210,h/297])
def line(x1,y1,x2,y2,col=light,lw=.7):fig.add_artist(plt.Line2D([x1/210,x2/210],[1-y1/297,1-y2/297],transform=fig.transFigure,color=col,lw=lw))
def panel(letter,x,y,title,sub=None):
 txt(x,y,letter,11.3,fontweight='bold');txt(x+7,y,title,9,fontweight='bold')
 if sub:txt(x+7,y+4.7,sub,6.8,color=gray)
def style(a):
 a.spines[['top','right']].set_visible(False);a.spines[['bottom','left']].set_color(gray);a.tick_params(length=2.5,width=.6,color=gray,pad=2)
txt(10,10,'TYMP–SPP1 macrophage programs and clinical associations',12,fontweight='bold')
txt(10,16,'Single-cell immune programs  /  SPP1-related signature  /  C4 state',7.4,color=gray)
line(10,20,200,20)
panel('A',10,26,'Shared metabolism, distinct immune programs','Obradovic RCC | 11 untreated patients | schematic of observed associations')
# Native vector biological illustration; connectors are associations, never causal arrows.
a=axmm(10,35,190,42);a.set(xlim=(0,190),ylim=(0,42));a.axis('off')
theta=np.linspace(0,2*np.pi,220);rad=1+.045*np.sin(5*theta)+.025*np.cos(9*theta)
xx=95+25*rad*np.cos(theta);yy=23+16*rad*np.sin(theta)
a.fill(xx,yy,color='#F6EAE1',ec='#C9A994',lw=1.2,zorder=1)
a.add_patch(Ellipse((96,27),17,11,angle=-20,fc='#D5C9DE',ec='#AA95B7',lw=.8))
a.add_patch(Ellipse((98,27),5,4,fc='#B7A4C6',ec='none'))
for x,y,r in [(79,28,1.8),(85,35,1.1),(113,29,1.6),(105,13,1.3),(80,17,1.5),(111,21,1.2)]:
 a.add_patch(Circle((x,y),r,fc='#E8B79A',ec='#BF886D',lw=.6))
 for dx,dy in [(-.45,.1),(.4,.55),(.15,-.45)]:a.add_patch(Circle((x+dx,y+dy),.2,fc='#AC765F',ec='none'))
# Surface antigen-presentation motifs and extracellular peptides.
for x,y,s in [(74,32,-1),(120,29,1)]:
 a.plot([x,x+s*3],[y,y+2],color='#829693',lw=1)
 a.plot([x+s*3,x+s*5],[y+2,y+1],color='#829693',lw=1)
 a.plot([x+s*3,x+s*3],[y+2,y+4],color='#829693',lw=1)
 a.plot([x+s*3.4,x+s*5],[y+3,y+3.6],color='#D09653',lw=1.7)
for x,y,col,label in [(86,21,blue,'TYMP'),(105,21,red,'SPP1')]:
 a.add_patch(FancyBboxPatch((x-8,y-3),16,6,boxstyle='round,pad=.3,rounding_size=1.3',fc='white',ec=col,lw=1))
 a.text(x,y,label,ha='center',va='center',fontsize=8.3,color=col,fontweight='bold',zorder=5)
# Left / right observed programs.
left=[(35,'IFN-γ response'),(25,'MHC-I antigen presentation'),(15,'Phagocytosis')]
right=[(35,'IFN-γ response'),(20,'MHC-II antigen presentation')]
for y,label in left:
 a.plot([78,66,60],[21,y,y],color=blue,lw=.85,zorder=0)
 a.text(2,y,label,fontsize=8,va='center',color=ink)
 a.text(57,y,'+',ha='center',va='center',fontsize=11,color=blue,fontweight='bold')
for y,label in right:
 a.plot([113,124,126],[21,y,y],color=red,lw=.85,zorder=0)
 a.text(134,y,label,fontsize=8,va='center',color=ink)
 a.text(129,y,'−',ha='center',va='center',fontsize=11,color=red,fontweight='bold',bbox=dict(fc='white',ec='none',pad=.3))
a.plot([86,86,94],[18,5,5],color=teal,lw=.8);a.plot([105,105,97],[18,5,5],color=teal,lw=.8)
a.text(95,3.5,'+  Glycolysis',fontsize=8,ha='center',va='center',color=teal,bbox=dict(fc='white',ec='none',pad=1))
a.text(2,3,'+ / −  Positive / negative association',fontsize=6.6,color=gray,va='center')
a.text(189,3,'RNA programs',fontsize=6.6,color=gray,ha='right',va='center')
panel('B',10,83,'Patient-level immune landscape','16 programs | same 11 patients in both gene blocks | broad macrophage definition')
cmap=LinearSegmentedColormap.from_list('corr',['#3C78A1','#A5C6D9','#F8F7F3','#E0A99F','#B34D57']);norm=Normalize(-.4,.4)
for gene,x,col in [('TYMP',70,blue),('SPP1',121,red)]:
 arr=sel[sel.source==gene].pivot(index='target',columns='patient',values='rho').loc[targets,patients].values
 m=axmm(x,100,46,56);m.imshow(arr,cmap=cmap,norm=norm,aspect='auto',interpolation='nearest');m.set(xticks=[],yticks=[])
 m.set_xticks(np.arange(-.5,11,1),minor=True);m.set_yticks(np.arange(-.5,16,1),minor=True);m.grid(which='minor',color='white',lw=.5);m.tick_params(which='minor',bottom=False,left=False)
 for sp in m.spines.values():sp.set_visible(False)
 txt(x+23,95,gene,8.2,fontweight='bold',color=col,ha='center');line(x,98,x+46,98,col,1.2)
 for j,pat in enumerate(patients):txt(x+(j+.5)*46/11,159,pat.split('Patient')[1],6.7,ha='center')
 for k,t in enumerate(targets):
  row=sm[(sm.source==gene)&(sm.target==t)].iloc[0]
  assert np.isclose(np.median(arr[k]),row.median_rho)
  xt=178 if gene=='TYMP' else 193
  txt(xt,100+(k+.5)*3.5,('0.00' if abs(row.median_rho)<.005 else f'{row.median_rho:+.2f}'),6.6,ha='center',fontweight='bold' if row.FDR<.05 else 'normal',color=ink if row.FDR<.05 else gray)
for i,label in enumerate(labels):txt(66,100+(i+.5)*3.5,label,7.1,ha='right')
txt(185.5,91,'Median ρ',6.8,ha='center');txt(178,96,'TYMP',6.5,ha='center',color=blue);txt(193,96,'SPP1',6.5,ha='center',color=red)
line(68,152.5,199,152.5,'#849498',.75)
txt(70,164,'Patient',6.8,color=gray)
cax=axmm(93,166,49,2.3);cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax,orientation='horizontal',ticks=[-.4,0,.4]);cb.outline.set_visible(False);cb.ax.tick_params(labelsize=6.5,length=2,pad=1)
txt(10,168,'Partial rank correlation (ρ)',6.8,color=gray)
txt(199,166,'Bold median: FDR < 0.05',6.6,ha='right',color=gray)
line(10,174,200,174)
panel('C',10,180,'SPP1 TAM score and survival','CheckMate nivolumab | n = 181; 123 events')
panel('D',110,180,'C4 state and PFS over time','CheckMate nivolumab | n = 164; 145 events')
km=read('KM_POINTS.tsv');kr=read('KM_RISK_TABLE.tsv');ks=read('KM_STATISTICS.tsv').iloc[0]
c=axmm(24,190,75,39);style(c)
for group,col,count in [('Low',blue,97),('High',red,84)]:
 z=km[km.group==group].sort_values('time');c.step(z.time,z.survival*100,where='post',color=col,lw=1.1,label=f'{group} ({count})');c.fill_between(z.time,z.lower*100,z.upper*100,step='post',color=col,alpha=.11,lw=0)
 zz=z[z.n_censor>0];c.plot(zz.time,zz.survival*100,'+',ms=2.6,mew=.5,color=col)
c.set(xlim=(0,ks.xmax),ylim=(0,103),ylabel='Overall survival (%)',xticks=[0,24,48,72,84],yticks=[0,25,50,75,100]);c.legend(loc='upper right',frameon=False,fontsize=6.4,handlelength=1.1,labelspacing=.25)
c.text(.035,.09,f'Log-rank P = {ks.logrank_P:.4f}\nFDR = {ks.logrank_FDR:.4f}',fontsize=6.6,transform=c.transAxes,linespacing=1.4)
txt(61.5,238,'Time since treatment (months)',7,ha='center')
txt(24,243,'Number at risk',6.5,color=gray)
for yy,g,col in [(247,'Low',blue),(251,'High',red)]:
 txt(20,yy,g,6.5,color=col,ha='right')
 for row in kr[(kr.group==g)&kr.time.isin([0,24,48,72,84])].itertuples():txt(24+row.time/ks.xmax*75,yy,str(row.n_risk),6.5,ha='center')
cv=read('C4_TIME_HR_CURVE.tsv');rr=read('C4_TIME_RISK_TABLE.tsv');md=read('C4_TIME_MODEL.tsv').iloc[0]
d=axmm(125,190,74,39);style(d);d.set_yscale('log');d.axvspan(12,24,color='#F0F2F3',zorder=0);d.fill_between(cv.month,cv.CI_low,cv.CI_high,color=teal,alpha=.18,lw=0);d.plot(cv.month,cv.HR,color=teal,lw=1.4);d.axhline(1,color=gray,ls='--',lw=.65)
d.set(xlim=(0,24),ylim=(.25,3.2),ylabel='HR per SD of C4 score',xticks=[0,6,12,18,24]);d.yaxis.set_major_locator(FixedLocator([.25,.5,1,2]));d.yaxis.set_major_formatter(FixedFormatter(['0.25','0.5','1','2']));d.yaxis.set_minor_locator(NullLocator())
d.text(.97,.95,f'Joint P = {md.joint_P:.6f}\nFDR = {md.joint_FDR:.6f}',ha='right',va='top',transform=d.transAxes,fontsize=6.6,linespacing=1.4)
txt(162,238,'Time since treatment (months)',7,ha='center');txt(125,243,'Number at risk',6.5,color=gray)
for row in rr.itertuples():txt(125+row.month/24*74,247,str(row.n_risk),6.5,ha='center')
txt(125,251,'Shading: 95% CI; grey: reduced late risk set',6.3,color=gray)
line(10,255,200,255)
panel('E',10,261,'SPP1 TAM score: treatment interaction across clinical endpoints')
bi=read('LITERATURE_INTERACTIONS.tsv');e=axmm(68,268,57,16);style(e);e.spines['left'].set_visible(False);e.set_xscale('log');e.set(xlim=(.6,2.4),ylim=(-.5,2.5),yticks=[]);e.xaxis.set_major_locator(FixedLocator([.6,1,2]));e.xaxis.set_major_formatter(FixedFormatter(['0.6','1','2']));e.xaxis.set_minor_locator(NullLocator());e.axvline(1,color=gray,ls='--',lw=.7)
for x,t in [(133,'HR ratio (95% CI)'),(176,'P'),(193,'FDR')]:txt(x,267,t,6.4,ha='left' if x==133 else 'center',fontweight='bold')
for i,(ds,ep) in enumerate([('CheckMate','OS'),('CheckMate','PFS'),('JAVELIN','PFS')]):
 row=bi[(bi.dataset==ds)&(bi.endpoint==ep)].iloc[0];y=2-i;col=red if ep=='OS' else '#8A969D';e.errorbar(row.HR_ratio,y,xerr=[[row.HR_ratio-row.CI_low],[row.CI_high-row.HR_ratio]],fmt='o',ms=4,color=col,capsize=2,lw=1)
 yy=268+(2.5-y)/3*16
 txt(17,yy,('CM025' if ds=='CheckMate' else 'JAVELIN')+f'  {ep}  (n = {row.n})',7)
 txt(133,yy,f'{row.HR_ratio:.2f} ({row.CI_low:.2f}–{row.CI_high:.2f})',6.6)
 txt(176,yy,f'{row.P:.4f}',6.5,ha='center');txt(193,yy,f'{row.FDR:.4f}',6.5,ha='center')
txt(96.5,292,'Interaction HR ratio',6.8,ha='center')
name='26210_immune_clinical_A4_v1'
fig.savefig(P/'delivery'/f'{name}.pdf');fig.savefig(P/'delivery'/f'{name}.png',dpi=300);fig.savefig(P/'figures'/f'{name}.svg');fig.savefig(P/'qc/preview.png',dpi=150)
fig.canvas.draw();renderer=fig.canvas.get_renderer();outside=[]
for obj in fig.findobj(matplotlib.text.Text):
 if obj.get_visible() and obj.get_text():
  b=obj.get_window_extent(renderer);fb=fig.bbox
  if b.x0 < -1 or b.y0 < -1 or b.x1>fb.x1+1 or b.y1>fb.y1+1:outside.append(obj.get_text())
assert not outside,outside
(P/'qc/BUILD.json').write_text(json.dumps({'size_mm':[210,297],'panels':5,'minimum_font_pt':6.3,'programs':16,'patients':11,'heatmap_cells':352,'new_statistical_tests':False,'outside_canvas':outside,'inputs':identities},indent=2))
print('FIGURE_READY',str(P/'delivery'/f'{name}.pdf'))
