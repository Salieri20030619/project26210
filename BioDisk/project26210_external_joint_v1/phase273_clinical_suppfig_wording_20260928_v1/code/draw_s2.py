from pathlib import Path
import pandas as pd,numpy as np,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator,FixedLocator,FixedFormatter
p=Path(__file__).resolve().parents[1]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8.5,'xtick.labelsize':7.5,'ytick.labelsize':7.5,'axes.linewidth':.7,'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
blue='#347A91';red='#B95059';teal='#317F78';gray='#89939D';dark='#202B33'
fig=plt.figure(figsize=(210/25.4,297/25.4),facecolor='white')
fig.text(.055,.963,'Supplementary Fig. S2 | Clinical associations of SPP1 and C4 states',fontsize=12.3,fontweight='bold',color=dark)
fig.text(.055,.94,'Published SPP1 TAM signature and a project-defined C4 transcriptional state',fontsize=8.3,color='#58636B')
fig.add_artist(plt.Line2D([.055,.955],[.921,.921],transform=fig.transFigure,color='#C7D0D5',lw=.8))
def title(letter,x,y,text,sub):
 fig.text(x,y,letter,fontsize=14,fontweight='bold',color=dark)
 fig.text(x+.036,y,text,fontsize=9.4,fontweight='bold',color=dark)
 fig.text(x+.036,y-.021,sub,fontsize=7.7,color='#647078')
def style(ax):
 ax.spines[['top','right']].set_visible(False)
 ax.spines[['left','bottom']].set_color('#53626D')
 ax.tick_params(length=3,width=.6,color='#53626D')
 ax.set_axisbelow(True)
title('A',.055,.89,'SPP1 TAM score and overall survival','CheckMate | nivolumab')
km=pd.read_csv(p/'data/KM_POINTS.tsv',sep='\t');kr=pd.read_csv(p/'data/KM_RISK_TABLE.tsv',sep='\t');ks=pd.read_csv(p/'data/KM_STATISTICS.tsv',sep='\t').iloc[0]
ax=fig.add_axes([.105,.645,.35,.205]);style(ax)
for group,col,count in [('Low',blue,97),('High',red,84)]:
 z=km[km.group==group].sort_values('time')
 ax.step(z.time,z.survival*100,where='post',color=col,lw=1.6,label=f'{group} score (n={count})')
 ax.fill_between(z.time,z.lower*100,z.upper*100,step='post',color=col,alpha=.10,linewidth=0)
 zz=z[z.n_censor>0];ax.plot(zz.time,zz.survival*100,'+',ms=4,mew=.7,color=col)
ax.set(xlim=(0,ks.xmax),ylim=(0,103),ylabel='Overall survival (%)')
ax.set_xticks(kr.time.unique());ax.set_yticks([0,25,50,75,100])
ax.legend(loc='upper right',frameon=False,fontsize=7.4,handlelength=1.8)
ax.text(.025,.11,f'Log-rank P = {ks.logrank_P:.4f}\nFDR = {ks.logrank_FDR:.4f}',transform=ax.transAxes,fontsize=7.6,linespacing=1.5,color=dark)
fig.text(.28,.616,'Time since treatment (months)',ha='center',fontsize=8)
ra=fig.add_axes([.105,.553,.35,.045]);ra.set(xlim=(0,ks.xmax),ylim=(-.5,1.8));ra.axis('off')
ra.text(0,1.65,'Number at risk',fontsize=7.5,color=dark)
for yy,group,col in [(1,'Low',blue),(0,'High',red)]:
 ra.text(-.045,yy,group,transform=ra.get_yaxis_transform(),ha='right',va='center',color=col,fontsize=7.5)
 for row in kr[kr.group==group].itertuples():ra.text(row.time,yy,str(row.n_risk),ha='center',va='center',fontsize=7.3)
title('B',.525,.89,'SPP1 TAM score: treatment interaction','Clinical adjustment | continuous score')
bi=pd.read_csv(p/'data/LITERATURE_INTERACTIONS.tsv',sep='\t')
order=[('CheckMate','OS'),('CheckMate','PFS'),('JAVELIN','PFS')]
bx=fig.add_axes([.705,.655,.15,.174]);style(bx);bx.spines['left'].set_visible(False)
bx.set_xscale('log');bx.set(xlim=(.55,2.5),ylim=(-.45,2.65),yticks=[])
bx.xaxis.set_major_locator(FixedLocator([.6,1,2]));bx.xaxis.set_major_formatter(FixedFormatter(['0.6','1','2']));bx.xaxis.set_minor_locator(NullLocator());bx.axvline(1,color='#79858D',ls='--',lw=.8)
for ii,(ds,ep) in enumerate(order):
 row=bi[(bi.dataset==ds)&(bi.endpoint==ep)].iloc[0];y=2-ii;col=red if ep=='OS' else gray
 bx.errorbar(row.HR_ratio,y,xerr=[[row.HR_ratio-row.CI_low],[row.CI_high-row.HR_ratio]],fmt='o',ms=5,color=col,ecolor=col,capsize=3,lw=1.15)
 ypos=bx.transData.transform((1,y))[1]/fig.bbox.height
 fig.text(.56,ypos,('CM025' if ds=='CheckMate' else 'JAVELIN')+' | '+ep+'\n'+f'n={row.n}',va='center',fontsize=7.7,linespacing=1.5,color=dark)
 fig.text(.88,ypos,f'{row.P:.4f}\n{row.FDR:.4f}',va='center',fontsize=7.3,linespacing=1.5,color=dark)
fig.text(.88,.845,'P / FDR',fontsize=7.5,fontweight='bold',color=dark)
bx.set_xlabel('Interaction HR ratio (95% CI)',labelpad=7,fontsize=8)
fig.text(.56,.6,'CM025: nivolumab vs everolimus\nJAVELIN: avelumab + axitinib vs sunitinib',fontsize=7.2,color='#647078',linespacing=1.6)
fig.add_artist(plt.Line2D([.055,.955],[.529,.529],transform=fig.transFigure,color='#DCE2E5',lw=.7))
title('C',.055,.499,'C4 state: time-varying PFS association','CheckMate | nivolumab | full adjustment')
cv=pd.read_csv(p/'data/C4_TIME_HR_CURVE.tsv',sep='\t');rr=pd.read_csv(p/'data/C4_TIME_RISK_TABLE.tsv',sep='\t');md=pd.read_csv(p/'data/C4_TIME_MODEL.tsv',sep='\t').iloc[0]
cx=fig.add_axes([.105,.248,.35,.211]);style(cx);cx.set_yscale('log')
cx.axvspan(12,24,color='#F1F3F4',zorder=0)
cx.fill_between(cv.month,cv.CI_low,cv.CI_high,color=teal,alpha=.18,linewidth=0)
cx.plot(cv.month,cv.HR,color=teal,lw=1.8)
cx.axhline(1,color='#79858D',ls='--',lw=.8)
cx.set(xlim=(0,24),ylim=(.25,3.2),ylabel='HR per SD of baseline C4 score')
cx.set_xticks([0,6,12,18,24]);cx.yaxis.set_major_locator(FixedLocator([.25,.5,1,2]));cx.yaxis.set_major_formatter(FixedFormatter(['0.25','0.5','1','2']));cx.yaxis.set_minor_locator(NullLocator())
cx.text(.97,.92,'n=164; 145 events\nJoint P=0.000717\nFDR=0.001435',ha='right',va='top',transform=cx.transAxes,fontsize=7.4,color=dark,linespacing=1.55)
fig.text(.28,.217,'Time since treatment (months)',ha='center',fontsize=8)
ra=fig.add_axes([.105,.16,.35,.035]);ra.set(xlim=(0,24),ylim=(-.4,1.25));ra.axis('off')
ra.text(0,1.0,'Number at risk',fontsize=7.5,color=dark)
for row in rr.itertuples():ra.text(row.month,0,str(row.n_risk),ha='center',va='center',fontsize=7.5)
title('D',.525,.499,'C4 state: model robustness','Conditional HR at month 12')
ss=pd.read_csv(p/'data/C4_SENSITIVITY.tsv',sep='\t')
dx=fig.add_axes([.705,.279,.15,.164]);style(dx);dx.spines['left'].set_visible(False)
dx.set_xscale('log');dx.set(xlim=(.35,1.15),ylim=(-.45,2.6),yticks=[])
dx.xaxis.set_major_locator(FixedLocator([.4,.6,1]));dx.xaxis.set_major_formatter(FixedFormatter(['0.4','0.6','1']));dx.xaxis.set_minor_locator(NullLocator());dx.axvline(1,color='#79858D',ls='--',lw=.8)
for ii,label in enumerate(['Full adjustment','Trial stratified','CM025 only']):
 row=ss[ss.model==label].iloc[0];y=2-ii
 dx.errorbar(row.HR,y,xerr=[[row.HR-row.CI_low],[row.CI_high-row.HR]],fmt='o',ms=5,color=teal,ecolor=teal,capsize=3,lw=1.15)
 ypos=dx.transData.transform((1,y))[1]/fig.bbox.height
 textlabel={'Full adjustment':'Full model','Trial stratified':'Trial-stratified','CM025 only':'CM025 only'}[label]
 fig.text(.56,ypos,textlabel+'\n'+f'n={row.n}',va='center',fontsize=7.6,linespacing=1.5,color=dark)
 fig.text(.88,ypos,f'{row.P:.4f}\n{row.FDR:.4f}',va='center',fontsize=7.3,linespacing=1.5,color=dark)
fig.text(.88,.456,'Joint P / FDR',fontsize=7.0,fontweight='bold',ha='left',color=dark)
dx.set_xlabel('12-month HR (95% CI)',labelpad=7,fontsize=8)
fig.text(.56,.22,'Same frozen C4 score; no cutpoint selection\nP values test the full time-dependent effect.',fontsize=7.2,color='#647078',linespacing=1.6)
fig.add_artist(plt.Line2D([.055,.955],[.132,.132],transform=fig.transFigure,color='#C7D0D5',lw=.8))
fig.text(.055,.109,'SPP1 TAM score: SPP1, FN1, CLEC5A, OLR1, IL1RN.',fontsize=7.8,color='#53616B')
fig.text(.055,.09,'C4: a frozen 19-gene state score excluding TYMP and all five SPP1 TAM signature genes.',fontsize=7.8,color='#53616B')
fig.text(.055,.067,'Shading: pointwise 95% CI. In C, the grey background marks later follow-up with fewer patients at risk.',fontsize=7.6,color='#53616B')
fig.text(.055,.049,'Clinical adjustment and test families are detailed in the accompanying legend and statistical tables.',fontsize=7.6,color='#53616B')
pdf=p/'figures/Supplementary_Figure_S2.pdf';fig.savefig(pdf)
fig.savefig(p/'figures/Supplementary_Figure_S2.png',dpi=300)
fig.savefig(p/'qc/S2_preview.png',dpi=150)
plt.close(fig)
(p/'qc/FIGURE_BUILD.json').write_text(json.dumps({'panels':4,'page_mm':[210,297],'minimum_font_pt':7,'pdf':str(pdf),'source_tables':[str(f) for f in (p/'data').glob('*.tsv')],'hypotheses_changed':False},indent=2))
print('FIGURE_DONE',pdf.stat().st_size)
