from pathlib import Path
import json,subprocess
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
w=Path(__file__).resolve().parents[1];q=w/'qc';out=w/'figures'
assert (q/'R_INDEPENDENT_VALIDATION.txt').read_text().startswith('PASS')
read=lambda p:pd.read_csv(p,sep='\t')
a=read(w/'data/TYMP_ACTIVITY_PATIENT_EFFECTS.tsv');s=read(q/'TYMP_ACTIVITY_SUMMARY.tsv');b=read(w/'data/TERMINAL_ENZYME_PATIENT_EFFECTS.tsv');bs=read(q/'TERMINAL_ENZYME_SUMMARY.tsv')
blue='#397E9C';red='#C2576D';gray='#BBC3CB'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'savefig.facecolor':'white'})
def fmt(x):return f'{x:.3g}'
def save(fig,name,title,foot):
 fig.suptitle(title,x=.04,y=.974,ha='left',fontsize=16,fontweight='bold')
 fig.text(.04,.03,foot,fontsize=8.5,ha='left',linespacing=1.5)
 for ext in ['png','pdf']:fig.savefig(out/(name+'.'+ext),dpi=360)
 plt.close(fig);subprocess.run(['pdftoppm','-r','125','-singlefile','-png',str(out/(name+'.pdf')),str(q/(name+'_pdf'))],check=True)
 print('SAVED',name,flush=True)
# A: all activity targets plus the separately labeled RNA endpoint, before and after technical adjustment.
targets=['CASP9','CHMP1A','PRSS8','PXN','RRM2B_RNA'];fig,axs=plt.subplots(1,2,figsize=(13.5,7.8));fig.subplots_adjust(left=.12,right=.985,top=.81,bottom=.19,wspace=.40)
bound=np.ceil(abs(a.effect).max()*10)/10+.1
for ax,model,title in zip(axs,['unadjusted','technical_adjusted'],['A  Within-patient association','B  Technical-covariate sensitivity']):
 for i,target in enumerate(targets):
  d=a[(a.target==target)&(a.model==model)].sort_values('patient');r=s[(s.target==target)&(s.model==model)].iloc[0]
  jitter=np.linspace(-.12,.12,len(d));ax.scatter(d.effect,i+jitter,s=20,c=blue,alpha=.65,zorder=2)
  ax.plot([r.median_CI_low,r.median_CI_high],[i,i],c=red,lw=2,zorder=3);ax.scatter([r.median_effect],[i],marker='D',s=48,c=red,edgecolors='white',lw=.6,zorder=4)
  stats=f'rho={r.median_effect:+.3f}; positive {int(r.positive)}/{int(r.n)}; P={fmt(r.P)}'
  stats+=f'; q={fmt(r.FDR)}' if np.isfinite(r.FDR) else ' (single RNA endpoint)'
  ax.text(-bound*.98,i+.34,stats,fontsize=7.7,va='center')
 ax.axvline(0,c='#8F9CA6',lw=.75,ls='--');ax.axhline(3.64,c='#AAB3BC',lw=.7)
 ax.set_xlim(-bound,bound);ax.set_ylim(4.65,-.5);ax.set_yticks(range(5),['CASP9 activity','CHMP1A activity','PRSS8 activity','PXN activity','RRM2B RNA'],fontsize=9)
 ax.set_title(title,loc='left',fontsize=11,pad=15);ax.set_xlabel('Association with TYMP RNA (patient rho)')
save(fig,'01_section2_5_TYMP_activity_support_v1','TYMP expression links to several inferred protein activities','20 independent patients | dots = patient effects; diamond = median; line = patient-bootstrap 95% interval for median\nSensitivity adjusts ranked library size and detected genes. Protein BH family = 4 per model; RRM2B is RNA, not protein activity.')
# B: full 13-gene x 3-method terminal screen, with all outcomes visible.
genes=['TYMP','RRM2B','TK1','SLC29A1','SLC29A2','ODC1','SRM','SMS','MTAP','CHKA','PCYT1A','PCYT2','CEPT1'];methods=['Monocle3','Slingshot','scTour']
mat=bs.pivot(index='gene',columns='method',values='median_effect').loc[genes,methods].values;lim=np.ceil(abs(mat).max()*10)/10
fig,ax=plt.subplots(figsize=(11.3,12));fig.subplots_adjust(left=.15,right=.86,top=.89,bottom=.12)
im=ax.imshow(mat,cmap='RdBu_r',vmin=-lim,vmax=lim,aspect='auto')
for i,gene in enumerate(genes):
 for j,method in enumerate(methods):
  r=bs[(bs.gene==gene)&(bs.method==method)].iloc[0]
  ax.text(j,i,f'{r.median_effect:+.3f} | {int(r.positive)}/{int(r.n)} positive\nP={fmt(r.P)}; q={fmt(r.FDR)}',ha='center',va='center',fontsize=8.2,color='white' if abs(r.median_effect)>.5 else '#263442')
  if r.P<.05:ax.add_patch(Rectangle((j-.48,i-.48),.96,.96,fill=False,ec='#CA8B2E',lw=2.3))
ax.set_xticks(range(3),methods,fontsize=11);ax.set_yticks(range(13),genes,fontsize=10);ax.tick_params(length=0)
cb=fig.colorbar(im,ax=ax,fraction=.04,pad=.025);cb.set_label('Median terminal-minus-early log2(CPM + 1)',fontsize=9)
save(fig,'02_section2_7_enzyme_terminal_complete_v1','Enzyme RNA changes vary across pseudotime methods','Each cell summarizes 20 independent patients | terminal >=80th percentile versus early <=20th percentile within patient\nPseudobulk RNA; orange border = nominal P<0.05. All 39 gene-method endpoints retained; BH across the full 39-test family.')
# C: expose every nominal endpoint, including negative directions, with patient pairs.
nom=bs[bs.P<.05].sort_values(['method','gene']);assert len(nom)==4
fig,axs=plt.subplots(2,2,figsize=(11,9));fig.subplots_adjust(left=.09,right=.97,top=.85,bottom=.15,hspace=.58,wspace=.3)
for ax,r in zip(axs.flat,nom.itertuples()):
 d=b[(b.method==r.method)&(b.gene==r.gene)].sort_values('patient')
 for rr in d.itertuples():ax.plot([0,1],[rr.early_log2_CPM1,rr.late_log2_CPM1],c=gray,lw=.75,zorder=1)
 ax.scatter(np.zeros(len(d)),d.early_log2_CPM1,s=24,c=blue,zorder=2);ax.scatter(np.ones(len(d)),d.late_log2_CPM1,s=24,c=red,zorder=2)
 ax.scatter([0,1],[d.early_log2_CPM1.median(),d.late_log2_CPM1.median()],marker='_',s=350,c='#202C36',linewidths=2.2,zorder=3)
 ax.set_xlim(-.3,1.3);ax.set_xticks([0,1],['Early','Terminal']);ax.set_ylabel('Pseudobulk log2(CPM + 1)')
 ax.set_title(r.method+' | '+r.gene,loc='left',fontweight='bold',fontsize=11,pad=27)
 ax.text(0,1.01,f'delta={r.median_effect:+.3f}; positive {r.positive}/{r.n}\nP={fmt(r.P)}; q={fmt(r.FDR)}',transform=ax.transAxes,fontsize=8,va='bottom')
save(fig,'03_section2_7_nominal_patient_pairs_v1','Nominal terminal changes show patient-level heterogeneity','One line per independent patient | all four nominal endpoints from the 39-test screen are shown, including negative directions\nThese are exploratory, method-specific RNA findings; full FDR and the complete screen accompany this figure.')
