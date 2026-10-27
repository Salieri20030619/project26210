from pathlib import Path
import json,hashlib,subprocess,sys
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,ListedColormap
from matplotlib.patches import Rectangle
from scipy.cluster.hierarchy import linkage,dendrogram,fcluster
w=Path(__file__).resolve().parents[1];old=w.parent/'phase127_section2_4_pyscenic_redraw_remote_v1';q=w/'qc'
methods=['GENIE3','GRNBoost2'];factors=['F'+str(i) for i in range(1,8)];candidates=['ILF2(+)','HOXD3(+)','HNF4A(+)']
paths=[old/'data/PATIENT_FACTOR_AUC_MEDIANS.tsv',old/'qc/GENIE3_PATIENT_BALANCED_AUC.tsv',old/'qc/GRNBoost2_PATIENT_BALANCED_AUC.tsv',old/'qc/HEATMAP_SELECTION.tsv',old/'scripts/draw_section24.py',old/'figures/03_B_cNMF_patient_balanced_regulon_heatmap_v1.pdf']
origmanifest={r['path']:r for r in json.loads((old/'qc/SOURCE_MANIFEST.json').read_text())}
for p in paths:
 if str(p) in origmanifest:assert hashlib.sha256(p.read_bytes()).hexdigest()==origmanifest[str(p)]['sha256']
sel=pd.read_csv(paths[3],sep='\t');selected=list(dict.fromkeys(sel.regulon.tolist()))
for c in candidates:
 if c not in selected:selected.append(c)
pf=pd.read_csv(paths[0],sep='\t');Z={};errors={}
for method in methods:
 gm=pf[pf.method==method].groupby('factor')[selected].median().reindex(factors).T
 orig=pd.read_csv(old/'qc'/(method+'_PATIENT_BALANCED_AUC.tsv'),sep='\t',index_col=0).loc[selected,factors]
 errors[method]=float(np.max(abs(gm.values-orig.values)));assert errors[method]<1e-12
 sd=gm.std(axis=1,ddof=0);Z[method]=gm.sub(gm.mean(axis=1),axis=0).div(sd.where(sd>0,1),axis=0)
Zall=np.column_stack([np.column_stack([Z[method][fac].values for method in methods]) for fac in factors])
tree=linkage(Zall,method='average',metric='euclidean',optimal_ordering=True)
leaves=dendrogram(tree,no_plot=True)['leaves'];regs=[selected[i] for i in leaves];n=len(regs)
expected=['TCF4(+)','FOXN2(+)','TFAP4(+)','ETV1(+)','NR1H4(+)','HNF4A(+)','ILF2(+)','WT1(+)','NR2F6(+)','DBP(+)','JUND(+)','HOXD3(+)','FOXC2(+)','BACH1(+)','ZBTB7A(+)','STAT3(+)','ZFX(+)','PPARG(+)','FOSL1(+)','ATF1(+)','IRF7(+)','ESRRG(+)','HOXB6(+)','MAFK(+)']
assert regs==expected,'Row ordering differs from original figure'
cluster_raw=fcluster(tree,t=4,criterion='maxclust')[leaves]
ids=list(dict.fromkeys(cluster_raw));clusters=np.array([ids.index(c)+1 for c in cluster_raw]);sizes={int(c):int(np.sum(clusters==c)) for c in np.unique(clusters)}
for c in sizes:
 ix=np.flatnonzero(clusters==c);assert np.all(np.diff(ix)==1)
nprof=pf[pf.method=='GENIE3'].groupby('factor').patient.nunique().reindex(factors)
assert nprof.tolist()==[9,12,2,4,4,2,4] and n==24
assert pf.patient.nunique()==21 and len(pf)==74
fc=dict(zip(factors,['#BA7868','#7295BE','#79A697','#C9AD63','#A58BB4','#8DB6C3','#C78FAD']))
mc={'GENIE3':'#AB7764','GRNBoost2':'#5E91AB'};cc={1:'#596D77',2:'#BA865B',3:'#4F91A6',4:'#9583B0'};cellcolor='#7A9D8C'
div=LinearSegmentedColormap.from_list('regulon_div',['#5E936F','#F7F6F2','#927093'])
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.linewidth':.65,'pdf.fonttype':42,'savefig.facecolor':'white'})
fig=plt.figure(figsize=(14,10));left,bottom,width,height=.20,.155,.46,.64
ax=fig.add_axes([left,bottom,width,height]);im=ax.imshow(Zall[leaves],cmap=div,vmin=-2,vmax=2,aspect='auto',interpolation='nearest')
ax.set_yticks(np.arange(n),regs,fontsize=9);ax.yaxis.tick_right();ax.tick_params(axis='y',length=0,pad=6)
ax.set_xticks(np.arange(14),['G3','GB2']*7,fontsize=8);ax.tick_params(axis='x',length=0,pad=6)
for s in ax.spines.values():s.set_visible(False)
for tick,reg in zip(ax.get_yticklabels(),regs):
 if reg in candidates:tick.set_fontweight('bold')
for k in np.arange(1.5,13,2):ax.axvline(k,color='white',lw=2)
for k in np.flatnonzero(np.diff(clusters)!=0)+.5:ax.axhline(k,color='white',lw=2)
# Original row dendrogram, supplemented with explicit cluster labels and matching strip.
ta=fig.add_axes([.09,bottom,.077,height]);dendrogram(tree,orientation='left',no_labels=True,color_threshold=0,above_threshold_color='#888888',ax=ta);ta.invert_yaxis();ta.axis('off')
strip=fig.add_axes([.178,bottom,.014,height]);strip.imshow((clusters-1)[:,None],cmap=ListedColormap(list(cc.values())),vmin=0,vmax=3,aspect='auto',interpolation='nearest');strip.set_xticks([]);strip.set_yticks([])
for s in strip.spines.values():s.set_visible(False)
for c,count in sizes.items():
 ix=np.flatnonzero(clusters==c);y=bottom+height*(1-(ix.mean()+.5)/n)
 fig.text(.025,y,f'C{c} ({count})',ha='left',va='center',fontsize=10,fontweight='bold',color=cc[c])
 y1=bottom+height*(1-(ix[0])/n);y0=bottom+height*(1-(ix[-1]+1)/n)
 fig.lines.append(plt.Line2D([.076,.083,.083,.076],[y0,y0,y1,y1],transform=fig.transFigure,color=cc[c],lw=1))
fig.text(.025,.809,'Cluster',ha='left',va='bottom',fontsize=9,fontweight='bold')
# Explicit annotation rows; cell type is homogeneous in this malignant-only analysis.
def annotation(y,h,values,palette,label):
 aa=fig.add_axes([left,y,width,h]);aa.imshow(np.array(values)[None,:],cmap=ListedColormap(palette),vmin=0,vmax=max(1,len(palette)-1),aspect='auto',interpolation='nearest');aa.set_xticks([]);aa.set_yticks([])
 for k in np.arange(1.5,13,2):aa.axvline(k,color='white',lw=1)
 for s in aa.spines.values():s.set_visible(False)
 fig.text(left-.008,y+h/2,label,ha='right',va='center',fontsize=9)
 return aa
cellax=annotation(.87,.023,[0]*14,[cellcolor],'Cell type')
cellax.text(6.5,0,'Malignant epithelial',ha='center',va='center',fontsize=10,color='white',fontweight='bold')
progax=annotation(.838,.024,np.repeat(np.arange(7),2),list(fc.values()),'cNMF program')
for j,fac in enumerate(factors):progax.text(j*2+.5,0,fac,ha='center',va='center',fontsize=10)
algax=annotation(.808,.020,[0,1]*7,list(mc.values()),'Algorithm')
for j in range(14):algax.text(j,0,['G3','GB2'][j%2],ha='center',va='center',fontsize=7,color='white')
for j,fac in enumerate(factors):
 fig.text(left+width*(j*2+1)/14,.111,'n='+str(int(nprof[fac])),ha='center',fontsize=9)
fig.text(left-.008,.111,'Patients',ha='right',fontsize=9)
# Right-side complete legends, each maps to an actual annotation in the heatmap.
lx=.80
ca=fig.add_axes([lx,.862,.16,.018]);fig.colorbar(im,cax=ca,orientation='horizontal',ticks=[-2,0,2],extend='both')
ca.set_title('Within-method row z-score',fontsize=9,pad=8);ca.tick_params(labelsize=8,length=3)
fig.text(lx,.815,'Cell type',fontweight='bold',fontsize=10)
def key(y,color,label):
 fig.add_artist(Rectangle((lx,y-.006),.012,.015,facecolor=color,edgecolor='#495660',lw=.5,transform=fig.transFigure))
 fig.text(lx+.02,y+.001,label,va='center',fontsize=9)
key(.79,cellcolor,'Malignant epithelial')
fig.text(lx,.747,'cNMF program',fontweight='bold',fontsize=10)
for i,fac in enumerate(factors):key(.719-i*.025,fc[fac],fac)
fig.text(lx,.52,'Algorithm',fontweight='bold',fontsize=10)
for i,method in enumerate(methods):key(.491-i*.026,mc[method],['G3: GENIE3','GB2: GRNBoost2'][i])
fig.text(lx,.412,'Regulon cluster',fontweight='bold',fontsize=10)
for i,(c,count) in enumerate(sizes.items()):key(.383-i*.026,cc[c],f'C{c} ({count} regulons)')
fig.text(lx,.245,'Bold gene labels',fontweight='bold',fontsize=10)
fig.text(lx,.222,'ILF2, HOXD3, HNF4A\nHistorical candidates',fontsize=9,va='top',linespacing=1.5)
fig.text(.025,.961,'B  Patient-balanced regulon activity across cNMF programs',fontsize=16,fontweight='bold',ha='left')
fig.text(.025,.93,'Malignant epithelial cells | 24 regulons | 21 patients contributing 37 patient-program profiles',fontsize=10,color='#53616D')
fig.text(.20,.061,'Median AUCell per patient-program (>=10 cells), then median across patients; rows standardized within each algorithm.',fontsize=8.5)
fig.text(.20,.041,'C1-C4: display clusters from the original row dendrogram (k=4); labels show regulon counts. No new cell subtypes or statistical tests.',fontsize=8.5)
name='03_B_cNMF_patient_balanced_regulon_heatmap_annotations_v2'
for ext in ['png','pdf']:fig.savefig(w/'figures'/(name+'.'+ext),dpi=360)
plt.close(fig)
# Save all labels and values used to create the revised panel for reproducibility.
pd.DataFrame(Zall[leaves],index=regs,columns=[f'{fac}_{method}' for fac in factors for method in methods]).to_csv(w/'data/HEATMAP_VALUES_UNCHANGED.tsv',sep='\t')
pd.DataFrame({'regulon':regs,'display_cluster':['C'+str(c) for c in clusters]}).to_csv(q/'ROW_CLUSTER_ANNOTATIONS.tsv',sep='\t',index=False)
pd.DataFrame([{'column':2*i+j,'cell_type':'Malignant epithelial','cnmf_program':fac,'algorithm':method,'n_patients':int(nprof[fac])} for i,fac in enumerate(factors) for j,method in enumerate(methods)]).to_csv(q/'COLUMN_ANNOTATIONS.tsv',sep='\t',index=False)
validation={'status':'PASS','regulons':24,'unique_patients':21,'patient_program_profiles':37,'patient_counts':nprof.to_dict(),'same_regulons_and_row_order':True,'same_column_order':True,'same_values_max_errors':errors,'row_cluster_counts':sizes,'row_cluster_rule':'original average/Euclidean dendrogram; maxclust k=4; labels in display order','cell_type':'Malignant epithelial','cell_type_note':'F1-F7 are cNMF programs, not independent cell types','new_statistical_tests':False,'models_rerun':False}
(q/'VALIDATION.json').write_text(json.dumps(validation,indent=2));print(json.dumps(validation),flush=True)
manifest=[]
for p in paths+[Path(__file__),w.parent/'local_records_migration_20260908_150455/references/image12.png']:
 manifest.append({'path':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
(q/'SOURCE_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
(q/'SOFTWARE_VERSIONS.json').write_text(json.dumps({'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__,'matplotlib':matplotlib.__version__},indent=2))
pdf=w/'figures'/(name+'.pdf');info=subprocess.run(['pdfinfo',str(pdf)],capture_output=True,text=True,check=True).stdout
assert any(x.startswith('Pages:') and x.split(':')[1].strip()=='1' for x in info.splitlines())
subprocess.run(['pdftoppm','-r','125','-singlefile','-png',str(pdf),str(q/'B_v2_pdf')],check=True)
(q/'PDF_INFO.txt').write_text(info)
print('B_ANNOTATIONS_RENDERED',flush=True)
