from pathlib import Path
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=r/'phase197_top20_network_dynamics_20260917_v4'
p=d/'scripts/redraw.py';s=p.read_text()
a=s.index("# Compact source-share");b=s.index("# Curves replace",a)
replacement="""# Source-share annotation: two compact side-by-side composition strips.
co=pd.read_csv(r/'phase177_L2_myeloid_state_integration_20260916_v1/source_tables/c4_contribution.tsv',sep='\\t')
fr=co[['C4_fraction','C4_SPP1_RNA_share']].median().values
fig,axs=plt.subplots(1,2,figsize=(12.5,1.65))
fig.subplots_adjust(left=.01,right=.99,bottom=.22,top=.94,wspace=.12)
for ax,value,label in zip(axs,fr,['Macrophage cell abundance','Captured SPP1 RNA contribution']):
 ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
 ax.text(.02,.86,label,fontsize=10,weight='bold')
 ax.add_patch(Rectangle((.02,.33),.76,.27,fc='#E1E5E5',ec='none'))
 ax.add_patch(Rectangle((.02,.33),.76*value,.27,fc=purple,ec='none'))
 ax.text(.82,.465,f'{value:.1%}',ha='left',va='center',fontsize=18,weight='bold',color=purple)
fig.text(.5,.03,'Purple: C4 | Gray: other macrophages | Patient medians, n = 16',ha='center',fontsize=8,color='#555555')
panel(fig,'C4_source_share')
"""
s=s[:a]+replacement+s[b:]
s=s.replace("fig,axs=plt.subplots(1,5,figsize=(12.5,3.0),layout='constrained');records=[]","fig,axs=plt.subplots(2,3,figsize=(6.8,5.2),layout='constrained');axs=axs.ravel();axs[-1].axis('off');records=[]")
s=s.replace("fig.legend(handles=","axs[-1].legend(handles=").replace("loc='outside upper center',ncol=2","loc='center',ncol=1")
p.write_text(s)
# Assemble only changed Fig2 and Fig3; prior Fig1 remains unchanged.
s=(r/'phase195_mainfigures_compact_20260917_v2/scripts/assemble_compact.py').read_text()
s=s.replace("d=r/'phase195_mainfigures_compact_20260917_v2'","d=r/'phase197_top20_network_dynamics_20260917_v4';old195=r/'phase195_mainfigures_compact_20260917_v2'")
s=s.replace("d/'panels/network_compact.pdf'","old195/'panels/network_compact.pdf'").replace("(d/'panels'/f'{name}.pdf'","(old195/'panels'/f'{name}.pdf'")
s=s.replace("# Reconstruct tightly", """for key in ['top20_pathways','gene_function_network','C4_source_share','expression_dynamics_curves']:
 pan[key]=(d/'panels'/f'{key}.pdf',[0,0,1,1])
pan['epi']=(next((r/'phase191_discovery_order_report_20260917_v1/figures').glob('Figure7_*.pdf')),[.27,.027,.68,.303])
# Reconstruct tightly""")
start=s.index('figs=[]');end=s.index('records=[]',start)
new="""figs=[]
yd=24+max(height('epi',71),height('2B',72))+12
ye=max(24+height('top20_pathways',165),yd+height('2E',152))+12
pp=[('A','Top 20 TYMP-associated metabolic pathways','top20_pathways',8,24,165),
('B','Malignant epithelial trajectory','epi',185,24,71),
('C','MTA biochemical context','2B',260,24,72),
('D','Discovery and localization of C4','2E',180,yd,152),
('E','Measured MTA and tissue expression','2G',8,ye,324)]
figs.append(('Main_Figure_2_top20_v4','Metabolic context and discovery of an SPP1-enriched state',pp))
yc=24+max(height('gene_function_network',155),height('3B',165))+12
yb=yc+height('C4_source_share',324)+12
pp=[('A','Functional genes and pathway membership','gene_function_network',8,24,155),
('B','Functional expression contours','3B',167,24,165),
('C','C4 abundance and contribution to SPP1 RNA','C4_source_share',8,yc,324),
('D','Original partitioned myeloid trajectories','3C',8,yb,155),
('E','Expression dynamics across myeloid states','expression_dynamics_curves',167,yb,165)]
figs.append(('Main_Figure_3_network_curves_v4','C4 functional identity, SPP1 contribution and state ordering',pp))
"""
s=s[:start]+new+s[end:]
(d/'scripts/assemble.py').write_text(s)
