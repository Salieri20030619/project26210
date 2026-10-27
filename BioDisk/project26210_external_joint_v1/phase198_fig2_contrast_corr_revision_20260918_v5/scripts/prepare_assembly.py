from pathlib import Path
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=r/'phase198_fig2_contrast_corr_revision_20260918_v5'
s=(r/'phase197_top20_network_dynamics_20260917_v4/scripts/assemble.py').read_text()
s=s.replace("d=r/'phase197_top20_network_dynamics_20260917_v4'","d=r/'phase198_fig2_contrast_corr_revision_20260918_v5'")
a=s.index('# Reconstruct tightly')
s=s[:a]+"""pan['2E']=(d/'panels/C4_no_nested_letters.pdf',[.023,.12,.993,.955])
pan['corr']=(d/'panels/correlation_consistency.pdf',[0,0,1,1])
pan['tissue_compact']=(r/'phase197_top20_network_dynamics_20260917_v4/panels/tissue_compact.pdf',[0,0,1,1])
pan={k:pan[k] for k in ['top20_pathways','epi','2E','corr','tissue_compact']}
"""+s[a:]
a=s.index('figs=[]');b=s.index('records=[]',a)
s=s[:a]+"""figs=[]
yd=24+max(height('epi',59),height('corr',104))+12
ye=max(24+height('top20_pathways',150),yd+height('2E',165))+12
pp=[('A','Top 20 TYMP-associated metabolic pathways','top20_pathways',8,24,150),
('B','Epithelial trajectory','epi',167,24,59),
('C','F1-pathway associations','corr',230,24,104),
('D','Discovery and localization of C4','2E',167,yd,165),
('E','Measured MTA and tissue expression','tissue_compact',8,ye,324)]
figs.append(('Main_Figure_2_contrast_consistency_v5','Metabolic context and discovery of an SPP1-enriched state',pp))
"""+s[b:]
(d/'scripts/assemble.py').write_text(s)
