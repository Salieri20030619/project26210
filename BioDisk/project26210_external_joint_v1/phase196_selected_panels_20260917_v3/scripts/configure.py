from pathlib import Path
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=r/'phase196_selected_panels_20260917_v3';p=d/'scripts/draw_column.py'
p.write_text(p.read_text().replace("dot_facets(paths,'85_pathways_column','','')","dot_facets(paths,'85_pathways_column','')"))
s=(r/'phase195_mainfigures_compact_20260917_v2/scripts/assemble_compact.py').read_text()
s=s.replace("d=r/'phase195_mainfigures_compact_20260917_v2'","d=r/'phase196_selected_panels_20260917_v3';old195=r/'phase195_mainfigures_compact_20260917_v2'")
s=s.replace("d/'panels/network_compact.pdf'","old195/'panels/network_compact.pdf'").replace("(d/'panels'/f'{name}.pdf'","(old195/'panels'/f'{name}.pdf'")
s=s.replace("# Reconstruct tightly", """pan['1D']=(d/'inputs/H_patient_program_expression_blocks_v2.pdf',[.05,0,1,1])
pan['2A']=(d/'panels/85_pathways_column.pdf',[0,0,1,1])
pan['epi']=(next((r/'phase191_discovery_order_report_20260917_v1/figures').glob('Figure7_*.pdf')),[.27,.027,.68,.303])
# Reconstruct tightly""")
s=s.replace("'Program genes and F1 expression'","'Patient-program expression structure'").replace("Main_Figure_1_compact_v2","Main_Figure_1_selected_v3")
start=s.index("y2=24+max(height('2A'")
end=s.index("records=[]",start)
new="""# The 85-pathway panel occupies one full-height left column.
yr=24+max(height('epi',90),height('2B',95))+12
ye=yr+height('2E',191)+12
pp=[('A','All 85 pathways by TYMP status','2A',8,24,205),
('B','Malignant epithelial trajectory','epi',221,24,90),
('C','Biochemical route to MTA','2B',317,24,95),
('D','Discovery and localization of C4','2E',221,yr,191),
('E','Measured MTA and tissue expression','2G',221,ye,191)]
figs.append(('Main_Figure_2_column_v3','Metabolic context and discovery of an SPP1-enriched state',pp))
"""
s=s[:start]+new+s[end:]
s=s.replace("W=340*mm","W=(420 if name.startswith('Main_Figure_2') else 340)*mm")
s=s.replace("print(name,'340 x',round(H/mm,1),'mm',flush=True)","print(name,round(W/mm,1),'x',round(H/mm,1),'mm',flush=True)")
(d/'scripts/assemble_selected.py').write_text(s)
