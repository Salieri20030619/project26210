from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
from scipy.stats import spearmanr
w=Path(__file__).resolve().parents[1];r=w.parent
b=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1')
read=lambda p:pd.read_csv(p,sep='\t')
checks=[];metrics={}
def ck(name,value):
 checks.append({'check':name,'pass':bool(value)})
 if not value:raise AssertionError(name)
meta=json.loads((w/'qc/FIGURE_METADATA.json').read_text())
bounds=json.loads((w/'qc/RENDER_TEXT_BOUNDS.json').read_text())
ck('A4_physical_size',meta['width_mm']<=210 and meta['height_mm']<=297)
ck('11_panels_A_to_K',[x['panel']for x in meta['panels']]==list('ABCDEFGHIJK'))
ck('75_library_scope',meta['libraries']==75)
ck('all_text_inside_page',not bounds['outside_page'])
ck('font_at_least_6_5pt',bounds['minimum_font_size']>=6.5)
c=read(w/'source_data/CNV_UMAP_POINTS.tsv.gz')
ck('CNV_all_cells_retained',len(c)==38670 and c.cell_id.is_unique)
ck('CNV_classes',c['class'].value_counts().to_dict()=={'Unresolved':26222,'Malignant':11950,'Nonmalignant':498})
pt=read(r/'phase123_section2_2_reference_redraw_remote_v1/qc/F1_TYMP_patients.tsv')
s=read(r/'phase123_section2_2_reference_redraw_remote_v1/qc/F1_TYMP_formal_summary.tsv').iloc[0]
ck('F1_TYMP_median_and_n',len(pt)==20 and (pt.effect>0).sum()==17 and abs(pt.effect.median()-s.median_effect)<1e-12)
rp=read(r/'phase127_section2_4_pyscenic_redraw_remote_v1/qc/INDEPENDENT_CANDIDATE_PATIENT_RECHECK.tsv')
rs=read(w/'source_data/REGULON_PANEL_STATISTICS.tsv')
for v in rs.itertuples():
 a=rp[(rp.method==v.method)&(rp.tf==v.tf)&rp.continuous_eligible]
 ck('reg_'+v.method+'_'+v.tf,len(a)==20 and abs(a.rho.median()-v.median_patient_effect)<1e-12)
act=pd.read_csv(r/'phase128_section2_5_metaviper_redraw_remote_v1/data/SELECTED_METAVIPER_ACTIVITY.tsv',sep='\t',index_col=0)
m=pd.read_csv(b/'05_cnmf/input_v1/cnmf_cell_metadata_v1.tsv.gz',sep='\t',index_col=0).loc[act.index]
u=pd.read_csv(b/'05_cnmf/full_ranks4_8_v1/consensus_usage_rank7_v1.tsv.gz',sep='\t',index_col=0).loc[act.index]
fs=read(w/'source_data/F1_ACTIVITY_PANEL_STATISTICS.tsv')
for v in fs.itertuples():
 vals=[]
 for patient,ix in m.groupby('curated_patient_id').groups.items():vals.append(spearmanr(u.loc[ix,'F1'],act.loc[ix,v.protein]).statistic)
 ck('F1_density_patient_rho_'+v.protein,len(vals)==20 and abs(np.median(vals)-v.median_patient_effect)<1e-12)
gp=read(r/'phase133_sections2_5_2_7_targeted_support_remote_v1/data/TYMP_ACTIVITY_PATIENT_EFFECTS.tsv')
gs=read(w/'source_data/TYMP_ACTIVITY_PANEL_STATISTICS.tsv')
for v in gs.itertuples():
 d=gp[(gp.target==v.target)&(gp.model==v.model)]
 ck('TYMP_activity_'+v.target+'_'+v.model,len(d)==20 and abs(d.effect.median()-v.median_effect)<1e-12 and int((d.effect>0).sum())==v.positive)
hp=read(r/'phase129_sections2_6_2_7_redraw_remote_v1/qc/original_four_pathway_AUCell_patient_deltas_v1.tsv')
hs=read(w/'source_data/METABOLIC_PANEL_STATISTICS.tsv')
ck('all8_metabolic_endpoints',len(hs)==8)
for v in hs.itertuples():
 d=hp[(hp.pathway==v.pathway)&(hp.contrast==v.contrast)]
 ck('metabolic_'+v.pathway+'_'+v.contrast,len(d)==14 and abs(d.delta.median()-v.median_patient_delta)<1e-12)
kp=read(r/'phase133_sections2_5_2_7_targeted_support_remote_v1/data/TERMINAL_ENZYME_PATIENT_EFFECTS.tsv')
ks=read(w/'source_data/TERMINAL_PANEL_STATISTICS.tsv')
ck('39_terminal_endpoints',len(ks)==39 and kp.shape[0]==780)
err=0
for v in ks.itertuples():
 d=kp[(kp.gene==v.gene)&(kp.method==v.method)]
 ck('terminal_'+v.method+'_'+v.gene,len(d)==20 and abs(d.effect.median()-v.median_effect)<1e-12 and (d.early_cells>=5).all() and (d.late_cells>=5).all())
 err=max(err,float(np.max(np.abs((d.late_log2_CPM1-d.early_log2_CPM1)-d.effect))))
ck('terminal_patient_difference_identity',err<1e-12)
ck('nominal_only_star_semantics',(ks.FDR>=.05).all() and (ks.P<.05).sum()==4)
heat=np.load(r/'phase132_section2_7_gradient_heatmaps_remote_v1/data/GRADIENT_HEATMAP_VALUES.npz')
labels=read(w/'source_data/HEATMAP_LABEL_POSITIONS.tsv');ordered=heat['genes'][heat['leaves']]
ck('2511_genes_all_rows',len(ordered)==2511 and heat['Monocle3'].shape==(2511,200))
ck('12_correct_label_rows',len(labels)==12 and all(ordered[int(v.matrix_row)]==v.gene for v in labels.itertuples()))
ck('label_spacing',np.diff(labels.label_y_mm.to_numpy()).min()>2.8)
ck('cluster_sizes',[(heat['clusters']==x).sum()for x in [1,2,3,4]]==[408,741,340,1022])
legend=(w/'figures/Figure1_TYMP_integrated_A4_v1_legend_zh.txt').read_text()
ck('11_legends_with_axes_results',legend.count('横轴：')==11 and legend.count('最主要结果：')==11)
manifest=json.loads((w/'qc/SOURCE_MANIFEST.json').read_text())
ck('source_hashes_unchanged',all(hashlib.sha256(Path(x['path']).read_bytes()).hexdigest()==x['sha256']for x in manifest))
metrics['maximum_terminal_delta_error']=err
(w/'qc/INDEPENDENT_VALIDATION.json').write_text(json.dumps({'status':'PASS','checks':checks,'metrics':metrics},indent=2))
print('PASS',len(checks),'checks',metrics)
