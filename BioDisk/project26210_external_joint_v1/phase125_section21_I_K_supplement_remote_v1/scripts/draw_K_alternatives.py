from pathlib import Path
import pandas as pd,numpy as np,anndata as ad,json,hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
w=Path(__file__).resolve().parents[1];b=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1')
# Fixed shared panel, selected from earlier 2.1/2.2 figures before looking at group differences.
genes=['TYMP','TSPAN1','CLDN4','AOC1','NAPSA','NDUFA4L2','CA9','MET','SLC17A3','VEGFA']
calls=pd.read_csv(b/'03_malignancy_cnv/malignancy_cnv_summary_v1/tumor_epithelial_malignancy_calls_v1.tsv.gz',sep='\t');calls=calls[calls.malignancy_status_v1.isin(['malignant_consensus','malignant_infercnv_only','non_malignant_consensus'])].copy();calls['CNV_group']=np.where(calls.malignancy_status_v1=='non_malignant_consensus','Non-malignant','Malignant');assert len(calls)==12448 and calls.exact_cell_id.is_unique
metadata=pd.read_csv(b/'07_pyscenic_v1/input/cell_metadata_v1.tsv.gz',sep='\t');contr=pd.read_csv(b/'04_subclones_tymp/cnv_clone_tymp_cnmf_association_v1/patient_factor_high_clone_tymp_contrasts_v1.tsv',sep='\t');contr=contr[(contr.distance_threshold==.25)&(contr.factor=='F1')];assert len(contr)==14
high=contr.set_index('infercnv_unit').factor_high_clone_id;metadata['selected_high_clone']=metadata.curated_patient_id.map(high);meta=metadata[metadata.selected_high_clone.notna() & metadata.clone_size_class_v1.eq('primary_clone_ge20_v1')].copy();meta['clone_group']=np.where(meta.formal_inferred_cnv_clone_v1==meta.selected_high_clone,'F1-high clone','Other clones')
clone_counts=meta.groupby(['study','curated_patient_id','clone_group']).size().rename('cells').reset_index();assert clone_counts.curated_patient_id.nunique()==14 and len(clone_counts)==28 and clone_counts.cells.min()>=20
oldqc=pd.read_csv(b/'05_cnmf/cnv_clone_program_function_v1/pseudobulk_sample_qc_v1.tsv',sep='\t');mapping={'F1_high_primary_clone':'F1-high clone','other_primary_clones':'Other clones'};oldqc['clone_group']=oldqc.comparison_group.map(mapping);match=clone_counts.merge(oldqc[['curated_patient_id','clone_group','n_cells']],on=['curated_patient_id','clone_group'],validate='one_to_one');assert (match.cells==match.n_cells).all()
ct=calls.groupby(['study','curated_patient_id','CNV_group']).size().unstack(fill_value=0);eligible=ct[(ct.Malignant>=20)&(ct['Non-malignant']>=20)].index.get_level_values('curated_patient_id').tolist();assert len(eligible)==3
ct['included_paired_radar']=(ct.Malignant>=20)&(ct['Non-malignant']>=20);ct.to_csv(w/'qc/CNV_ALL_PATIENT_GROUP_COUNTS.tsv',sep='\t');clone_counts.to_csv(w/'qc/CLONE_PATIENT_GROUP_COUNTS.tsv',sep='\t',index=False)
a=ad.read_h5ad(b/'02_integration_annotation/cell_annotation_v1/part1_scvi_integrated_cell_annotated_v1.h5ad',backed='r');idx=a.obs_names.get_indexer(calls.exact_cell_id);assert (idx>=0).all();order=np.argsort(idx);xx=a[idx[order],genes].X.toarray();values=np.empty_like(xx);values[order]=xx;a.file.close();assert np.isfinite(values).all();expr=pd.DataFrame(values,index=calls.exact_cell_id,columns=genes);expr.to_csv(w/'data/marker_expression_server_only.tsv.gz',sep='\t')
# Equal patient weighting: mean log1p RNA within each patient/group, then median across paired patients.
cl=meta[['exact_cell_id','study','curated_patient_id','clone_group']].merge(expr,left_on='exact_cell_id',right_index=True,validate='one_to_one');clmean=cl.groupby(['study','curated_patient_id','clone_group'])[genes].mean();cn=calls[calls.curated_patient_id.isin(eligible)][['exact_cell_id','study','curated_patient_id','CNV_group']].merge(expr,left_on='exact_cell_id',right_index=True,validate='one_to_one');cnmean=cn.groupby(['study','curated_patient_id','CNV_group'])[genes].mean()
clmedian=clmean.groupby('clone_group').median().reindex(['F1-high clone','Other clones']);cnmedian=cnmean.groupby('CNV_group').median().reindex(['Malignant','Non-malignant']);clmean.to_csv(w/'qc/K1_PATIENT_MEANS.tsv',sep='\t');cnmean.to_csv(w/'qc/K2_PATIENT_MEANS.tsv',sep='\t');clmedian.to_csv(w/'qc/K1_GROUP_MEDIANS.tsv',sep='\t');cnmedian.to_csv(w/'qc/K2_GROUP_MEDIANS.tsv',sep='\t')
# Independent median check by directly aggregating each patient slice.
for frame,groupfield,expected in [(cl,'clone_group',clmedian),(cn,'CNV_group',cnmedian)]:
 for group in expected.index:
  v=frame[frame[groupfield]==group];direct=np.array([p[genes].to_numpy().mean(axis=0) for _,p in v.groupby('curated_patient_id')]);assert np.max(np.abs(np.median(direct,axis=0)-expected.loc[group].to_numpy()))<1e-5
maxrad=float(np.ceil(max(clmedian.to_numpy().max(),cnmedian.to_numpy().max())*2)/2);angles=np.linspace(0,2*np.pi,len(genes),endpoint=False);closed=np.r_[angles,angles[0]];colors=['#C95D6A','#4479A9']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'axes.linewidth':.7})
def draw(df,title,n,name,subtitle):
 fig=plt.figure(figsize=(7.1,6.1));ax=fig.add_axes([.14,.20,.66,.68],projection='polar');ax.set_theta_offset(np.pi/2);ax.set_theta_direction(-1);ax.set_xticks(angles,genes,fontsize=10,fontstyle='italic');ax.tick_params(axis='x',pad=9);ax.set_ylim(0,maxrad);ticks=np.linspace(0,maxrad,5)[1:];ax.set_yticks(ticks,[f'{t:g}' for t in ticks],fontsize=8,color='#777777');ax.set_rlabel_position(18);ax.grid(color='#C8CCD1',lw=.6);ax.spines['polar'].set_color('#9AA0A6');ax.spines['polar'].set_linewidth(.7)
 for (group,row),color in zip(df.iterrows(),colors):
  vals=row[genes].to_numpy();ax.plot(closed,np.r_[vals,vals[0]],color=color,lw=1.9,marker='o',markersize=4.2,label=group);ax.fill(closed,np.r_[vals,vals[0]],color=color,alpha=.10)
 fig.text(.025,.96,'K',fontsize=19,fontweight='bold',va='top');fig.text(.14,.956,title,fontsize=12,fontweight='bold',va='top');fig.text(.14,.912,f'n = {n} paired patients'+subtitle,fontsize=9,color='#555555');fig.legend(*ax.get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.49,.086),ncol=2,frameon=False,fontsize=10);fig.text(.14,.035,'Radius: median of patient-level mean log1p RNA expression\nSame ten genes and radial scale in both versions.',fontsize=8.5,color='#555555')
 for ext in ['png','pdf']:fig.savefig(w/'figures'/(name+'.'+ext),dpi=350,bbox_inches='tight',pad_inches=.12)
 plt.close(fig)
draw(clmedian,'F1-high clone vs other clones',14,'K1_F1_high_vs_other_clones_radar_v1','')
draw(cnmedian,'Malignant vs non-malignant epithelial cells',3,'K2_malignant_vs_nonmalignant_radar_v1','; exploratory')
result={'status':'PASS','gene_panel':genes,'gene_selection':'fixed from previous 2.1/2.2 panels; identical between both radars, no effect-based selection','K1_patients':14,'K1_cells':len(cl),'K1_formal_sample_counts_match':True,'K2_patients':3,'K2_cells':len(cn),'K2_scope':'paired patients with at least 20 cells in each CNV group; GSE207493 only','radial_max':maxrad,'summary':'median across paired patients of within-patient arithmetic mean log1p RNA; same units and axes for both','new_P_or_FDR':False,'independent_summary_recheck':'PASS','I':'exact copy of previously delivered Phase122 I marker violin/points; n=11950 malignant and498 non-malignant cells'}
(w/'qc/VALIDATION.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))
