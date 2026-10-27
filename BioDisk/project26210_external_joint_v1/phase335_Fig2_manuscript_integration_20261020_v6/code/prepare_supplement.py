from pathlib import Path
import json,shutil,hashlib,tarfile
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=r/'phase336_supporting_materials_20261020_v1';base=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/05_cnmf')
for x in ['inputs','code','qc','delivery']:(d/x).mkdir(parents=True,exist_ok=True)
specs={
'S3_genes.pdf':r/'phase333_Fig2A_candidates_20261019_v1/delivery/Candidate_3_cNMF_representative_genes.pdf',
'D01_gene_loadings.tsv':r/'phase333_Fig2A_candidates_20261019_v1/source_data/frozen_gene_relative_loadings.tsv',
'D02_rank_qc.tsv':base/'full_ranks4_8_v1/rank_qc_v1.tsv',
'D02_study_counts.tsv':r/'phase123_section2_2_reference_redraw_remote_v1/qc/program_study_cell_counts.tsv',
'D03_clone_DE.tsv.gz':base/'cnv_clone_program_function_v1/F1_high_clone_vs_other_paired_pseudobulk_DE_v1.tsv.gz',
'D03_patient_deltas.tsv.gz':base/'cnv_clone_program_function_v1/top_DE_patient_logCPM_deltas_v1.tsv.gz',
'D03_pathway_ORA.tsv':base/'cnv_clone_program_function_v1/F1_high_clone_DE_pathway_ORA_v1.tsv',
'D06_full_regulon_targets.tsv':r/'phase127_section2_4_pyscenic_redraw_remote_v1/data/CANDIDATE_REGULON_TARGETS.tsv',
'Q_cnmf_validation.json':r/'phase123_section2_2_reference_redraw_remote_v1/qc/VALIDATION.json',
'Q_cnmf_independent_validation.json':r/'phase123_section2_2_reference_redraw_remote_v1/qc/INDEPENDENT_VALIDATION.json',
'Q_immune_validation.json':r/'phase321_TYMP_SPP1_conditional_gene_robustness_20261005_v1/qc/FINAL_VALIDATION.json',
'Q_spatial_validation.json':r/'phase219_spatial_reference_identifiability_20260920_v1/FINAL_VALIDATION.json',
'Q_spatial_notes.md':r/'phase219_spatial_reference_identifiability_20260920_v1/REFERENCE_AND_PROGRAM_CHECK_CN.md',
}
for p in (r/'phase321_TYMP_SPP1_conditional_gene_robustness_20261005_v1/tables').iterdir():
 if p.is_file():specs['D04_'+p.name]=p
for name in ['PROGRAM_CHECK_SUMMARY.tsv','PATIENT_EFFECTS.tsv','SECTION_EFFECTS.tsv','IDENTIFIABILITY_RESULTS.tsv','REFERENCE_COVERAGE.tsv']:
 specs['D05_'+name]=r/'phase219_spatial_reference_identifiability_20260920_v1/delivery'/name
manifest=[]
for name,p in specs.items():
 assert p.exists(),p
 dst=d/'inputs'/name;shutil.copy2(p,dst)
 manifest.append({'name':name,'source':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
(d/'qc/SOURCE_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
with tarfile.open(d/'inputs_bundle.tar.gz','w:gz') as t:
 t.add(d/'inputs',arcname='inputs');t.add(d/'qc/SOURCE_MANIFEST.json',arcname='qc/SOURCE_MANIFEST.json')
# Close out manuscript stage before supplement work.
m=r/'phase335_Fig2_manuscript_integration_20261020_v6';v=json.loads((m/'qc/PAGE_AND_LINK_CHECKS.json').read_text())
assert len(v['identical_pages'])==47 and len(v['links'])==79 and not v['invalid_links']
(m/'qc/COMPLETE.json').write_text(json.dumps({'status':'PASS','pages':54,'tables':41,'changed_pages':v['changed_pages'],'pixel_identical_pages':47,'valid_internal_links':79,'visual_review':'All seven changed pages reviewed; no clipping or overlap; other 47 pages pixel identical to v5','pdf_sha256':hashlib.sha256((m/'delivery/26210_manuscript_Fig2_updated_v6.pdf').read_bytes()).hexdigest()},indent=2),encoding='utf-8')
print('Supplement frozen inputs:',len(manifest),'bundle bytes',(d/'inputs_bundle.tar.gz').stat().st_size)
