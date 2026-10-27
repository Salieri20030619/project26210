from pathlib import Path
import pandas as pd,json
A=Path('/mnt/BioDisk/wangmingyuan/project26210_phase98_public_paired_metabolomics_v1')
for f in ['analysis_v2/PHASE98_PAIRED_MTA_SUMMARY_v2.md','analysis_v2/paired_sample_counts_v2.tsv','analysis_v2/cross_cohort_meta_analysis_v2.tsv','analysis_v2/cohort_associations_v2.tsv','analysis_phase99_v1/meta_tymp_mtap_associations_v1.tsv']:
 p=A/f;print('SOURCE',f);print(p.read_text()[:13500])
