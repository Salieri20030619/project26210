from pathlib import Path
import json
A=Path('/mnt/BioDisk/wangmingyuan/project26210_phase98_public_paired_metabolomics_v1')
B=Path('/mnt/BioDisk/wangmingyuan/project26210_phase77_mta_perturbation_signature_v1')
for sub in ['analysis_v2','analysis_phase99_v1','analysis_validation_v4','analysis_phase99_validation_v1']:
 p=A/sub;print('DIR',str(p));print([x.name for x in p.iterdir() if x.is_file()])
for name in ['phase77_weighted_signatures_v1.tsv','phase77_signature_counts_v1.tsv']:
 p=B/'training_v1'/name;print('FILE',str(p));print('\n'.join(p.read_text().splitlines()[:5]))
