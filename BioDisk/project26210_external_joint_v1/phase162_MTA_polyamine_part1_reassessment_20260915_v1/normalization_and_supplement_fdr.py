from pathlib import Path
import pandas as pd,numpy as np,json
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');O=P/'phase162_MTA_polyamine_part1_reassessment_20260915_v1'
C=pd.read_csv(O/'PART1_TARGET_CELL_METADATA.tsv.gz',sep='\t')
diff=np.abs(C['expr:TYMP']-C.TYMP_log1p_CP10K_v1)
print('ORIGINAL_TYMP_NORMALIZATION maxabs',float(diff.max()),'median',float(diff.median()))
checks=[dict(check='Part1_original_normalization',passed=bool(diff.max()<1e-6))]
D=pd.read_csv(O/'PART1_GENE_CONTRAST_SUMMARY.tsv',sep='\t');R=pd.read_csv(O/'PART1_TYMP_GENE_CORRELATION_SUMMARY.tsv',sep='\t')
di=D.index[D.gene=='SLC16A3'];ri=R.index[R.gene=='SLC16A3']
q=multipletests(np.r_[D.loc[di,'P'],R.loc[ri,'P']],method='fdr_bh')[1]
D.loc[di,'FDR_SLC16A3_supplement5']=q[:len(di)];R.loc[ri,'FDR_SLC16A3_supplement5']=q[len(di):]
D.to_csv(O/'PART1_GENE_CONTRAST_SUMMARY.tsv',sep='\t',index=False);R.to_csv(O/'PART1_TYMP_GENE_CORRELATION_SUMMARY.tsv',sep='\t',index=False)
pd.DataFrame(checks).to_csv(O/'NORMALIZATION_CHECK.tsv',sep='\t',index=False)
assert all(r['passed'] for r in checks)
s=json.loads((O/'STATUS.json').read_text());s['total_checks']+=1;s['normalization_checks']=1;(O/'STATUS.json').write_text(json.dumps(s,indent=2))
print(D[D.gene=='SLC16A3'].to_string(index=False));print(R[R.gene=='SLC16A3'].to_string(index=False))
