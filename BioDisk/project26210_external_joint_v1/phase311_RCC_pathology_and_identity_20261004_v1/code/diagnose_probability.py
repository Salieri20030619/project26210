from pathlib import Path
import pandas as pd
p=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase311_RCC_pathology_and_identity_20261004_v1')
f=p/'resources/proseg_region6_v1/transcript-metadata.csv.gz'
d=pd.read_csv(f,nrows=250000)
print(d.columns.tolist());print(d.head(3).to_string(index=False));print(d.iloc[:,11].describe().to_string());bad=d[~d.iloc[:,11].between(0,1)];print('BAD',len(bad));print(bad.head(10).to_string(index=False))
