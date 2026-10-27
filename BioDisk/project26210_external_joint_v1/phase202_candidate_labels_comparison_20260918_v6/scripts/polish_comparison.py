from pathlib import Path
p=Path('phase202_candidate_labels_comparison_20260918_v6/scripts/compare.py');s=p.read_text().replace('[14,110,120]','[14,120,145]');p.write_text(s)
head="from pathlib import Path\nimport pandas as pd,numpy as np,matplotlib\nmatplotlib.use('Agg')\nimport matplotlib.pyplot as plt\nr=Path.cwd();d=r/'phase202_candidate_labels_comparison_20260918_v6'\ns=pd.read_csv(d/'tables/CANDIDATE_VS_WHOLE_SUMMARY.tsv',sep='\\t')\n"
exec(compile(head+s[s.index('plt.rcParams.update'):],'polish','exec'))
import shutil
shutil.copy2(d/'figures/Candidate_states_vs_whole_macrophages_v6.pdf',r/'phase201_three_candidates_C4_main_report_20260918_v4/figures/Figure9_candidate_states_vs_whole.pdf')
