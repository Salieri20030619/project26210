from common import *
import json
src=R/'part2_rebuild_v1/scripts'
# Reuse frozen cell-level expression scores, graph, roots and pseudotimes.
s=(src/'20_trajectory_statistics_v2.py').read_text()
a=s.index("a=ad.read_h5ad(");b=s.index("m=pd.read_csv(T/'ALL_CELL",a)
s=s[:a]+"scores=pd.read_csv(R/'07_trajectory/dynamics_v2/CELL_EXPRESSION_PROGRAM_SCORES.tsv.gz',sep='\\t',index_col=0)\n"+s[b:]
for label,members in [('U',['4','10','15']),('C10',['10']),('C15',['15'])]:
 t=s.replace("O=R/'07_trajectory/dynamics_v2'",f"O=Path('{D}/trajectory/{label}')")
 t=t.replace("c4=c.eq('4')",f"c4=c.isin({members!r})")
 t=t.replace("other=~c.isin(['4','1','Monocyte'])",f"other=~c.isin({members+['1','Monocyte']!r})")
 t=t.replace('C4_minus_',label+'_minus_').replace('C4_cells',label+'_cells')
 (D/'scripts'/f'run_trajectory_{label}.py').write_text(t)
# Frozen pySCENIC network is unchanged; only state contrasts and eligible trajectory partitions change.
u=(src/'25_new_regulon_statistics_v2.py').read_text()
u=u.replace("O=R/'08_regulons/statistics_v1'",f"O=Path('{D}/regulons/union_legacy_method')")
u=u.replace("inside=z.index[c=='4']","inside=z.index[c.isin(['4','10','15'])]")
u=u.replace("~c.isin(['4','1','Monocyte'])","~c.isin(['4','10','15','1','Monocyte'])")
u=u.replace("(z.candidate=='4').any()","z.candidate.isin(['4','10','15']).any()")
u=u.replace("R/'04_markers_annotation/run_v4_macrophage_only/full_marker_tables/cluster_4_all_genes.tsv.gz'",f"Path('{D}/core/UNION_MARKER_LEGACY_SCHEMA.tsv.gz')")
u=u.replace('C4_minus_','U_minus_').replace('C4_cells','U_cells')
pd.read_csv(D/'core/contrasts/all_U_genes.tsv.gz',sep='\t').rename(columns={'P':'patient_wilcoxon_p','median_delta':'median_patient_log2CPM_difference'}).to_csv(D/'core/UNION_MARKER_LEGACY_SCHEMA.tsv.gz',sep='\t',index=False)
(D/'scripts/run_union_regulons.py').write_text(u)
print('ADAPTED_CODE_READY')
