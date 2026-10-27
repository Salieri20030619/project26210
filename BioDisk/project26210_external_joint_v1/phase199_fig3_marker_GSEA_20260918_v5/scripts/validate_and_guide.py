from pathlib import Path
import numpy as np,pandas as pd,json,hashlib,sys
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=r/'phase199_fig3_marker_GSEA_20260918_v5'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(r/'phase192_sections28_212_figures_20260917_v1/figures/03_macrophage_marker_heatmap.pdf')=='96e33e28031c58474e3b7671bfc5d4c4329bdae2cf51725f3e6fa827e1e170cf'
assert sha(d/'inputs/user_reference_r.09_GSEA.r')=='a84457cb84acf3db94feaaf1802f00f9fe1c702b282f3d9c8634291bac4c1c07'
rank=pd.read_csv(d/'tables/RANKED_PATIENT_GENES.tsv',sep='\t');assert len(rank)==12048 and rank.gene.is_unique
z=np.load(r/'part2_rebuild_v1/04_markers_annotation/run_v4_macrophage_only/paired_patient_effects/cluster_4.npz',allow_pickle=True)
idx={str(g):i for i,g in enumerate(z['genes'])};v=z['effects'][:,[idx[g] for g in rank.gene]]
expect=v.mean(0)/(v.std(0,ddof=1)/np.sqrt(26))
assert np.allclose(rank.ranking_score,expect,rtol=1e-10) and rank.ranking_score.is_unique
sets={}
for line in (d/'inputs/KEGG_2021_Human.gmt').read_text().splitlines():
 fields=line.split('\t');sets[fields[0]]=set(x for x in fields[2:] if x)
primary=pd.read_csv(d/'tables/GSEA_FOUR_PROGRAMS_primary.tsv',sep='\t')
sensitivity=pd.read_csv(d/'tables/GSEA_FOUR_PROGRAMS_exclude_TYMP_SPP1.tsv',sep='\t')
checks=[]
for mode in ['primary','exclude_TYMP_SPP1']:
 a=rank if mode=='primary' else rank.loc[~rank.gene.isin(['TYMP','SPP1'])]
 full=pd.read_csv(d/'tables'/f'GSEA_KEGG_FULL_FAMILY_{mode}.tsv',sep='\t')
 assert len(full)==299 and full.pvalue.notna().all()
 eligible={k for k,s in sets.items() if 10<=len(s&set(a.gene))<=500}
 assert set(full.ID)==eligible
 order=np.argsort(full.pvalue.to_numpy());q=np.empty(len(full))
 q[order]=np.minimum.accumulate((full.pvalue.to_numpy()[order]*len(full)/np.arange(1,len(full)+1))[::-1])[::-1].clip(0,1)
 assert np.allclose(q,full['p.adjust'],rtol=1e-9,atol=1e-14)
 tab=primary if mode=='primary' else sensitivity
 for row in tab.itertuples():
  hit=np.array([g in sets[row.ID] for g in a.gene]);N=len(hit);Nh=int(hit.sum())
  weight=np.abs(a.ranking_score.to_numpy())
  walk=np.cumsum(np.where(hit,weight/weight[hit].sum(),-1/(N-Nh)))
  es=walk.max() if walk.max()>=-walk.min() else walk.min()
  assert abs(walk[-1])<1e-10 and abs(es-row.enrichmentScore)<1e-10
  assert Nh==row.setSize
  assert set(str(row.core_enrichment).split('/'))<=sets[row.ID]&set(a.gene)
  checks.append({'mode':mode,'term':row.ID,'independent_ES':float(es),'R_ES':row.enrichmentScore,'setSize':Nh})
pd.DataFrame(checks).to_csv(d/'qc/INDEPENDENT_ES_VALIDATION.tsv',sep='\t',index=False)
old=pd.read_csv(r/'part2_rebuild_v1/06_functional_programs/run_v1/C4_FOCUS_PROGRAMS.tsv',sep='\t')
old=old[(old.library=='KEGG_2021_Human')&(old.version=='complete')&old.term.isin(primary.ID)]
compare=primary.merge(old,left_on='ID',right_on='term',suffixes=('_GSEA','_patient'))
compare.to_csv(d/'tables/GSEA_AND_PRIOR_PATIENT_PROGRAMS.tsv',sep='\t',index=False)
print(compare[['ID','NES','p.adjust','median_paired_program_difference','P','FDR']].to_string(index=False))
lines=["# Fig3：marker热图与GSEA修订 v5","",
"Fig1及Fig2已按用户意见暂时定稿。本次只改Fig3，仍为独立图件，未插入正文。",
"",
"## A｜保留用户指定的marker热图",
"直接复用03_macrophage_marker_heatmap.pdf的原始矢量内容；本地指定文件与服务器原图SHA256一致。替换上一版A的基因–功能归属网络。横轴为12个冻结巨噬候选状态，纵轴为38个marker；颜色为原患者等权表达的行Z分数，紫框定位C4。该图显示SPP1/GPNMB/CCL18及脂质处理marker共同富集的状态，SPP1并非C4特异。仅裁去原章节标题以统一主图编号，未重算热图。",
"",
"## B｜功能表达等高线",
"保留上一版SPP1、胆固醇代谢、溶酶体和抗原呈递的四面板等高线。坐标是冻结转录嵌入，紫色虚线为C4边界；图为转录表达/程序分布，不是组织空间。",
"",
"## C｜GSEA替换原数量与RNA来源比例带",
"参考用户r.09_GSEA.r的gseaplot2三层形式：上层运行富集分数曲线，中层通路基因命中位置，下层排序统计量。横轴是12,048个合格基因的排序位置；左端更偏C4、右端更偏其他成熟巨噬。原比例带移出新版，其数值和历史图仍保留。",
"",
"### 数据与计算",
"沿用冻结C4与其他成熟巨噬的26名患者内配对log2(CPM+1)差值。对每个基因计算mean(delta)/(SD(delta)/sqrt(26))作为连续排序统计量；没有将细胞数作为样本量，也没有新生成基因级t检验P值。纳入原先test_eligible的全部12,048基因，不按显著性筛选。",
"基因集为已冻结KEGG_2021_Human，按symbol匹配；与排序列表交集大小限定10–500，得到299条合格通路。clusterProfiler::GSEA使用fgsea后端，exponent=1、eps=0、nPermSimple=10000；脚本设定种子20260918并启用GSEA内部seed=TRUE，软件版本和完整参数保存。每个版本在全部299通路内BH校正，而非仅对图中四条校正。",
"四条展示通路在计算前固定：胆固醇代谢、溶酶体、吞噬体、抗原呈递，不是按结果筛选的Top4。另做去除TYMP/SPP1的敏感性，保持同一基因集库与参数。",
"",
"| 通路 | 主NES | 主P | 主FDR | 去除TYMP/SPP1 NES | 敏感性FDR |",
"|---|---:|---:|---:|---:|---:|"]
for _,row in primary.iterrows():
 sr=sensitivity.loc[sensitivity.ID==row.ID].iloc[0]
 lines.append(f"| {row.ID} | {row.NES:.3f} | {row.pvalue:.4g} | {row['p.adjust']:.4g} | {sr.NES:.3f} | {sr['p.adjust']:.4g} |")
lines += ["","### 结果及证据定位",
"胆固醇代谢、溶酶体和吞噬体基因集在C4偏高端富集；抗原呈递向相反方向，但本次GSEA未达显著。去除TYMP/SPP1后四条方向均保持。",
"GSEA的P/FDR来自给定基因排序下的基因集富集检验，不是重新对患者标签做置换，也不是同一数据中C4发现的独立验证。原患者级整体程序检验是不同终点，完整保留，不能用GSEA阳性覆盖原程序阴性。",
"",
"| 原患者级程序比较 | 配对中位差 | 患者n | P | FDR |",
"|---|---:|---:|---:|---:|"]
for _,z in old.iterrows():
 lines.append(f"| {z.term} | {z.median_paired_program_difference:.4g} | {int(z.n_patients)} | {z.P:.4g} | {z.FDR:.4g} |")
lines += ["",
"特别保留：原吞噬体全程序患者级比较不支持一致增强；当前GSEA提示其基因集合在差异排序中的富集，二者不合并写成吞噬功能整体已增强。抗原呈递原患者级下降与本次GSEA负方向并列呈现，但不把本次GSEA称为显著。",
"",
"### 排序质量及留痕",
"首次参考式sign(median delta)×−log10(Wilcoxon P)排序出现大量并列；全部输出保存在wilcoxon_rank_sensitivity。因排序技术问题改用连续配对t统计量，主排序无并列。原排序下抗原呈递FDR达到阈值、连续主排序未达到，故该通路GSEA显著性依赖排序方式，不作稳健显著结论。",
"",
"## D/E｜原轨迹及表达动态",
"D保留5研究的原分区轨迹、根节点及C4定位，不重拟合或补连。主分析C4有效覆盖442/2674，敏感性1743/2674。",
"E保留Phase197表达动态曲线；原40个秩位置、主分析14患者及敏感性24患者不变，没有新增平滑或动态检验。",
"",
"## 交付与核验",
"主图仅保留A这一张色块热图。两张PDF均配SVG、600dpiPNG及JPG预览；独立GSEA图便于放大。脚本、完整299通路结果、排序列表、R对象、原始版本、敏感性及审计均保存在服务器。",
"独立Python复算8条ES与R一致，299通路BH复核通过，配对t排序由原26患者矩阵复算一致，热图与参考脚本哈希核对通过。最后以实际PDF渲染检查版式后交付。",
""]
(d/'REVISION_GUIDE_CN.md').write_text('\n'.join(lines),encoding='utf8')
(d/'qc/STATISTICAL_VALIDATION.json').write_text(json.dumps({'rank_genes':12048,'patients':26,'rank_ties':0,'full_family':299,'nonestimable':0,'independent_ES_checks':8,'BH_recomputed':True,'source_heatmap_hash_match':True,'reference_script_hash_match':True,'four_programs_preselected':True},indent=2))
print('STATISTICS_VALIDATED')
