from pathlib import Path
import json,shutil,datetime,hashlib
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=r/'phase199_fig3_marker_GSEA_20260918_v5'
now=datetime.datetime.now().isoformat()
for n in ['PROJECT_STATUS.md','DECISIONS.md','CODEX_HANDOFF.md','RESULTS_MANIFEST.tsv']:
 assert not (d/'record_backup'/n).exists()
 shutil.copy2(r/n,d/'record_backup'/n)
note=f"""
## phase199 | Fig2暂定稿；Fig3保留marker热图并以GSEA替换比例带 | {now}
用户明确Fig2暂时定稿；Fig1/2保持现有版本，仅修改Fig3。A改为用户指定03_macrophage_marker_heatmap.pdf原矢量热图（本地/服务器SHA256一致），移出旧基因–功能网络；C移出原C4细胞/RNA比例带，换为参考r.09_GSEA.r三层形式的四通路GSEA。B去掉内部A–D，四张局部均值矩阵与原值完全一致；D/E原轨迹和动态保留。
新GSEA使用26患者C4对其他成熟巨噬的原配对差值，12,048合格基因，连续paired t统计量排序无并列。299条KEGG同家族BH，无非估计条目。四条绘图通路计算前固定，去除TYMP/SPP1敏感性同方向：胆固醇NES3.064/FDR8.64e-12，溶酶体3.732/4.53e-34，吞噬体2.054/1.01e-5，抗原呈递-1.330/P0.0861/FDR0.153。抗原呈递本次GSEA未显著，原吞噬体整体程序患者级P0.542/FDR0.563仍完整保留；GSEA不能替代这项患者级结论。
最初signed Wilcoxon P排序存在大量并列，原输出完整保存在wilcoxon_rank_sensitivity；其中抗原呈递显著性与连续主排序不同，已记入说明，不作稳健GSEA显著结论。新GSEA为给定患者排序的基因集富集，不等同患者标签置换或独立C4验证。
独立复算8条ES及两版299通路BH、配对t排序、热图源哈希通过；两张最终PDF实际渲染视觉检查通过，GSEA图例/底部统计颜色一致。PDF/SVG/600dpiPNG/JPG及中文说明已交付，本地9/9哈希匹配。
用户新授权SSH config中已配置服务器与本地可随时双向传文件，已写入本地全局C:/Users/31937/.codex/AGENTS.md及本项目AGENTS.md；原参考脚本按新授权成功归档。
服务器：{d}/delivery；本地：C:/Users/31937/Downloads/26210/29_Fig3热图与GSEA修订_v5。
当前断点：Fig3 v5待用户审阅；Fig1/2暂时定稿，Phase191正文未改，尚未插入最新三主图；无活跃分析/绘图任务。
"""
for n in ['PROJECT_STATUS.md','DECISIONS.md','CODEX_HANDOFF.md']:
 with (r/n).open('a',encoding='utf8') as f:f.write(note)
with (r/'RESULTS_MANIFEST.tsv').open('a') as f:f.write(f'\nphase199_fig3_marker_gsea_v5\tphase199\tfigures_and_GSEA\t{d}/delivery\t{now}\tPASS\t26 paired patients; 299 KEGG BH; 8 independent ES; 2 PDF visual QA; local hash match\tFig1/2 accepted; Fig3 revised\n')
p=d/'STATUS.json';s=json.loads(p.read_text());s.update(status='complete',local_delivery='C:/Users/31937/Downloads/26210/29_Fig3热图与GSEA修订_v5',local_hash_matches=9,completed=now);p.write_text(json.dumps(s,indent=2))
(d/'qc/LOCAL_DELIVERY_VALIDATION.json').write_text(json.dumps({'files':9,'hash_mismatches':0,'checked_on':'Windows destination','checked_at':now},indent=2))
print('PHASE199_COMPLETE')
