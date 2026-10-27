from pathlib import Path
import json,shutil,datetime
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=r/'phase198_fig2_contrast_corr_revision_20260918_v5'
now=datetime.datetime.now().isoformat()
for n in ['PROJECT_STATUS.md','DECISIONS.md','CODEX_HANDOFF.md','RESULTS_MANIFEST.tsv']:
 assert not (d/'record_backup'/n).exists()
 shutil.copy2(r/n,d/'record_backup'/n)
note=f"""
## phase198 | Fig1暂定稿，Fig2色阶、去重相关图及内部编号修订 | {now}
用户明确Fig1可暂时定稿，保留Main_Figure_1_selected_v3，记录原PDF SHA256。仅修订Fig2，Fig3及正文不变。
A保留Top20原排序、原患者等权AUCell及点面积，颜色改为每通路两组全部合格程序均值共同计算的行Z分数；色标明确标注。颜色用于通路内相对差异，不作跨通路绝对值比较。全局非线性色阶草稿仍在本阶段脚本中留痕。
原C生化示意图从主图移出，逐字节原图单独交付。新C采用用户指定相关图：左下保留20患者中位Spearman rho，右上换为原逐患者rho>0人数/20与比例带，去除重复相关系数及恒等对角圆；原四通路P/FDR不变。10对相关及一致性与原数据核对通过。D复用原密度代码和冻结数据，仅删除内部A/B。
3张PDF实际渲染逐图目视核验，Top20原值/点面积/排序和行Z计算通过；无新推断检验。PDF、SVG、600dpiPNG、JPG及中文说明已交付，本地13/13哈希匹配。
服务器：{d}/delivery；本地：C:/Users/31937/Downloads/26210/28_Fig2色阶与相关图修订_v5。
当前断点：Fig2 v5待用户审阅；Fig1暂定稿，后续不自动修改；主图继续独立交付，未插入Phase191正文。无活跃绘图任务。
"""
for n in ['PROJECT_STATUS.md','DECISIONS.md','CODEX_HANDOFF.md']:
 with (r/n).open('a',encoding='utf8') as f:f.write(note)
with (r/'RESULTS_MANIFEST.tsv').open('a') as f:f.write(f'\nphase198_fig2_v5\tphase198\tfigures\t{d}/delivery\t{now}\tPASS\t3 PDF visual review; raw-value and consistency checks; local hash match\tFig1 provisionally accepted; Fig2 only; standalone schematic\n')
p=d/'STATUS.json';s=json.loads(p.read_text());s.update(status='complete',local_delivery='C:/Users/31937/Downloads/26210/28_Fig2色阶与相关图修订_v5',local_hash_matches=13,completed=now);p.write_text(json.dumps(s,indent=2))
print('PHASE198_COMPLETE')
