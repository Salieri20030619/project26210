from common import *
import json,shutil,datetime,hashlib
now=datetime.datetime.now().isoformat()
assert json.loads((D/'qc/VALIDATION_STATUS.json').read_text())['status']=='PASS'
assert json.loads((D/'qc/ATTRIBUTION_VALIDATION.json').read_text())['status']=='PASS'
for n in ['PROJECT_STATUS.md','DECISIONS.md','CODEX_HANDOFF.md','RESULTS_MANIFEST.tsv']:
 target=D/'record_backup'/n
 assert not target.exists()
 shutil.copy2(R/n,target)
branch=R/'part2_rebuild_v1/STATUS.json'
shutil.copy2(branch,D/'record_backup/PART2_STATUS_before_union.json')
b=json.loads(branch.read_text())
b['latest_union_reassessment']={'phase':200,'path':str(D),'status':'COMPLETE_NUMERICAL_REASSESSMENT','definition':'fixed original C4+C10+C15 union; preserve individual groups','cells':4502,'paired_all_tissue_patients':28,'original_C4_preserved':True,'figures_created':False,'completed':now}
branch.write_text(json.dumps(b,ensure_ascii=False,indent=2))
note=f"""
## phase200 | C4+C10+C15联合候选群全面数值复评完成 | {now}
用户明确授权此前依赖C4定义的分析从U=C4∪C10∪C15角度重看，保留原C4及三群拆分。原成员固定，联合4502细胞，出现37患者；all-tissue配对U28/C4 26，肿瘤配对U26/C4 24。同文库上游固定16患者，群内表达共同12患者。无新图，无聚类/轨迹/网络重训，旧结果只读。
已完成：原始counts及SPP1份额、全基因marker、三库功能/ORA、KEGG GSEA、原轨迹覆盖/动态、冻结pySCENIC/regulon及三TF、固定116候选/38代谢成分/4MTA代理632个比例桥接模型、2784个群内SPP1表达终点模型、受体/STAT3/LGALS3/IL6/IL10/SLC16A3及组成调整、同患者共同背景/共同基因对比、去成分/主导患者/研究敏感性。原C4 marker/功能/regulon复现；独立数值验证2301项通过，归因追加核验通过。
主要结果：同16患者SPP1原始RNA份额中位54.59%→59.63%，细胞占比7.41%→9.76%。U28患者SPP1 25正、中位差1.821/P2.93e-5/FDR4.17e-4；共同26患者逐人U−C4效应增量中位-0.0516，10正16负/P0.980，不称普遍增强。
新线索：U糖酵解程序中位+0.1228/FDR0.00543，C4为-0.1332；G2-M/E2F转正。CollecTRI STAT3代理U26患者20正、差+0.0231/P0.0111/FDR0.0333；同24患者同nonU参照，U−C4增量+0.01856，24/24正。去C10/主导患者保持方向，去C15后糖酵解/G2-M及STAT3整体升高消失，明确为C15参与的组合状态信号。
原始TYMP—比例rho0.468→0.291，U P0.274/FDR0.411；研究调整U rho0.107/P0.725。恶性糖酵解—比例rho0.565→0.682，rho差分层bootstrap95%CI0.0089–0.275，但研究调整仅名义显著，原50程序FDR未过。MTA响应RNA代理仍为负向，不能据新STAT3信号确认MTA串联。主轨迹覆盖442/2674(16.53%)→463/4502(10.28%)；无应激/增殖模型去掉C15，不能当完整U验证。原C4脂质/溶酶体保留但多基因特征稀释。
建议联合群用于SPP1总体来源/丰度与上游关系，C4/C15/C10分别解释脂质溶酶体、增殖/新STAT3糖酵解、患者集中候选；没有按显著性挑选最优组合。与C4定义无关的公共实测MTA、旧SPP1空间状态/临床及第一部分保留原定义，未偷换为U证据。
权威：{D}；本地：C:/Users/31937/Downloads/26210/30_C4_C10_C15联合复评_v1。摘要/完整报告+统计表，21/21文件哈希一致。项目AGENTS及part2 STATUS记录新增联合分析授权与完成状态，旧C4定义不变。
当前断点：数值复评已交付，尚未据联合视角重画Fig3或改Phase191正文。Fig1/2仍暂定稿，Fig3 phase199为原C4视角版本；所有本轮计算结束，无活跃任务。
"""
for n in ['PROJECT_STATUS.md','DECISIONS.md','CODEX_HANDOFF.md']:
 with (R/n).open('a',encoding='utf8') as f:f.write(note)
with (R/'RESULTS_MANIFEST.tsv').open('a',encoding='utf8') as f:
 f.write(f'\nphase200_union_reassessment\tphase200\tnumerical_analysis_and_reports\t{D}/delivery\t{now}\tPASS\t2301 independent checks; attribution checks; 21 local hash matches\tFixed C4+C10+C15; original C4 retained; no new figures\n')
p=D/'STATUS.json';z=json.loads(p.read_text());z.update(status='COMPLETE',completed=now,local_hash_matches=21,local_delivery='C:/Users/31937/Downloads/26210/30_C4_C10_C15联合复评_v1');p.write_text(json.dumps(z,ensure_ascii=False,indent=2))
(D/'qc/LOCAL_DELIVERY_VALIDATION.json').write_text(json.dumps({'files':21,'hash_mismatches':0,'checked_at':now},indent=2))
print('PHASE200_COMPLETE')
