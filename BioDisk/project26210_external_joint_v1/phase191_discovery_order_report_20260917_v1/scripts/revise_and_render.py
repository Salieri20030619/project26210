from pathlib import Path
import json,copy,re,subprocess,os
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');old=r/'phase190_mechanism_focused_report_20260917_v1';d=r/'phase191_discovery_order_report_20260917_v1'
for sub in ['report','scripts','qc','figures','supplement']:(d/sub).mkdir(parents=True,exist_ok=True)
a=json.loads((old/'report/report_materialized.json').read_text());b=copy.deepcopy(a)
def ix(start):return next(i for i,x in enumerate(b) if x.get('type')=='h2' and x.get('text','').startswith(start))
i8,i9,i10,i11=[ix(s) for s in ['2.8 ','2.9 ','2.10 ','2.11 ']]
s8=b[i8:i9];link=b[i9:i10];discovery=b[i10:i11]
s8[0]['text']='2.8 蛋氨酸与多胺代谢特征提供候选代谢物线索'
s8[1]['text']='TYMP相关恶性上皮状态中，半胱氨酸/蛋氨酸和精氨酸/脯氨酸程序升高，并伴随ODC1、SRM等多胺相关酶的表达变化（图6、表5及附表S5）。这些结果共同涉及SAM供给与多胺合成环节，提示可沿其具体代谢反应寻找连接细胞状态与局部代谢环境的候选产物。'
s8.insert(2,{'type':'p','text':'在多胺合成中，SRM和SMS分别参与亚精胺及精胺生成，其氨丙基转移反应同时产生5′-甲硫腺苷（5′-methylthioadenosine，MTA）；MTA随后由MTAP参与的分解过程进入蛋氨酸补救循环[1-3]。这一生成与清除关系将前述蛋氨酸/多胺表达特征连接到具体代谢物。因此，我们进一步利用配对组织RNA与实测代谢组，检验MTA与合成及清除背景的对应。'})
discovery[0]['text']='2.9 无监督髓系鉴定发现SPP1富集的C4巨噬状态'
discovery[1]['text']='为解析肿瘤髓系细胞的组成及状态异质性，我们在同一单细胞资源中进行谱系纯化、无监督细分和完整marker鉴定。在保留的12个成熟巨噬候选中，C4包含2,674个细胞，出现在37名患者及五项研究。GPNMB、CCL18、CYP27A1、LPL、APOE与SPP1的共同富集构成该群的多基因身份；SPP1由分群后的标志鉴定识别。'
discovery.append({'type':'p','text':'这一独立状态发现将SPP1确定为本项目髓系细胞中值得进一步追踪的分子特征。基于已定位的SPP1富集巨噬状态，下一步将其与肿瘤蛋氨酸/多胺代谢所引出的候选MTA联系起来，检验组织层关联并核查直接髓系作用依据。'})
link[0]['text']='2.10 围绕已发现的SPP1富集状态考察MTA相关联系'
link[1]['text']='鉴于SPP1在本项目巨噬候选中呈集中表达，我们进一步将这一发现与前述MTA代谢线索对应。在255例配对组织RNA与代谢组数据中，实测MTA与组织SPP1表达在四队列均呈正向，未调整合并rho=0.13894，P=0.02926；原多候选检验家族FDR=0.09755，定位为名义正向关联。髓系背景调整及加入TYMP条件后的完整模型结果见附表S11。'
link[-1]['text']='由独立鉴定的SPP1富集巨噬状态出发，组织层MTA-SPP1名义关联与直接MTA髓系扰动资料提供了进一步研究两者联系的依据。随后结合该群的功能程序和状态排序，描述SPP1信号所处的细胞背景。'
b=b[:i8]+s8+discovery+link+b[i11:]
# Swap marker table 9 and link evidence table 8 everywhere in the main-text numbering.
for x in b:
 for k in ['text','caption']:
  if k in x:
   x[k]=re.sub(r'(?<!附)表\s*8(?!\d)','表 __SWAP__',x[k])
   x[k]=re.sub(r'(?<!附)表\s*9(?!\d)','表 8',x[k])
   x[k]=x[k].replace('表 __SWAP__','表 9')
b[1]['text']='项目 26210　主线发现顺序修订版　2026-09-17'
b[5]['text']='在肿瘤代谢程序及具体生化反应的基础上，本研究将MTA纳入候选代谢物考察。独立的无监督髓系细分识别出2,674个细胞的SPP1富集C4；其SPP1在24/26名有效患者升高，并伴随胆固醇处理和溶酶体相关程序增强。围绕这一细胞状态发现，公共配对组织分析进一步探索性提示MTA与组织SPP1的名义正向关联，直接MTA髓系扰动资料提供生物学对应。报告由肿瘤代谢特征展开至候选代谢物、巨噬状态发现及其联系，客户TYMP敲低代谢图置于实验对应段。'
for x in b:
 t=x.get('text','')
 if t.startswith('本报告以恶性上皮细胞TYMP相关状态为起点'):
  x['text']='研究以恶性上皮细胞TYMP相关状态为起点，考察CNV组织背景、共识转录程序、调控与蛋白活性及核苷和氨基酸代谢；依据具体代谢反应引出MTA候选。同时，对同一资源中的髓系细胞进行独立细分与多基因鉴定，在发现SPP1富集巨噬状态后，进一步考察MTA与SPP1的组织联系，并解析该状态的功能程序及排序特征。'
 elif t.startswith('基于上述结果，本研究提出的工作模型是'):
  x['text']='基于上述结果，本研究提出的工作模型是：TYMP相关恶性细胞状态伴随核苷及蛋氨酸/多胺代谢组织变化，MTA作为其中具有明确生成和清除途径的候选代谢物，参与肿瘤代谢环境与髓系状态联系的研究。单细胞分析独立发现SPP1富集C4后，公共MTA-SPP1名义关联与直接MTA髓系作用依据，进一步支持将这两类发现置于同一机制框架中。'
 elif t.startswith('本研究显示TYMP相关恶性上皮状态具有核苷'):
  x['text']='本研究显示TYMP相关恶性上皮状态具有核苷、蛋氨酸及部分多胺代谢特征，并据具体生化反应引出MTA候选。独立髓系鉴定发现SPP1富集、脂质处理和溶酶体相关的C4巨噬状态；围绕该发现，公共组织关联及直接髓系作用资料进一步提供MTA-SPP1联系的研究依据。该框架构成后续免疫功能与空间组织研究的共同基础。'
main=b[:next(i for i,x in enumerate(b) if x.get('type')=='h1' and x.get('text','').startswith('4 '))]
assert not any('第一部分' in x.get('text','') for x in main)
oldtables=[x for x in a if x['type']=='table'];newtables=[x for x in b if x['type']=='table']
assert len(oldtables)==len(newtables)==31
assert all(any(x['headers']==y['headers'] and x['rows']==y['rows'] for y in newtables) for x in oldtables)
(d/'report/report_materialized.json').write_text(json.dumps(b,ensure_ascii=False,indent=2))
(d/'qc/CONTENT_VALIDATION.json').write_text(json.dumps({'tables_preserved':31,'figures_preserved':sum(x['type']=='figure' for x in b),'first_part_wording_removed_from_main':True,'sequence':['metabolic_clue_MTA','independent_SPP1_discovery','MTA_SPP1_link'],'table_numbers_swapped':[8,9]},indent=2))
src=(old/'scripts/render_report.py').read_text().replace('phase190_mechanism_focused_report_20260917_v1','phase191_discovery_order_report_20260917_v1').replace('macrophage_integrated_report_v2','macrophage_integrated_report_v3')
(d/'scripts/render_report.py').write_text(src)
env=os.environ.copy();env['PYTHONPATH']=str(r/'phase189_integrated_TYMP_MTA_macrophage_report_20260917_v1/vendor')
subprocess.run(['python3',str(d/'scripts/render_report.py')],env=env,check=True)
qa=(old/'scripts/check_and_render.py').read_text().replace('phase190_mechanism_focused_report_20260917_v1','phase191_discovery_order_report_20260917_v1').replace('macrophage_integrated_report_v2','macrophage_integrated_report_v3')
(d/'scripts/check_and_render.py').write_text(qa)
subprocess.run(['python3',str(d/'scripts/check_and_render.py')],env=env,check=True)
