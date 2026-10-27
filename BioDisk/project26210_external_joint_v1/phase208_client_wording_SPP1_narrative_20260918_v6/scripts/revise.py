from pathlib import Path
import json,re,copy,shutil
r=Path.cwd();d=r/'phase208_client_wording_SPP1_narrative_20260918_v6';old=r/'phase207_manuscript_figures_crossrefs_20260918_v5';b=json.loads((old/'report/report_materialized.json').read_text());orig=copy.deepcopy(b)
for n in range(1,5):shutil.copy2(old/'figures'/f'Figure{n}.pdf',d/'figures'/f'Figure{n}.pdf')
repl=[('直接 TYMP 关联对技术协变量调整的敏感程度不同','直接 TYMP 关联在技术协变量调整后的变化幅度不同'),('0.225 和 0.275 的阈值敏感性分析','0.225 和 0.275 两种距离阈值下的比较'),('逐患者的四蛋白敏感性结果','逐患者的四蛋白对应结果'),('技术敏感性分析','文库量及检出基因数调整分析'),('技术敏感性模型','文库量及检出基因数调整模型'),('广义敏感性代理','广义功能代理'),('去应激/增殖敏感性模型','排除应激和增殖细胞后的模型'),('去应激/周期敏感性模型','排除应激和周期细胞后的模型'),('对邻近距离阈值敏感','在邻近距离阈值下的克隆划分发生变化'),('靶向敏感性分析','候选调控子的补充评估'),('主分析和技术敏感性','未调整与文库量及检出基因数调整结果'),('敏感性阈值','其他距离阈值'),('阈值敏感性','不同距离阈值比较'),('主克隆敏感性','不同阈值下的主克隆划分'),('敏感性端点','不同模型设定下的端点'),('敏感性模型','排除应激和增殖细胞后的模型'),('虚线为敏感性24例','虚线为排除应激和增殖细胞后的模型24例'),('敏感性','去应激/增殖模型'),('18 端点审核家族','18 端点检验家族'),('谱系审核','谱系鉴定'),('技术审计和完整候选材料','分析过程记录和完整候选材料'),('使用足以复核的数值精度','保留完整数值精度'),('主要反向结果不因本次文字重组而删除','各端点的方向和统计结果分别呈现')]
def sub(v):
 if isinstance(v,str):
  for a,z in repl:v=v.replace(a,z)
  return v
 if isinstance(v,list):return [sub(x) for x in v]
 if isinstance(v,dict):return {k:sub(x) for k,x in v.items()}
 return v
b=sub(b)
for x in b:
 if x.get('text','').startswith('患者内 CNV 状态显示'):x['text']=x['text'].replace('5 位患者的克隆划分在邻近距离阈值下的克隆划分发生变化','5 位患者在邻近距离阈值下的克隆划分发生变化')
 if x['type']=='note' and x.get('text','').startswith('E：原五类动态表达曲线'):x['text']=x['text'].replace('实线为主分析14例，虚线为排除应激和增殖细胞后的模型24例','Model 1实线为主分析14例，Model 2虚线为排除应激和增殖细胞后的模型24例')
start=next(i for i,x in enumerate(b) if x['type']=='h2' and x['text'].startswith('2.9 '));end=next(i for i in range(start+1,len(b)) if b[i]['type']=='h2')
new=[dict(type='h2',text='2.9 巨噬细胞状态鉴定揭示多个SPP1富集群'),dict(type='p',text='上述结果将肿瘤细胞的核苷、蛋氨酸和多胺代谢特征连接到MTA这一候选代谢物。接下来，我们转向同一单细胞资源中的髓系细胞，解析与这些代谢线索相衔接的细胞状态。对髓系细胞进行无监督细分，并结合巨噬谱系标志完成鉴定后，在39名患者的22,444个成熟巨噬细胞中识别出12个候选状态（图3D）。'),dict(type='p',text='对这12个状态进行标志基因比较时，SPP1在C4、C10和C15中呈现较高表达，同时伴随不同的分子特征：C4共同富集GPNMB、CCL18和CYP27A1，C10呈USP53/ALOX5AP特征，C15呈MKI67/TOP2A增殖特征（图4A）。因此，SPP1由分群后的标志鉴定自然进入后续研究，它对应的是分布于多个巨噬状态的表达信号。三群在图谱中的位置及各自表达水平见图3D、图3E。'),dict(type='p',text='为描述这一信号的来源范围，将C4、C10和C15合计为三群联合，共4,502个细胞。在原16名同文库配对患者中，联合群占成熟巨噬细胞的比例中位数为9.76%，C4单群为7.41%；联合群贡献捕获SPP1原始RNA的份额中位数为59.63%，C4为54.59%（图3E）。三群联合由此补充了SPP1在不同巨噬状态中的来源信息，后续再根据各群的分子身份和患者支持范围确定重点研究对象。'),dict(type='p',text='C4出现在37名患者和五项研究，26名患者满足群内与其余成熟巨噬均至少20个细胞的配对条件；C10仅2名患者满足该条件，且97.64%的细胞来自同一患者；C15有9名有效配对患者，并具有明确的增殖特征（有效配对患者数见图3E）。结合患者覆盖、多基因身份及研究内方向一致性，后续以C4开展主要功能、调控和状态排序分析，同时保留C10和C15用于解释SPP1来源的异质性。')]
b[start:end]=new
# Keep 2.10 focused on the independent evidence for choosing C4, rather than repeating the discovery step.
for x in b:
 if x.get('text','').startswith('基于上述候选比较，C4作为后续主分析状态。'):
  x['text']='C4的多基因表达组合进一步明确了这一群体的分子身份。该群包含2,674个细胞，出现在37名患者及五项研究，GPNMB、CCL18、CYP27A1、LPL、APOE与SPP1共同富集（图4A）。在此基础上，我们进一步按患者比较C4与其他成熟巨噬细胞的表达差异。'
b[1]['text']='项目 26210　四图整合与SPP1发现叙事修订版　2026-09-18'
assert not re.search(r'敏感|sensitivity|稳健性|鲁棒性|审计|审核',json.dumps(b,ensure_ascii=False),re.I)
ot=[x for x in orig if x['type']=='table'];nt=[x for x in b if x['type']=='table'];assert len(ot)==len(nt)==32
for a,z in zip(ot,nt):
 assert len(a['rows'])==len(z['rows'])
 for ar,zr in zip(a['rows'],z['rows']):
  for av,zv in zip(ar,zr):assert re.findall(r'-?\d+(?:\.\d+)?(?:e[+-]?\d+)?',str(av),re.I)==re.findall(r'-?\d+(?:\.\d+)?(?:e[+-]?\d+)?',str(zv),re.I)
(d/'report/report_materialized.json').write_text(json.dumps(b,ensure_ascii=False,indent=2))
changes=[dict(before=a,after=z) for a,z in zip(orig[:start],b[:start]) if a!=z];(d/'qc/WORDING_CHANGES.json').write_text(json.dumps(changes,ensure_ascii=False,indent=2))
s=(old/'scripts/render_report.py').read_text().replace('phase207_manuscript_figures_crossrefs_20260918_v5','phase208_client_wording_SPP1_narrative_20260918_v6').replace('26210_TYMP_MTA_SPP1_four_figures_crossrefs_v5','26210_TYMP_MTA_SPP1_four_figures_narrative_v6');(d/'scripts/render_report.py').write_text(s)
print('TEXT_READY; numerical table content unchanged')
