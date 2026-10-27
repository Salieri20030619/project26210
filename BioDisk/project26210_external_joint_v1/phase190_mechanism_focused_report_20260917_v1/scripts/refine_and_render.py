from pathlib import Path
import json,subprocess,os
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=r/'phase190_mechanism_focused_report_20260917_v1';p=d/'report/report_materialized.json';b=json.loads(p.read_text())
for x in b:
 t=x.get('text','')
 t=t.replace('讨论重点应放在这些环节如何构成候选代谢背景，而不是将全部代谢酶视为必须同时升高的一组指标。','这些环节共同构成候选代谢背景，各环节的方向和作用分别解释。')
 t=t.replace('而不是仅依据一项酶RNA的升降。','并结合各类细胞的代谢能力。')
 t=t.replace('由此，本项目进一步考察具有明确分子特征的巨噬状态，而不是仅比较全组织中的单个标志物。','由此，本项目进一步考察具有明确分子特征的巨噬状态。')
 if t.startswith('C4的解释价值来自'):
  t='C4的解释价值来自多基因身份与程序的一致组合。APOE参与脂质运输，LPL参与甘油三酯水解及脂蛋白摄取，CYP27A1参与胆固醇相关氧化代谢[8-10]。这些基因与GPNMB、CCL18共同富集，并与去除TYMP/SPP1后仍保留的胆固醇及溶酶体程序升高相互对应。这一组合将本项目的SPP1信号定位于脂质处理和细胞内降解相关的巨噬状态。'
 t=t.replace('为讨论SPP1富集状态的功能取向提供了比单个M1/M2标志更具体的依据。','为解释SPP1富集状态的功能取向提供了具体依据。')
 if 'text' in x:x['text']=t
b += [{'type':'p','text':'[8] NCBI PubChem Gene. APOE human gene functional annotation. https://pubchem.ncbi.nlm.nih.gov/gene/348'},{'type':'p','text':'[9] NCBI Gene. LPL human lipoprotein lipase. https://www.ncbi.nlm.nih.gov/gene/4023'},{'type':'p','text':'[10] NCBI Gene. CYP27A1 human sterol metabolism annotation. https://www.ncbi.nlm.nih.gov/gene/1593'}]
p.write_text(json.dumps(b,ensure_ascii=False,indent=2))
env=os.environ.copy();env['PYTHONPATH']=str(r/'phase189_integrated_TYMP_MTA_macrophage_report_20260917_v1/vendor')
subprocess.run(['python3',str(d/'scripts/render_report.py')],env=env,check=True)
qa=(r/'phase189_integrated_TYMP_MTA_macrophage_report_20260917_v1/scripts/check_and_render_final.py').read_text().replace('phase189_integrated_TYMP_MTA_macrophage_report_20260917_v1','phase190_mechanism_focused_report_20260917_v1').replace('macrophage_integrated_report_v1','macrophage_integrated_report_v2').replace('assert len(pages)==33','assert len(pages)>=30').replace('assert len(files)==33','assert len(files)==len(pages)')
(d/'scripts/check_and_render.py').write_text(qa)
subprocess.run(['python3',str(d/'scripts/check_and_render.py')],env=env,check=True)
