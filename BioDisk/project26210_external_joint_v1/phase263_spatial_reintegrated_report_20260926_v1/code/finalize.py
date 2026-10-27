from pathlib import Path
import json,hashlib,shutil,fcntl
from datetime import datetime,timezone,timedelta
P=Path(__file__).resolve().parents[1];R=P.parent
now=datetime.now(timezone(timedelta(hours=8))).isoformat()
m=json.loads((P/'DELIVERY_MANIFEST.json').read_text())
assert m['sha256']=='558342d7b8f523bdac5c5ea8c65b2dff3c206c88e4b160e74b205d5622d55c8e'
assert m['bytes']==27556278
local='C:/Users/31937/Downloads/26210/Phase263_空间回归正文_v1/26210_integrated_spatial_report_v1.pdf'
(P/'qc/LOCAL_DELIVERY_CONFIRMED.json').write_text(json.dumps(dict(path=local,sha256=m['sha256'],bytes=m['bytes'],status='MATCH',time=now),ensure_ascii=False,indent=2))
summary=f"""## {now} Phase263 空间方法结果讨论并入前文并继承原编号完成
用户最新明确要求“把之前单列的空间部分重新放回到前文中，序号继承之前的”。本条覆盖Phase252/253仅追加、整体重组待审的当前执行限制；以Phase208完整排版源文件和Phase262已审空间正文为输入，按章节整合，不采用Phase251全文重写内容。
新版共39页：方法1.12–1.18并入原1.11后，结果2.14–2.21并入原2.13及图4图注后，讨论3.8–3.11并入原3.7后；引用11–18接原1–10合并至4.12。既有2.3排除与其余章节/表号不变，图5保持既有A–I及原生A4矢量，位于第24页，图注第25页。原32张表、全部统计和研究结论保留。原文仅更新版次副标题及“四图”为“五图”的图件说明；空间开场段移至2.14前，不重写摘要及结论。
333个正文/标题/图注文本块完整性、全部空间段落、原表对象、原输入哈希及页码检查通过；63个PDF图号链接有效，其中15个指向图5。五图均原尺寸矢量嵌入。实际检查39页，修整文献网址断行后复查变化的1/37/38/39页，其余35页像素一致。BUILD_EXIT_STATUS=0、VALIDATION_EXIT_STATUS=0、qc/FINAL_VALIDATION.json=PASS。首次校验的“摘要”子串误匹配已改用精确PDF书签；未改任何数据或统计。
正式交付：{m['path']}
本地交付：{local}。大小{m['bytes']}字节，SHA256={m['sha256']}，本地和服务器一致。当前完整报告以Phase263为准；Phase262及其他所有历史版本保留。
本轮仅文稿结构及排版，新增统计/模型/图件=0；免疫临床和参考3G分析未启动。无本轮运行任务。进度仅服务器。

"""
(P/'SESSION_SUMMARY_CN.md').write_text(summary)
(P/'README_REPRODUCE.md').write_text("""# Phase263 reproduction
Use the frozen report/report_materialized.json and unchanged figures/Figure1.pdf through Figure5.pdf.
Runtime: /mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
PYTHONPATH: ../phase198_fig2_contrast_corr_revision_20260918_v5/deps:../phase206_all_mainfigures_A4_20260918_v10/deps (absolute resolved paths recommended).
Run code/render_report.py followed by code/validate_render.py. Inspect rendered pages before code/prepare_delivery.py.
code/assemble.py records the initial structural assembly; the final renderer includes reference formatting polish. Run final code/render_report.py directly against the frozen materialized JSON for exact reproduction.
Original input hashes, editorial changes, software versions, page differences and checks are in qc. No new statistical analyses.
""")
backup=P/'qc/record_backup';backup.mkdir(exist_ok=True)
names=['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md','RESULTS_MANIFEST.tsv','AGENTS.md']
with (R/'phase263_record_update.lock').open('w') as lock:
 fcntl.flock(lock,fcntl.LOCK_EX)
 for name in names:
  if not (backup/name).exists():shutil.copy2(R/name,backup/name)
 for name in names[:3]:
  f=R/name;t=f.read_text()
  if not t.startswith(summary):f.write_text(summary+t)
 agents=R/'AGENTS.md'
 note="""
## 2026-09-26 用户要求空间部分回归前文并继承原编号
当前用户明确要求把此前单列的空间部分重新放回前文。本条覆盖此前当前阶段“仅追加、整体重组待审”的结构限制。
Phase263为当前完整报告：方法1.12–1.18接1.11、结果2.14–2.21接2.13、讨论3.8–3.11接3.7；图1–5、原表号及统计保持，空间文献11–18接原文献1–10。以Phase208和Phase262内容整合，不采用Phase251全文重写；不据本次结构授权启动免疫临床新分析、3G或重算。39页PDF在phase263_spatial_reintegrated_report_20260926_v1/delivery。进度仅服务器，旧版本保留。
"""
 if '## 2026-09-26 用户要求空间部分回归前文并继承原编号' not in agents.read_text():agents.write_text(agents.read_text()+note)
 with (R/'RESULTS_MANIFEST.tsv').open('a') as f:
  for rid,kind,path in [('phase263_integrated_report','report',Path(m['path'])),('phase263_final_validation','validation',P/'qc/FINAL_VALIDATION.json'),('phase263_materialized_report','source',P/'report/report_materialized.json')]:
   f.write('\t'.join([rid,P.name,kind,str(path),now,'PASS',str(P/'qc/FINAL_VALIDATION.json'),'Spatial sections reintegrated; original numbering and statistics retained; 39 pages; native five figures'])+'\n')
(P/'STATUS.json').write_text(json.dumps(dict(status='COMPLETE_VALIDATED_DELIVERED',time=now,pages=39,sha256=m['sha256'],running_tasks=0),ensure_ascii=False,indent=2))
print('COMPLETE_VALIDATED_DELIVERED; records updated; local SHA256 MATCH',flush=True)

