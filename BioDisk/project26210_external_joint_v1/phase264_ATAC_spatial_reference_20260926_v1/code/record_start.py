from pathlib import Path
from datetime import datetime,timezone,timedelta
import json,shutil
R=Path.cwd();P=R/'phase264_ATAC_spatial_reference_20260926_v1';Q=R/'phase265_spatial_reference_preview_20260926_v1'
now=datetime.now(timezone(timedelta(hours=8))).isoformat()
note=f"""## {now} Phase264 ATAC足迹空间投射试点进行中 Phase265 RNA对照预览完成
用户要求按《分析待做（新）》4.1参考图先出图，并追加尝试参考图原方法的ATAC结合活性到空间投射。已识别原图Kuppe et al Nature 2022 Fig1h-j（doi10.1038/s41586-022-05060-x）：cell2location、PROGENy及HINT-ATAC/JASPAR2020结合活性按细胞类型比例加权。
现有pySCENIC确实完成，但属于RNA网络，不是ATAC结合矩阵。Phase265已完成A4 24图RNA对照：两病例×4细胞丰度、4既有ssGSEA代谢通路、4空间AUCell（HNF4A/ILF2/JUNB/CEBPD）。HOXD3目标覆盖43/54未达事先80%门槛，资格记录保留；无统计筛选或新推断P。全部7097spot坐标匹配，AUCell独立公式最大误差2.78e-16；PDF实际视觉检查PASS，PDF/PNG/SVG/中文图注4文件本地SHA256全部一致。目录phase265_spatial_reference_preview_20260926_v1；本地Downloads/26210/Phase265_空间参考图_RNA对照_v1。当前完整报告Phase263未修改。
Phase264正从GSE207493公开ATAC包顺序提取首个样本RCC100（物理归档顺序选择，不按信号），仅取该样本，zip逐成员CRC验证。原包为gzip/tar内嵌38.25GB ZIP；首项fragments.tsv.gz为2.155GB。直连约20KB/s，保留direct_attempt.partial；已转专用代理127.0.0.1:17897 -> Windows7897，本地隧道PID30472。下载PID见download.pid，代码extract_RCC100_proxy.py，完成须核对DOWNLOAD_EXIT_STATUS和DOWNLOAD_STATUS及成员CRC。
同时在本项目software/venv安装RGT1.0.2（prepare_env.py，不能改共享环境），下载官方GENCODE GRCh38 genome及v38 basic GTF（download_reference.py）。尚未生成ATAC TF矩阵、ATAC空间图，不能将RNA预览当ATAC完成。
恢复入口：查看本阶段各EXIT_STATUS、下载清单及logs；完成后核对ATAC输入、细胞类型对应与HINT输入资格再计算。若只有fragment端点，必须明确任何兼容表示及坐标/片段长度无损校验，不能伪造测序读段或把chromVAR分数称HINT结合活性。没有未知细胞类型的零值填充或假定新C4映射。
统计遵循用户2026-09-26新指令：有效检验先按原始P<0.05筛名义显著，FDR完整保留且不作唯一门槛；本次单ATAC样本及两空间病例仅描述，不补造P。Windows新规则tty=true/login=false，辅助进程隐藏；不修改全局终端配置。全部进度仅服务器。

"""
backup=P/'qc/root_records_before_start';backup.mkdir(parents=True,exist_ok=True)
for fn in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
 f=R/fn
 if not (backup/fn).exists():shutil.copy2(f,backup/fn)
 f.write_text(note+f.read_text())
(P/'STATUS.json').write_text(json.dumps(dict(status='ATAC_INPUT_PREPARATION_RUNNING',time=now,RNA_preview='COMPLETE_VALIDATED_DELIVERED',atac_activity_matrix=False,atac_projection=False),indent=2))
(Q/'STATUS.json').write_text(json.dumps(dict(status='COMPLETE_VALIDATED_DELIVERED',time=now,maps=24,atac_included=False),indent=2))
print('PROGRESS_RECORDED',flush=True)

