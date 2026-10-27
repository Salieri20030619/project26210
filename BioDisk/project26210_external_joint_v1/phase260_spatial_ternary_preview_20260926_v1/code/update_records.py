from pathlib import Path
import json,datetime,shutil,csv,fcntl,os
p=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase260_spatial_ternary_preview_20260926_v1');root=p.parent
qc=json.loads((p/'qc/FINAL_VALIDATION.json').read_text());assert qc['status']=='PASS'
assert (p/'build_EXIT_STATUS').read_text().strip()=='0' and (p/'finalize_EXIT_STATUS').read_text().strip()=='0'
now=datetime.datetime.now().astimezone().isoformat()
entry=f"""## {now} Phase260 参考Fig5D的独立三元空间距离预览完成
用户要求先生成单图看效果，本轮未加入空间大图。独立图{p}/delivery/26210_Ternary_spatial_preview_v1.pdf，160×170mm，上方GLMF1/2叠加、下方分别显示；蓝/珊瑚色区分病例。
参照为冻结恶性比例>=0.20的spot，三目标为SPP1-state q05/内皮q05/T-cell q05各自片内前20%区域。沿实际100um左右边长的冻结组织图计算最短图距离，不跨分量/空洞。三元坐标为d_i/三距之和；靠近顶点代表该项距离占比大即相对更远，不是细胞比例或更邻近。真实点全部显示，不抖动/平滑/挑点，只有2患者描述，无推断检验。
GLMF1参照1314、绘图1282，不可达11/三距皆零21；GLMF2参照3133、绘图3113，不可达0/三距皆零20。全部排除和原绝对距离逐spot保存，不加伪距离。两例分别对应Patient1/Patient7；不套原文TIB分组，T细胞不改称CD8，旧SPP1-state不改称C4。
独立heapq与scipy六组最短路全位置一致，三元反算、数量/输入哈希和PDF实际渲染验收通过。FINAL_VALIDATION PASS，构建/收尾退出0。PDF/SVG/300dpi PNG/600dpi TIFF/图注/绘图数据六文件已下载C:/Users/31937/Downloads/26210/Phase260_三元空间单图_v1/，SHA256全部匹配。
当前空间大图仍Phase259 v4；三元图仅预览待审。无新模型、旧RDS与报告不变，无运行本轮任务。全部进度仅服务器。
"""
lock=(root/'phase260_record_update.lock').open('a+');fcntl.flock(lock,fcntl.LOCK_EX)
backup=p/'qc/record_backup';backup.mkdir(exist_ok=True)
for name in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md']:
 f=root/name;old=f.read_text();assert 'Phase260 参考Fig5D的独立三元空间距离预览完成' not in old
 shutil.copyfile(f,backup/name);temp=f.with_name(name+'.phase260.tmp');temp.write_text(entry+'\n'+old,encoding='utf-8');os.replace(temp,f)
f=root/'RESULTS_MANIFEST.tsv';shutil.copyfile(f,backup/f.name)
with f.open('a',newline='') as h:
 w=csv.writer(h,delimiter='\t',lineterminator='\n')
 for i,m in enumerate(json.loads((p/'DELIVERY_MANIFEST.json').read_text()),1):
  w.writerow([f'phase260_ternary_preview_{i}',p.name,'standalone_figure_preview',m['path'],now,'PASS_DESCRIPTIVE',str(p/'qc/FINAL_VALIDATION.json'),'Two patients; raw graph-distance shares; not added to main figure; SHA256='+m['sha256']])
fcntl.flock(lock,fcntl.LOCK_UN)
(p/'DELIVERY_CONFIRMATION.json').write_text(json.dumps({'status':'PASS','time':now,'local_root':'C:/Users/31937/Downloads/26210/Phase260_三元空间单图_v1','files_verified':6,'verification':'Windows SHA256 matched server manifest','added_to_main_figure':False},ensure_ascii=False,indent=2))
(p/'qc/RECORD_UPDATE_VALIDATION.json').write_text(json.dumps({'status':'PASS','time':now,'updated':['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md','RESULTS_MANIFEST.tsv'],'entries':6},indent=2))
print('Standalone ternary preview delivered and registered.')
