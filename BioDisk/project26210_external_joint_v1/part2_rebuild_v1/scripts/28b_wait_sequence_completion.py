from pathlib import Path
import time,subprocess
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/part2_rebuild_v1');out=R/'logs/phase148_sequence_v2_exit.txt'
# Completion monitoring only: no scientific result inspection, once per15minutes.
while not out.exists():
 ps=subprocess.check_output(['ps','-eo','args'],text=True)
 alive=any(('17_export_trajectory_matrices.py' in s or '18b_run_partitioned_monocle3_header_fixed.R' in s or '27b_remaining_sequence_v2.sh' in s) and not s.startswith('grep ') for s in ps.splitlines())
 if not alive:
  print('PIPELINE_PROCESS_EXITED_WITHOUT_COMPLETION_MARKER',flush=True);break
 time.sleep(900)
if out.exists():print(out.read_text(),flush=True)
