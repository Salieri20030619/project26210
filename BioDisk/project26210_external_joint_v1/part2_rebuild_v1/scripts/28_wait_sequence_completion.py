from pathlib import Path
import time,subprocess
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/part2_rebuild_v1');out=R/'logs/phase148_sequence_exit.txt'
# Completion monitoring only: no scientific result inspection, once per15minutes.
while not out.exists():
 ps=subprocess.check_output(['ps','-eo','args'],text=True)
 alive=any(('17_export_trajectory_matrices.py' in s or '18_run_partitioned_monocle3.R' in s or '27_remaining_sequence_after_trajectory.sh' in s) and not s.startswith('grep ') for s in ps.splitlines())
 if not alive:
  print('PIPELINE_PROCESS_EXITED_WITHOUT_COMPLETION_MARKER',flush=True);break
 time.sleep(900)
if out.exists():print(out.read_text(),flush=True)
