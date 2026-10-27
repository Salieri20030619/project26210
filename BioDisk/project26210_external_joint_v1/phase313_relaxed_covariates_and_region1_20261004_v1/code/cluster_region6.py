import sys,runpy
from pathlib import Path
sys.argv=[str(Path(__file__).with_name('cluster_region.py')),'region6']
runpy.run_path(sys.argv[0],run_name='__main__')
