from pathlib import Path
import ast,subprocess,time
ssh=['ssh','-T','-F','C:/Users/31937/.ssh/config','-o','BatchMode=yes','-o','ClearAllForwardings=yes','-o','UserKnownHostsFile=C:/Users/31937/.ssh/known_hosts','-i','C:/Users/31937/.ssh/id_ed25519','sgs10-codex']
def call(code):
 for i in range(4):
  p=subprocess.run(ssh+['python3 -'],input=code.encode(),capture_output=True)
  if p.returncode!=255:break
  time.sleep(3)
 print(p.stdout.decode(errors='replace'),p.stderr.decode(errors='replace'),flush=True);p.check_returncode();return p.stdout
dest=call("from pathlib import Path\np=sorted(p for p in Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1').glob('local_records_migration_*') if (p/'MIGRATION_MANIFEST.json').exists())[-1];print(p)\n").decode().strip()
tree=ast.parse(Path('C:/codex_test/migrate_26210.py').read_text(encoding='utf8'))
node=next(n for n in tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='verify' for t in n.targets))
verify=eval(compile(ast.Expression(node.value),'<verification>','eval'),{'dest':dest})
verify=verify.replace('DEST',repr(dest))
time.sleep(3);call(verify)
