"""Unified public evidence, including accepted historical files and layouts."""
from pathlib import Path
import subprocess,zipfile,tarfile,json,hashlib,time
start=time.monotonic();root=Path('.');sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
# Git's tracked file list is the release boundary, excluding local caches/logs.
files=[Path(x) for x in subprocess.check_output(['git','ls-files'],text=True).splitlines()]
files=[p for p in files if p.is_file() and not p.name.endswith('.zip')]
# The release workflow creates its exact-commit CI receipt after checkout.
receipt=Path('reports/M5/CI-receipt.json')
if receipt.exists():files.append(receipt)
files=sorted(set(files));dest=Path('scpPower_unified_v1.0.0_evidence.zip')
source=json.loads(Path('reports/M3b/run-manifest.json').read_text())
entry_hashes={str(p):sha(p) for p in files}
with zipfile.ZipFile(dest,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in files:z.write(p,str(p))
 for prefix,path in [('reports/M4/',Path('reports/M4/first-assignments-and-layouts.tar.gz')),('reports/M4b/',Path('reports/M4b/first-assignments-and-layouts.tar.gz')),('',Path('reports/M5/accepted-history.tar.gz'))]:
  with tarfile.open(path,'r:gz') as t:
   for member in t:
    if not member.isfile():continue
    name=prefix+member.name;data=t.extractfile(member).read()
    if name not in entry_hashes:z.writestr(name,data);entry_hashes[name]=hashlib.sha256(data).hexdigest()
 z.writestr('BUNDLE_MANIFEST.json',json.dumps({'version':'1.0.0','seed':20261013,'versions':source['versions'],'file_sha256':entry_hashes,'scope':'Unified M0–M5 sources and evidence, including archived accepted stages'},indent=2)+'\n')
with zipfile.ZipFile(dest) as z:assert z.testzip() is None;count=len(z.namelist())
m={'archive':str(dest),'bytes':dest.stat().st_size,'files':count,'sha256':sha(dest),'seed':20261013,'versions':source['versions'],'input_hashes':{'M3b_manifest':sha(Path('reports/M3b/run-manifest.json')),'historical_archive':sha(Path('reports/M5/accepted-history.tar.gz'))},'wall_seconds':time.monotonic()-start}
Path('scpPower_unified_v1.0.0_evidence-manifest.json').write_text(json.dumps(m,indent=2)+'\n')
print(json.dumps(m))
