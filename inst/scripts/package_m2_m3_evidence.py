"""Package the published M2 revision and M3 evidence without private data."""
from pathlib import Path
import zipfile,json,hashlib,time
start=time.monotonic();dest=Path('scpPower_M2_M3_v0.2.0_evidence.zip')
root=Path('.');files=[]
for d in ['R','man','tests','inst/scripts','inst/presets','reports/M2-revised','reports/M3','.github']:
 for p in Path(d).rglob('*'):
  if p.is_file() and not any(a in p.parts for a in ['__pycache__','points','cluster-logs']):files.append(p)
files += list(map(Path,['DESCRIPTION','NAMESPACE','NEWS.md','README.md','renv.lock','scpPower_mandate_v1.1.md']))
with zipfile.ZipFile(dest,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in sorted(set(files)):z.write(p,str(p))
m=json.loads(Path('reports/M3/run-manifest.json').read_text())
receipt=dict(seed=20261013,versions=m['versions'],path=str(dest),sha256=hashlib.sha256(dest.read_bytes()).hexdigest(),bytes=dest.stat().st_size,input_hashes={'manifest':hashlib.sha256(Path('reports/M3/run-manifest.json').read_bytes()).hexdigest()},wall_seconds=time.monotonic()-start)
Path('scpPower_M2_M3_v0.2.0_evidence-manifest.json').write_text(json.dumps(receipt,indent=2)+'\n')
print('Evidence bundle:',dest,receipt['bytes'],'bytes')
