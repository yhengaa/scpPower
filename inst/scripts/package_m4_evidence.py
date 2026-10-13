"""Package the checked M4 evidence without source caches or private data."""
from pathlib import Path
import zipfile,json,hashlib,time,tarfile
start=time.monotonic();root=Path('.');files=[]
for directory in ['R','man','tests','inst/scripts','inst/presets','reports/M4','.github']:
 for p in Path(directory).rglob('*'):
  if p.is_file() and not any(a in p.parts for a in ['__pycache__','points','cluster-logs','assignments','layouts']):files.append(p)
files+=list(map(Path,['DESCRIPTION','NAMESPACE','NEWS.md','README.md','renv.lock','scpPower_mandate_v1.1.md','.Rbuildignore','.Rprofile','renv/activate.R']))
for p in Path('reports/M3').glob('reading-grid*'):files.append(p)
files+=[Path('reports/M3/power-report.md'),Path('reports/M3/power-report-v0.2.0.md')]
m=json.loads(Path('reports/M4/run-manifest.json').read_text())
version=next(line.split(':',1)[1].strip() for line in Path('DESCRIPTION').read_text().splitlines() if line.startswith('Version:'))
dest=Path('scpPower_M4_v'+version+'_evidence.zip')
with zipfile.ZipFile(dest,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in sorted(set(files)):z.write(p,str(p))
 with tarfile.open('reports/M4/first-assignments-and-layouts.tar.gz','r:gz') as tar:
  for item in tar:z.writestr('reports/M4/'+item.name,tar.extractfile(item).read())
with zipfile.ZipFile(dest) as z:assert z.testzip() is None
receipt={'seed':m['seed'],'versions':m['versions'],'path':str(dest),'sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),'bytes':dest.stat().st_size,'input_hashes':{'M4_manifest':hashlib.sha256(Path('reports/M4/run-manifest.json').read_bytes()).hexdigest()},'wall_seconds':time.monotonic()-start}
Path('scpPower_M4_v'+version+'_evidence-manifest.json').write_text(json.dumps(receipt,indent=2)+'\n')
print('M4 evidence packaged:',dest,receipt['bytes'],'bytes')
