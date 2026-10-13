"""Package J, M2b and M4b, preserving the accepted earlier validation."""
from pathlib import Path
import zipfile,tarfile,json,hashlib,time
start=time.monotonic();files=[]
for directory in ['R','man','tests','inst/scripts','inst/presets','reports/J','reports/M2b','reports/M4b','.github']:
 for p in Path(directory).rglob('*'):
  if p.is_file() and not any(a in p.parts for a in ['__pycache__','points','cluster-logs','assignments','layouts']):files.append(p)
files += [Path(p) for p in ['DESCRIPTION','NAMESPACE','NEWS.md','README.md','renv.lock','renv/activate.R','.Rprofile','.Rbuildignore','scpPower_mandate_v1.1.md','reports/M4/validation-report.md','reports/M4/validation-report-v0.3.0.md','reports/M4/validation-table.csv','reports/M4/dataset-verdicts.csv','reports/M4/run-manifest.json']]
version=next(l.split(':',1)[1].strip() for l in Path('DESCRIPTION').read_text().splitlines() if l.startswith('Version:'))
dest=Path('scpPower_TMT_round_v'+version+'_evidence.zip')
with zipfile.ZipFile(dest,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in sorted(set(files)):z.write(p,str(p))
 archive=Path('reports/M4b/first-assignments-and-layouts.tar.gz')
 if archive.exists():
  with tarfile.open(archive,'r:gz') as tar:
   for item in tar:z.writestr('reports/M4b/'+item.name,tar.extractfile(item).read())
with zipfile.ZipFile(dest) as z:assert z.testzip() is None
m=json.loads(Path('reports/J/run-manifest.json').read_text())
r={'seed':m['seed'],'versions':m['versions'],'archive':str(dest),'bytes':dest.stat().st_size,'sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),'input_hashes':{'diagnostic_J':hashlib.sha256(Path('reports/J/run-manifest.json').read_bytes()).hexdigest(),'M2b':hashlib.sha256(Path('reports/M2b/run-manifest.json').read_bytes()).hexdigest()},'wall_seconds':time.monotonic()-start}
Path('scpPower_TMT_round_v'+version+'_evidence-manifest.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps({k:r[k] for k in ['archive','bytes','sha256']}))
