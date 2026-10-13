"""Preserve every registered first assignment and native count matrix in one archive."""
from pathlib import Path
import gzip,tarfile,hashlib,json,time,io
start=time.monotonic();root=Path('reports/M4b');dest=root/'first-assignments-and-layouts.tar.gz'
files=sorted([*root.joinpath('assignments').glob('*.csv'),*root.joinpath('layouts').glob('*.csv')]);records={}
with dest.open('wb') as f,gzip.GzipFile(filename='',fileobj=f,mode='wb',mtime=0) as gz,tarfile.open(fileobj=gz,mode='w') as tar:
 for p in files:
  content=p.read_bytes();name=p.relative_to(root).as_posix();info=tarfile.TarInfo(name);info.size=len(content);info.mode=0o644;info.mtime=0
  tar.addfile(info,io.BytesIO(content));records[name]={'sha256':hashlib.sha256(content).hexdigest(),'bytes':len(content)}
with tarfile.open(dest,'r:gz') as tar:
 for item in tar:
  assert hashlib.sha256(tar.extractfile(item).read()).hexdigest()==records[item.name]['sha256']
m=json.loads(root.joinpath('run-manifest.json').read_text())
root.joinpath('layout-archive-manifest.json').write_text(json.dumps({'seed':m['seed'],'versions':m['versions'],'archive':dest.name,'sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),'members':records,'wall_seconds':time.monotonic()-start},indent=2)+'\n')
print('Consolidated and verified',len(files),'assignment and count-matrix files')
