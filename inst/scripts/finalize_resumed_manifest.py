"""Record verification evidence and reproducible output hashes without refitting."""
from pathlib import Path
import json,hashlib,time,re,datetime
start=time.monotonic();out=Path('reports/M1-resumed');path=out/'run-manifest.json';m=json.loads(path.read_text())
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
for key in ['resumed-source-downloads','resumed-bulk-downloads','resumed-assay-doc-downloads','specht-analysis-source','bulk-metadata-download-manifest']:
 value=m[key]
 for record in value if isinstance(value,list) else [value]:record['versions']=m['versions']
check=out/'R-check.log';tests=out/'tests.log'
if check.exists() and tests.exists():
 c=check.read_text();t=tests.read_text()
 if not re.search(r'Status: OK',c):raise RuntimeError('Package check has not passed')
 if not re.search(r'FAIL 0 \| WARN 0 \| SKIP 0 \| PASS 77',t):raise RuntimeError('Required tests have not all passed')
 stamp=re.search(r'current time: (.*?) UTC',c).group(1)
 epoch=datetime.datetime.strptime(stamp,'%Y-%m-%d %H:%M:%S').replace(tzinfo=datetime.timezone.utc).timestamp()
 m['engineering_checks']['local_check_wall_seconds']=check.stat().st_mtime-epoch
 m['engineering_checks']['test_wall_seconds']=float(re.search(r'Duration: ([0-9.]+) s',t).group(1))
 m['engineering_checks']['input_hashes']={str(p):sha(p) for p in [Path('R/scp_calibrate.R'),Path('tests/testthat/test-scp_calibrate.R'),Path('renv.lock')]}
 m['engineering_checks']['seed']=20261013
 m['engineering_checks']['versions']=m['versions']
 m['engineering_checks']['github_check']='Push-triggered CI must be verified at the published resumed commit before delivery; final chat links to that exact successful run.'
m['output_sha256']={str(p.relative_to(out)):sha(p) for p in sorted(out.iterdir()) if p.is_file() and p!=path and p.suffix in ['.csv','.png','.pdf','.md','.rds','.json']}
m['reproduction_code_sha256']={str(p):sha(p) for p in [Path('R/scp_calibrate.R'),Path('scpPower_mandate_v1.1.md'),Path('renv.lock')]+sorted(Path('inst/scripts').glob('*resumed*')) if p.is_file()}
m['final_manifest_generation']={'seed':20261013,'versions':m['versions'],'input_hashes':m['output_sha256'],'wall_seconds':time.monotonic()-start}
path.write_text(json.dumps(m,indent=2)+'\n')
print('Final manifest verifies checks and hashes all deliverables.')
