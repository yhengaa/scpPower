"""Hash a completed stage's public evidence and local validation logs."""
from pathlib import Path
import json,hashlib,sys,re,time,datetime
start=time.monotonic();out=Path('reports')/sys.argv[1];path=out/'run-manifest.json';m=json.loads(path.read_text())
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
check=out/'R-check.log';tests=out/'tests.log'
if check.exists():
 if 'Status: OK' not in check.read_text():raise RuntimeError('Local R CMD check is not clean')
 t=tests.read_text();match=re.search(r'FAIL 0 \| WARN 0 \| SKIP 0 \| PASS ([0-9]+)',t)
 if not match:raise RuntimeError('Required tests did not pass without warnings/skips')
 m['engineering_checks']=dict(local_R_CMD_check='Status: OK',assertions_passed=int(match.group(1)),tests_wall_seconds=float(re.search(r'Duration: ([0-9.]+) s',t).group(1)),seed=20261013,versions=m['versions'],input_hashes={'package_code':sha(Path('R/scp_calibrate.R')),'test_code':sha(Path('tests/testthat/test-scp_calibrate.R')),'renv_lock':sha(Path('renv.lock'))})
 m['engineering_checks']['input_hashes'].update({str(p):sha(p) for p in list(Path('R').glob('*.R'))+list(Path('tests/testthat').glob('*.R'))})
 text='\n'.join(p.read_text() for p in Path('tests/testthat').glob('*.R'))
 m['engineering_checks']['test_seeds']=sorted({int(a or b) for a,b in re.findall(r'set.seed\((\d+)\)|seed\s*=\s*(\d+)',text)})
 stamp=re.search(r'current time: (.*?) UTC',check.read_text()).group(1)
 epoch=datetime.datetime.strptime(stamp,'%Y-%m-%d %H:%M:%S').replace(tzinfo=datetime.timezone.utc).timestamp()
 m['engineering_checks']['local_check_wall_seconds']=check.stat().st_mtime-epoch
 m['engineering_checks']['remote_CRAN_lookup']='Disabled for the local check after a public CRAN archive timeout; hosted CI performs the full check.'
m['output_sha256']={str(p.relative_to(out)):sha(p) for p in sorted(out.iterdir()) if p.is_file() and p!=path and p.suffix in ['.csv','.json','.png','.pdf','.md']}
m['reproduction_code_sha256']={str(p):sha(p) for p in list(Path('R').glob('*.R'))+list(Path('inst/scripts').glob('*'))+[Path('scpPower_mandate_v1.1.md'),Path('renv.lock')] if p.is_file()}
m['final_manifest_generation']=dict(seed=20261013,versions=m['versions'],wall_seconds=time.monotonic()-start,input_hashes=m['output_sha256'])
path.write_text(json.dumps(m,indent=2)+'\n');print('Stage manifest finalized:',out)
