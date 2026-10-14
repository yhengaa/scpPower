"""Build the numerical vignette, run R CMD check, and retain run provenance."""
from pathlib import Path
import subprocess,time,json,hashlib,os
start=time.monotonic();out=Path('reports/M5');sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
env=os.environ.copy();env['RENV_CONFIG_AUTOLOADER_ENABLED']='FALSE';env['R_LIBS_USER']=str(Path('.local/R-library').resolve());env['_R_CHECK_CRAN_INCOMING_REMOTE_']='false'
for args,name in [(['R','CMD','build','.','--no-manual'],'R-build.log'),(['R','CMD','check','scpPower_1.0.0.tar.gz','--no-manual','--as-cran'],'R-check.log')]:
 with (out/name).open('w') as log:
  x=subprocess.run(args,stdout=log,stderr=subprocess.STDOUT,env=env)
 if x.returncode:raise RuntimeError(f'{name} failed; preserve the log')
log=Path('scpPower.Rcheck/00check.log').read_text()
import re
if re.search(r'[1-9][0-9]* (ERROR|WARNING)',log):raise RuntimeError('Check has errors or warnings')
context=json.loads((out/'worked-examples-manifest.json').read_text())
m={'seed':20261013,'versions':context['versions'],'input_hashes':{p:sha(p) for p in ['DESCRIPTION','renv.lock','vignettes/designing-scp.Rmd','R/scp_calibrate.R','R/scp_simulate.R','R/scp_power.R']},'output_hashes':{p.name:sha(p) for p in [out/'R-build.log',out/'R-check.log',Path('scpPower.Rcheck/00check.log')]},'vignette_built_and_checked':True,'wall_seconds':time.monotonic()-start,'remote_CRAN_incoming_disabled_locally_only':True}
(out/'local-check-manifest.json').write_text(json.dumps(m,indent=2)+'\n')
print('Release build and check passed, including the numerical vignette.')
