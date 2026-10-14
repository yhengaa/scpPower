"""Check release scope and unmodified registered numerical code/presets."""
from pathlib import Path
import csv,json,hashlib,subprocess,time,re
start=time.monotonic();sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
base='36f39044f6d5e128283dd4470b0ff1e82d4e38b7'
for p in ['R/scp_calibrate.R','R/scp_simulate.R','R/scp_power.R','renv.lock']+[str(p) for p in Path('inst/presets').glob('*.yaml')]:
 old=subprocess.check_output(['git','show',base+':'+p])
 assert hashlib.sha256(old).hexdigest()==sha(p),p
assert 'Version: 1.0.0' in Path('DESCRIPTION').read_text().splitlines()
news=Path('NEWS.md').read_text().split('# scpPower 0.3.1')[0]
assert len([line for line in news.splitlines() if line.startswith('* ')])==10
vignette=Path('vignettes/designing-scp.Rmd').read_text()
assert len(vignette.splitlines())<150
assert ' --no-build-vignettes' not in Path('.github/workflows/R-CMD-check.yaml').read_text()
assert 'VignetteBuilder: knitr' in Path('DESCRIPTION').read_text()
statement="Code was written with the assistance of GPT (OpenAI) and independently reviewed with Claude (Anthropic), both under the author's direction. All scientific decisions, registered criteria, data choices and verifications are the author's."
assert statement in Path('README.md').read_text()
blocks=re.findall(r'```r\n(.*?)```',Path('README.md').read_text(),re.S)
for block in blocks[1:4]:assert len(block.strip().splitlines())<=10
reg=list(csv.DictReader(open('reports/M3b/registered-grid.csv')))
original=list(csv.DictReader(open('reports/M3/registered-grid.csv')))
assert reg==original[:630]
assert len(list(csv.DictReader(open('reports/M3b/power-grid.csv'))))==945
old_tests=subprocess.check_output(['git','ls-tree','-r','--name-only',base,'tests'],text=True).splitlines()
for p in old_tests:assert subprocess.check_output(['git','show',base+':'+p])==Path(p).read_bytes(),p
m=json.loads(Path('reports/M3b/run-manifest.json').read_text())
out={'seed':20261013,'versions':m['versions'],'base_commit':base,'registered_R_code_and_presets_unchanged':True,'all_existing_tests_unchanged':True,'new_tests':'tests/testthat/test-pilot-conversion.R','TMT_grid_and_seeds_identical_to_original':True,'README_example_code_under_ten_lines':True,'vignette_source_lines':len(vignette.splitlines()),'input_hashes':{p:sha(p) for p in ['README.md','DESCRIPTION','NEWS.md','CITATION.cff','vignettes/designing-scp.Rmd','inst/doc/pilot-protocol.md']},'wall_seconds':time.monotonic()-start}
Path('reports/M5/release-audit.json').write_text(json.dumps(out,indent=2)+'\n')
print('Release scope and reproducibility audit passed.')
