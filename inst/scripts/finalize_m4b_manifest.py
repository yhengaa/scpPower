"""Audit the registered M4 run and complete its immutable evidence manifest."""
from pathlib import Path
import json,hashlib,time,csv,subprocess
start=time.monotonic();out=Path('reports/M4b')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
m=json.loads((out/'run-manifest.json').read_text());reg=json.loads((out/'registration.json').read_text());frozen=json.loads((out/'predictions-frozen.json').read_text())
for p,h in reg['input_hashes']['code'].items():assert sha(p)==h,('Changed registered code',p)
for d,h in reg['input_hashes']['datasets'].items():assert sha('inst/presets/'+d+'.yaml')==h['released_preset'],('Changed released preset',d)
rows=list(csv.DictReader((out/'registered-grid.csv').open()));counted=[r for r in rows if r['counted']=='TRUE']
assert len(counted)==141 and frozen['successful_predictions']==141
for g in counted:
 task=g['task'];prefix=out/'points'/f'{int(task):04d}'
 p=json.loads(Path(str(prefix)+'-prediction-manifest.json').read_text())
 o=json.loads(Path(str(prefix)+'-observation-manifest.json').read_text())
 assert p['replicates']==o['replicates']==200
 assert p['completed_UTC']<=frozen['frozen_UTC']<=o['completed_UTC'],task
 assert o['input_hashes']['predictions_frozen']==sha(out/'predictions-frozen.json'),task
 assert sha(Path(str(prefix)+'-predicted.csv'))==frozen['prediction_hashes'][task]['results'],task
 assert sha(Path(str(prefix)+'-prediction-manifest.json'))==frozen['prediction_hashes'][task]['manifest'],task
for f in ['power-grid.csv','realised-fdr.csv','instrument-hours.csv','cheapest-designs.csv','replicate-metrics.csv.gz','power-curves.png','testability-curves.png','power-tmt-sensitivity.png','testability-tmt-sensitivity.png']:
 p='reports/M3/'+f
 old=subprocess.check_output(['git','rev-parse','v0.2.0:'+p],text=True).strip()
 assert subprocess.check_output(['git','hash-object',p],text=True).strip()==old,p
for dataset in ['derks2022_timsTOF','derks2022_QExactive','brunner2022']:
 p='inst/presets/'+dataset+'.yaml'
 assert subprocess.check_output(['git','hash-object',p],text=True).strip()==subprocess.check_output(['git','rev-parse','v0.3.0:'+p],text=True).strip(),p
for p in ['reports/M4/validation-table.csv','reports/M4/dataset-verdicts.csv','reports/M4/replicate-metrics.csv.gz']:
 assert subprocess.check_output(['git','hash-object',p],text=True).strip()==subprocess.check_output(['git','rev-parse','v0.3.0:'+p],text=True).strip(),p
m['accepted_M4_numbers_and_non_TMT_presets_unchanged']=True
m['M2b_manifest_sha256']=sha('reports/M2b/run-manifest.json')
m['diagnostic_J_manifest_sha256']=sha('reports/J/run-manifest.json')
m['registration']=reg
m['predictions_frozen']=frozen
m['package_version']=next(x.split(':',1)[1].strip() for x in Path('DESCRIPTION').read_text().splitlines() if x.startswith('Version:'))
m['cluster_accounting']={'path':'cluster-accounting.txt','sha256':sha(out/'cluster-accounting.txt')}
accounting=list(csv.DictReader((out/'cluster-accounting.txt').open(),delimiter='|'))
assert any(r['JobID']==m['cluster_execution']['job'] and r['State']=='COMPLETED' for r in accounting),'Cluster allocation not completed'
m['wall_seconds']=m['cluster_execution']['wall_seconds']+reg['wall_seconds']+m['report_wall_seconds']
m['wall_seconds_definition']='Preparation plus live compute execution plus report; per-point and per-replicate timings are retained separately; queue and allocation elapsed are in Slurm accounting.'
m['half_presets']={p.stem:{'path':str(p),'sha256':sha(p),'calibration_manifest':json.loads((out/(p.stem+'-half-preset-manifest.json')).read_text())} for p in Path('inst/presets/m4b').glob('*.yaml')}
m['diagnostic_J_and_closure']={'manifest':'reports/J/report-manifest.json','sha256':sha('reports/J/report-manifest.json'),'diagnostic_report_sha256':sha('reports/J/diagnostic-report.md'),'original_report_sha256':sha('reports/M4/validation-report-v0.3.0.md')}
if (out/'R-check.log').exists():
 log=(out/'R-check.log').read_text();tests=(out/'tests.log').read_text()
 assert 'Status: OK' in log
 import re
 match=re.search(r'FAIL 0 \| WARN 0 \| SKIP 0 \| PASS (\d+)',tests);assert match
 m['engineering_checks']={'local_R_CMD_check':'Status: OK','assertions_passed':int(match.group(1)),'seed':reg['seed'],'versions':m['versions'],'input_hashes':{'R/scp_power.R':sha('R/scp_power.R'),'tests/testthat/test-scp_power.R':sha('tests/testthat/test-scp_power.R'),'renv_lock':sha('renv.lock')},'log_hashes':{'tests':sha(out/'tests.log'),'check':sha(out/'R-check.log')}}
 duration=re.search(r'Duration: ([0-9.]+) s',tests)
 if duration:m['engineering_checks']['tests_wall_seconds']=float(duration.group(1))
 born=int(subprocess.check_output(['stat','-c','%W',str(out/'R-check.log')],text=True).strip())
 m['engineering_checks']['local_check_wall_seconds']=max(0,(out/'R-check.log').stat().st_mtime-born)
 m['engineering_checks']['local_check_wall_definition']='Check log creation-to-final-modification time; filesystem birth time has one-second precision.'
 m['engineering_checks']['test_seed_record']='Individual test seeds are explicit in the hashed test source; all previous tests are retained.'
paths=[p for p in out.iterdir() if p.is_file() and p.name!='run-manifest.json']
m['output_sha256']={p.name:sha(p) for p in sorted(paths)}
scripts=['inst/scripts/m4b_helpers.R','inst/scripts/prepare_m4b.R','inst/scripts/run_m4b_point.R','inst/scripts/run_m4b_cluster.R','inst/scripts/run_m4b.sbatch','inst/scripts/report_m4b.R','inst/scripts/finalize_m4b_manifest.py','inst/scripts/reproduce_tmt_round.sh','inst/scripts/read_m3_grid.R','inst/scripts/render_reading_grid.py','inst/scripts/consolidate_m4b_layouts.py','inst/scripts/package_tmt_round_evidence.py']
m['reproduction_code_sha256']={p:sha(p) for p in scripts}
m['finalization']={'seed':reg['seed'],'versions':m['versions'],'input_hashes':{'registration':sha(out/'registration.json'),'predictions_frozen':sha(out/'predictions-frozen.json')},'wall_seconds':time.monotonic()-start}
(out/'run-manifest.json').write_text(json.dumps(m,indent=2)+'\n')
print('M4 manifest audited:',m['overall_verdict'],'141 counted points, 28,200 predictions and 28,200 observations')
