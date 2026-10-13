"""Add scheduler receipts and finalize the fully completed M3 evidence."""
from pathlib import Path
import json,hashlib,subprocess,time,csv,datetime
start=time.monotonic();out=Path('reports/M3');p=out/'run-manifest.json'
m=json.loads(p.read_text());sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
cmd=['sacct','-j','5491088,5491758,5492391','--format=JobIDRaw,State,ElapsedRaw,TotalCPU,MaxRSS,Submit,Start,End','-P']
s=subprocess.run(cmd,text=True,capture_output=True,check=True);(out/'cluster-accounting.txt').write_text(s.stdout)
rows=list(csv.DictReader(s.stdout.splitlines(),delimiter='|'))
live=[r for r in rows if r['JobIDRaw']=='5492391']
if len(live)!=1 or live[0]['State']!='COMPLETED':raise RuntimeError('Interactive cluster allocation is not complete')
m['cluster_execution']=dict(initial_array_job='5491088',initial_array_retained_points='Complete point manifests retained; cancelled pending points were executed with the same seeds in the live allocation',grouped_array_job='5491758',grouped_array_status='Cancelled while pending to use a live allocation',live_job='5492391',live_node='nc10107',live_cores=16,live_memory_gib=32,live_allocation_elapsed_seconds=int(live[0]['ElapsedRaw']),accounting_file='cluster-accounting.txt',accounting_sha256=sha(out/'cluster-accounting.txt'),scheduling_only_changes=True)
m['input_hashes']['registered_power_code']=m['grid_registration']['input_hashes']['power']
m['input_hashes']['registered_simulator_code']=m['grid_registration']['input_hashes']['simulator']
assert sha(Path('R/scp_power.R'))==m['input_hashes']['registered_power_code']
assert sha(Path('R/scp_simulate.R'))==m['input_hashes']['registered_simulator_code']
m['native_fidelity_replay_verification']=json.loads((out/'M2-replay-verification.json').read_text())
m['instrument_time_provenance']=json.loads((out/'instrument-preset-manifest.json').read_text())
m['replicate_file_sha256']=sha(out/'replicate-metrics.csv.gz')
m['grid_point_manifest_sha256']=sha(out/'grid-point-manifests.json.gz')
m['report_finalization']=dict(seed=20261013,versions=m['versions'],input_hashes={'scheduler':sha(out/'cluster-accounting.txt'),'report':sha(out/'power-report.md')},wall_seconds=time.monotonic()-start)
p.write_text(json.dumps(m,indent=2)+'\n')
subprocess.run(['python','inst/scripts/finalize_stage_manifest.py','M3'],check=True)
# Binary evidence is also checksummed in output_sha256.
m=json.loads(p.read_text());m['output_sha256'].update({x.name:sha(x) for x in out.glob('*.gz')});p.write_text(json.dumps(m,indent=2)+'\n')
print('M3 final manifest and cluster receipt complete.')
