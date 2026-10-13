"""Complete J and M2b audit receipts without modifying a scientific estimate."""
from pathlib import Path
import json,hashlib,time,csv,subprocess
start=time.monotonic();sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
for name in ['J','M2b']:
 root=Path('reports')/name;m=json.loads((root/'run-manifest.json').read_text())
 m['mandate_addendum']='1.1.5';m['output_sha256']={p.name:sha(p) for p in sorted(root.iterdir()) if p.is_file() and p.name!='run-manifest.json'}
 if name=='J':
  m['report_before_preset_change']=True;m['pipeline_audit_manifest_sha256']=sha(root/'pipeline-audit-manifest.json');m['closure_manifest_sha256']=sha(root/'closure-manifest.json')
  m['diagnostic_depth_definition']='Twenty fixed-validation-protein, native-layout G2 detection-only replicates per mask; independent p/q draws. Descriptive, not gated.'
  m['diagnostic_run_definition']='One comparison experiment per dataset; replicate seeds and hashes recorded within it, wall time recorded for each dataset experiment.'
  from datetime import datetime,timezone
  report_time=(root/'report-manifest.json').stat().st_mtime
  preset_time=Path('reports/M2b/preset-manifest.json').stat().st_mtime
  assert report_time<preset_time,'J report did not precede I preset generation'
  m['report_completed_UTC']=datetime.fromtimestamp(report_time,timezone.utc).isoformat()
  m['I_preset_generation_completed_UTC']=datetime.fromtimestamp(preset_time,timezone.utc).isoformat()
  m['ordering_time_provenance']='Original on-disk report and first-I preset manifest completion times, recorded before publication; J report content hash is also embedded in each I preset.'
  m['diagnostic_manifest_before_I_sha256']=sha(root/'diagnostic-manifest-before-I.json')
 else:
  summary=list(csv.DictReader((root/'fidelity-summary.csv').open()));m['both_M2b_pass']=len(summary)==2 and all(r['status']=='PASS' for r in summary)
  m['cluster_accounting_sha256']=sha(root/'cluster-accounting.txt');m['diagnostic_J_report_sha256']=sha('reports/J/diagnostic-report.md')
 m['audit']={'seed':m['seed'],'versions':m['versions'],'input_hashes':{'script':sha(__file__)},'wall_seconds':time.monotonic()-start}
 (root/'run-manifest.json').write_text(json.dumps(m,indent=2)+'\n')
print('J and M2b manifests completed; quantitative gate results unchanged')
