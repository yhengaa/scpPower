"""Finalize revision evidence after the fixed gate and package validation."""
from pathlib import Path
import hashlib,json,subprocess,time
start=time.monotonic();out=Path('reports/M2-revised');path=out/'run-manifest.json'
m=json.loads(path.read_text());sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
m['engineering_validation_scope']='M2 revision package 0.1.2: 210 assertions and clean R CMD check before M3 implementation. Whole-package 0.2.0 validation is recorded separately in M3 (250 assertions). Reproduction-code hashes record the final repository snapshot.'
m['metric_definitions']['half_cells_fraction']='Descriptive only, outside the gate: fraction of all proteins detected in at least half the native-design cells; observed versus mean and SD over the same twenty replicates.'
m['startup_correction']=json.loads((out/'startup-correction.json').read_text())
m['unchanged_calibration_evidence_sha256']={str(p):sha(p) for p in Path('reports/M1-final').iterdir() if p.is_file() and p.suffix in ['.csv','.json','.png','.pdf','.md']}
m['parameter_generation_manifest']=json.loads((out/'preset-manifest.json').read_text())
m['fallback_precedence']='Own finite component first. Structurally unidentified batch is zero under G1; a missing protein fit otherwise uses tier component. Unmapped donor variance uses tier component. Fallback indicators in the parameter CSV mark originally unavailable estimates before the explicit unidentified-batch zero override.'
m['reproduction_reference']='Original M2 at b34d9d4; accepted M1 at v0.1.0. Revised scripts download the same registered public inputs through the M1 reproduction chain, then generate empirical presets and run the unchanged gate.'
m['revision_finalization']=dict(seed=20261013,versions=m['versions'],input_hashes={'preset_manifest':sha(out/'preset-manifest.json'),'startup_correction':sha(out/'startup-correction.json')},wall_seconds=time.monotonic()-start)
path.write_text(json.dumps(m,indent=2)+'\n')
subprocess.run(['python','inst/scripts/finalize_stage_manifest.py','M2-revised'],check=True)
