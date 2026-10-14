"""Archive old TMT rows, compare cheapest designs and publish the combined grid."""
from pathlib import Path
import csv,json,hashlib,time,shutil,subprocess
start=time.monotonic();out=Path('reports/M3b');sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
read=lambda p:list(csv.DictReader(open(p)))
old=read(out/'archive-v0.2.0/power-grid.csv');new=read(out/'power-grid.csv')
assert len(new)==945
key=lambda r:(r['dataset'],r['n_per_condition'],r['m_cells_per_sample'],r['beta'],r['k'])
for a,b in zip(old,new):
 assert key(a)==key(b)
 if a['dataset'] not in ['specht2019v3','leduc2022_pSCoPE']:
  for field in a:
   try: assert abs(float(a[field])-float(b[field]))<1e-12 or a[field]==b[field]
   except ValueError: assert a[field]==b[field],(field,a[field],b[field])
for file in ['power-grid.csv','realised-fdr.csv','instrument-hours.csv','cheapest-designs.csv']:
 rows=read(out/'archive-v0.2.0'/file)
 rows=[r for r in rows if r['dataset'] in ['specht2019v3','leduc2022_pSCoPE']]
 with open(out/'archive-v0.2.0'/('TMT-'+file),'w') as f:
  w=csv.DictWriter(f,fieldnames=rows[0],lineterminator="\n");w.writeheader();w.writerows(rows)
a=read(out/'archive-v0.2.0/cheapest-designs.csv');b=read(out/'cheapest-designs.csv')
diffs=[]
for before,after in zip(a,b):
 assert all(before[k]==after[k] for k in ['dataset','beta','k'])
 changed=[k for k in before if before[k]!=after[k]]
 diffs.append(dict(dataset=after['dataset'],beta=after['beta'],k=after['k'],changed=bool(changed),changed_fields=';'.join(changed),before_status=before['status'],after_status=after['status'],before_n=before['n_per_condition'],after_n=after['n_per_condition'],before_m=before['m_cells_per_sample'],after_m=after['m_cells_per_sample'],before_hours=before['instrument_hours'],after_hours=after['instrument_hours']))
with open(out/'cheapest-design-diff.csv','w') as f:
 w=csv.DictWriter(f,fieldnames=diffs[0],lineterminator="\n");w.writeheader();w.writerows(diffs)
reached=lambda rows:{key(r) for r in rows if float(r['power'])>=.8}
added=reached(new)-reached(old);removed=reached(old)-reached(new)
change={'cheapest_rows_changed':sum(r['changed'] for r in diffs),'designs_newly_reaching_0_8':[list(x) for x in sorted(added)],'designs_no_longer_reaching_0_8':[list(x) for x in sorted(removed)],'no_adjustment':True}
(out/'target-design-change.json').write_text(json.dumps(change,indent=2)+'\n')
p=out/'power-report.md'
p.write_text(p.read_text().split('\n## Cheapest-design difference',1)[0].replace('](replicate-metrics.csv.gz)','](../M3b/replicate-metrics.csv.gz)').replace('](run-manifest.json)','](../M3b/run-manifest.json)'))
with p.open('a') as f:
 f.write('\n## Cheapest-design difference\n\n')
 f.write(f"{change['cheapest_rows_changed']} cheapest-design rows changed; {len(added)} designs newly reach 0.8 and {len(removed)} designs no longer reach 0.8. No quantity was adjusted. See [the full diff](cheapest-design-diff.csv).\n")
 f.write('\nM3b reran 630 TMT points with 200 replicates and the original M3 seeds. The combined table reuses 315 non-TMT points without rerunning them. Their conditional-power columns are unavailable because the v0.2.0 run did not record that quantity; worked examples have separate current-package manifests. The original TMT rows are [archived](archive-v0.2.0/TMT-power-grid.csv).\n')
# Replace the current M3 tables and plots while retaining the original snapshot.
for name in ['power-grid.csv','realised-fdr.csv','instrument-hours.csv','cheapest-designs.csv','power-curves.png','power-tmt-sensitivity.png','testability-curves.png','testability-tmt-sensitivity.png','power-report.md','reading-grid-diagnostics.csv','reading-grid-small-effects.csv','reading-grid-manifest.json','reading-grid.pdf','reading-grid-render-manifest.json']:
 shutil.copyfile(out/name,Path('reports/M3')/name)
# Current M3 numeric tables have a new provenance pointer, original receipt archived.
Path('reports/M3/README.md').write_text('# Current registered grid\n\nTMT rows use amendment I; non-TMT rows are unchanged. See [M3b report](../M3b/power-report.md), [combined grid](power-grid.csv), [M3b manifest](../M3b/run-manifest.json), and [original archive](../M3b/archive-v0.2.0/).\n')
m=json.loads((out/'run-manifest.json').read_text());m['target_design_change']=change
point_manifests=json.loads(__import__('gzip').open(out/'grid-point-manifests.json.gz','rt').read())
job=point_manifests['1']['cluster']['job']
accounting=subprocess.check_output(['sacct','-j',job,'--format=JobIDRaw,State,ElapsedRaw,TotalCPU,MaxRSS,Submit,Start,End','-P'],text=True)
(out/'cluster-accounting.txt').write_text(accounting)
rows=list(csv.DictReader(accounting.splitlines(),delimiter='|'))
allocation=next(r for r in rows if r['JobIDRaw']==job)
assert allocation['State']=='COMPLETED',allocation
m['cluster_execution']={'job':job,'cores':32,'memory_gib':64,'allocation_elapsed_seconds':int(allocation['ElapsedRaw']),'accounting_sha256':sha(out/'cluster-accounting.txt')}
m['input_hashes']['original_grid']=sha(out/'archive-v0.2.0/power-grid.csv')
m['output_hashes']={p.name:sha(p) for p in out.iterdir() if p.is_file() and p.name not in ['run-manifest.json','cluster-allocation.log']}
m['finalization_wall_seconds']=time.monotonic()-start
Path('reports/M3/current-grid-manifest.json').write_text(json.dumps({'versions':m['versions'],'seed':m['seed'],'current_numeric_results_manifest':'../M3b/run-manifest.json','original_manifest':'../M3b/archive-v0.2.0/run-manifest.json','non_TMT_results_unchanged':True,'TMT_endpoint_masks_updated':True},indent=2)+'\n')
(out/'run-manifest.json').write_text(json.dumps(m,indent=2)+'\n')
print(json.dumps(change))
