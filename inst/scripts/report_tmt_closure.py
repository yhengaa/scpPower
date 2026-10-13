"""Report accepted M4 diagnostics and K without altering any verdict or number."""
from pathlib import Path
import csv,json,hashlib,time,math
start=time.monotonic();root=Path('reports/J')
rows=list(csv.DictReader(Path('reports/M4/validation-table.csv').open()));summary=[]
for dataset in sorted(set(r['dataset'] for r in rows)):
 ds=[r for r in rows if r['dataset']==dataset and r['counted']=='TRUE']
 high=[r for r in ds if max(float(r['predicted_power']),float(r['observed_power']))>=.2]
 def gap(rs):
  a=[abs(float(r['predicted_conditional_power'])-float(r['observed_conditional_power'])) for r in rs if r['predicted_conditional_power']!='NA' and r['observed_conditional_power']!='NA']
  return max(a) if a else None
 conditional=[float(r['predicted_conditional_power'])-float(r['observed_conditional_power']) for r in ds if r['predicted_conditional_power']!='NA' and r['observed_conditional_power']!='NA']
 summary.append(dict(dataset=dataset,counted_points=len(ds),high_power_points=len(high),high_power_within_10_pp=sum(r['within_10_pp']=='TRUE' for r in high),mean_signed_high_gap=sum(float(r['predicted_power'])-float(r['observed_power']) for r in high)/len(high) if high else None,mean_signed_conditional_gap=sum(conditional)/len(conditional) if conditional else None,max_conditional_gap_all=gap(ds),max_conditional_gap_high=gap(high),max_unconditional_gap=max((float(r['absolute_gap']) for r in ds),default=None)))
versions=json.loads(Path('reports/M4/run-manifest.json').read_text())['versions']
for r in summary:r['software_versions']='; '.join(f'{k}={v}' for k,v in versions.items())
with root.joinpath('v0.3.0-power-window-diagnostics.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(summary[0]));w.writeheader();w.writerows(summary)
limitation=('In carrier-TMT presets, unconditional power is validated only up to observed power of about 0.3. At twofold, v0.3.0 predictions ran above observation by up to 20 percentage points, associated with over-predicted testability. Q Exactive and timsTOF could not be validated at useful power for lack of cells. No TMT design numbers or new design examples are presented in the README or a vignette during this round.')
correction=('The 1-pp conditional-power description holds for mean signed TMT gaps: -0.82 pp for Specht and +0.94 pp for Leduc across counted points. Brunner averages +1.30 pp and Q Exactive -0.46 pp. It is not a bound on every point: among points with maximum power at least 0.2, maximum conditional gaps are 1.53 pp for Specht, 4.56 pp for Leduc and 6.05 pp for Brunner. Across all counted points the maxima are 7.77, 11.13, 7.57 and 2.32 pp respectively. The high-power agreement counts do reproduce: 3/8 Specht, 11/14 Leduc and 12/12 Brunner; Q Exactive has none. Testability is the main TMT discrepancy, while the statement that the entire gap is testability with universal agreement within 1 pp is stronger than the recorded columns support. This clarification does not change any data, criterion or accepted verdict.')
original=Path('reports/M4/validation-report-v0.3.0.md');current=Path('reports/M4/validation-report.md')
if not original.exists():original.write_bytes(current.read_bytes())
addition='\n\n## M4 closure and limitation K — addendum v1.1.5\n\n'+limitation+'\n\n'+correction+'\n\n[Reproducible power-window and conditional-gap table](../J/v0.3.0-power-window-diagnostics.csv). The original v0.3.0 report is preserved in validation-report-v0.3.0.md; its original manifest remains the receipt for that snapshot. No M4 numerical result or verdict is changed.\n'
current.write_text(original.read_text()+addition)
p=Path('README.md');text=p.read_text()
if '```r\npower <- scp_power(' in text:
 begin=text.index('```r\npower <- scp_power(');end=text.index('```',begin+4)+3
 text=text[:begin]+'See the documented `scp_power()` return values in [the function documentation](man/scp_power.Rd).'+text[end:]
text=text.split('\n\nCarrier-TMT validation limitation (addendum v1.1.5):',1)[0]
text+='\n\nCarrier-TMT validation limitation (addendum v1.1.5):\n\n'+limitation+'\n\n'+correction+'\n\n[Diagnostic J](reports/J/diagnostic-report.md) determines the assay-consistency branch before any preset change. The registered v0.3.0 M4 verdict stands. No M5 work is included.\n'
p.write_text(text)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
root.joinpath('closure-manifest.json').write_text(json.dumps({'seed':20261013,'versions':versions,'input_hashes':{'original_M4_table':sha('reports/M4/validation-table.csv'),'original_M4_report':sha(original),'script':sha(__file__)},'output_hashes':{'validation_report':sha(current),'README':sha(p),'power_window_table':sha(root/'v0.3.0-power-window-diagnostics.csv')},'registered_M4_verdict_unchanged':True,'wall_seconds':time.monotonic()-start},indent=2)+'\n')
root.joinpath('limitation-K.txt').write_text(limitation+'\n\n'+correction+'\n')
print('K documented with verified conditional-gap correction; no original M4 numbers changed')
