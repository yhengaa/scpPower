"""Generate README numeric outputs exclusively from recorded package runs."""
from pathlib import Path
import csv,json,time,hashlib
start=time.monotonic();read=lambda p:list(csv.DictReader(open(p)));sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def table(rows,columns,labels=None):
 labels=labels or columns
 out=['| '+' | '.join(labels)+' |','| '+' | '.join(['---']*len(columns))+' |']
 for r in rows:out.append('| '+' | '.join(str(r[c]) for c in columns)+' |')
 return '\n'.join(out)
def num(x,d=3):
 return 'NA' if x in ['NA','nan',''] else f'{float(x):.{d}f}'
text='''# scpPower

scpPower estimates how many biological samples and cells per sample a two-condition single-cell proteomics experiment needs. Specify the fold change, instrument and labelling scheme, and compare power under an explicitly calibrated model. Each design also reports nominal instrument hours.

[![R CMD check](https://github.com/yhengaa/scpPower/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yhengaa/scpPower/actions/workflows/R-CMD-check.yaml)

## Install

```r
install.packages("remotes")
remotes::install_github("yhengaa/scpPower", ref = "v1.0.0")
library(scpPower)
```

For the recorded software environment, clone the repository and optionally run `renv::restore()` with the committed lockfile. Version differences warn and the package continues.

## Worked examples

The following calls ask whether more biological samples help detect a twofold change with a fixed number of cells per sample. Outputs are Monte Carlo estimates, reproduced by [run_release_examples.R](inst/scripts/run_release_examples.R) with [seeds, hashes, versions and wall times](reports/M5/worked-examples-manifest.json). These demonstration runs are separate from the registered grid. Instrument hours exclude unreported overhead.
'''
for name,label in [('leduc2022_pSCoPE','Carrier TMT: Leduc'),('derks2022_timsTOF','plexDIA: timsTOF'),('brunner2022','Label-free: Brunner')]:
 rows=read(f'reports/M5/{name}-worked-example.csv')
 text+=f'\n### {label}\n\nFor a twofold change, compare eight and twelve samples per condition, each with thirty-two cells.\n\n```r\nx <- scp_power("{name}", n_per_condition = c(8, 12),\n  m_cells_per_sample = 32, beta = log2(2), replicates = 200, seed = 20261013)\nx$results[, c("n_per_condition", "m_cells_per_sample", "power", "testability",\n  "conditional_power", "realised_fdr", "instrument_hours")]\n```\n\n'
 for r in rows:
  for c in ['power','testability','conditional_power']:r[c]=num(r[c])
  r['realised_fdr']=num(r['realised_fdr'],4);r['instrument_hours']=num(r['instrument_hours'],2)
 text+=table(rows,['n_per_condition','m_cells_per_sample','power','testability','conditional_power','realised_fdr','instrument_hours'],['Samples/condition','Cells/sample','Power','Testability','Conditional power','Realised FDR','Hours'])+'\n'
text+='''
## Cheapest registered design reaching 0.8 power

The combined [M3b grid](reports/M3b/power-grid.csv) contains updated TMT rows and unchanged non-TMT rows. The table below selects the least instrument time among designs whose mean simulated power reaches the target; it does not extrapolate beyond the registered grid. See the [report](reports/M3b/power-report.md) and [cheapest-design difference](reports/M3b/cheapest-design-diff.csv).

'''
rows=read('reports/M3b/cheapest-designs.csv')
for r in rows:
 r['fold_change']=f"{2**float(r['beta']):.2f}";r['k']=num(r['k'],2)
 if r['n_per_condition']=='NA':
  r['status']='no registered design reaches 0.8'
  for c in ['n_per_condition','m_cells_per_sample','power','instrument_hours']:r[c]='—'
 else:r['power']=num(r['power']);r['instrument_hours']=num(r['instrument_hours'],2)
text+=table(rows,['dataset','fold_change','k','n_per_condition','m_cells_per_sample','power','instrument_hours','status'],['Preset','Fold change','k','Samples/condition','Cells/sample','Power','Hours','Status'])+'\n'
text+='''
## Calibrate your own data

Use an unimputed QFeatures protein assay and explicit biological-sample, condition and batch metadata. Here `pilot` has an already-log2 `proteins` assay; `donor_means` is a matching protein-by-donor matrix from the reference condition and `donor_info$donor` identifies its columns. Run duration, usable single-cell slots and carrier ratio are measured pilot inputs. In carrier TMT, `raw_proteins` supplies aligned uncentred abundance values for the intensity curve while the preset's detection mask comes from the endpoint assay.

```r
cal <- scp_calibrate(pilot, assay = "proteins", log2_input = TRUE,
  sample_col = "sample", condition_col = "condition", batch_col = "plex",
  bulk = donor_means, bulk_metadata = donor_info,
  normalize_cells = TRUE, missingness_input = raw_proteins,
  carrier_ratio = pilot_carrier_ratio, compression_k = pilot_k, seed = 20261013)
source(system.file("scripts", "pilot_to_preset.R", package = "scpPower"))
preset <- pilot_to_preset(cal,
  SummarizedExperiment::assay(pilot[["proteins"]]),
  as.data.frame(SummarizedExperiment::colData(pilot[["proteins"]])), "plex",
  cells_per_run = pilot_single_cell_slots, minutes_per_run = pilot_run_minutes,
  carrier = TRUE)
sim <- scp_simulate(NULL, 4, 20, beta = log2(2), preset = preset, seed = 20261013)
power <- scp_power(preset, 4, 20, beta = log2(2), seed = 20261013)
```

Inspect component reason columns and overlap before using a preset; unavailable donor or cell variances are reported rather than invented. The [pilot protocol](inst/doc/pilot-protocol.md) explains controls and label-free inputs. The conversion script uses the existing parameter fallbacks and endpoint-mask rule; it fits no additional model.

## Limitations

1. Between-donor variance is mapped from a whole-PBMC bulk cohort. Transfer to single cells, including composition differences, is an assumption tested only by a pilot in the intended cell type.
2. Carrier compression k is an assumption: default 1, with sensitivities 0.75 and 0.5. A known-ratio standard is needed to estimate it for a lab.
3. TMT predictions are validated on real cells only up to observed power of about 0.35. In v0.3.0, twofold predictions exceeded observation by up to 20 percentage points through testability; amendment I corrected endpoint-mask consistency. The accepted correction to the conditional-power statement is recorded in the [TMT closure report](reports/J/diagnostic-report.md).
'''
depth=[r for r in read('reports/J/native-depth-summary.csv') if r['mask']=='aggregated']
ratios=[float(r['observed_depth_sd'])/float(r['simulated_mean_depth_sd']) for r in depth]
text+=f'4. Cell-depth dependence is not modelled. Real depth SD is roughly 5 to 7 times simulated SD (recorded endpoint ratios: {min(ratios):.2f} and {max(ratios):.2f}); this did not affect the registered validation verdict. See [Diagnostic J](reports/J/native-depth-summary.csv).\n'
text+='''5. Q Exactive and timsTOF could not be validated at useful power for lack of public cells. Q Exactive passed only at low power; timsTOF was not evaluable.
6. The primary endpoint tests per-sample means with limma-trend. Other tests will give other numbers.

## Validation

'''
m2=read('reports/M2-revised/fidelity-summary.csv');m4=read('reports/M4/dataset-verdicts.csv');m4b=read('reports/M4b/dataset-verdicts.csv')
summary=[{'stage':'M2 (amendment G)','datasets':len(m2),'result':f"{sum(r['status']=='PASS' for r in m2)}/{len(m2)} PASS",'report':'[Fidelity report](reports/M2-revised/README.md)'},
 {'stage':'M4, v0.3.0','datasets':sum(r['status']!='NOT EVALUABLE' for r in m4),'result':'; '.join(f"{r['dataset']}: {r['within_10_pp_points']}/{r['counted_points']}" for r in m4 if r['status']=='PASS')+'; overall PASS','report':'[Validation report](reports/M4/validation-report.md)'},
 {'stage':'M4b, v0.3.1','datasets':len(m4b),'result':'; '.join(f"{r['dataset']}: {r['within_10_pp_points']}/{r['counted_points']}" for r in m4b)+'; both PASS','report':'[Validation report](reports/M4b/validation-report.md); [M2b](reports/M2b/README.md)'}]
text+=table(summary,['stage','datasets','result','report'],['Stage','Evaluated datasets','Result','Evidence'])+'\n'
text+='''
M4 validates cell-level noise, missingness and the registered test. It cannot validate donor-variance transfer or k. The [vignette](vignettes/designing-scp.Rmd) explains testability, conditional power and compression sensitivity. Every reported data estimate has a public-input reproduction script in [inst/scripts](inst/scripts), including [release reproduction](inst/scripts/reproduce_release.sh); raw input hashes are retained in the stage manifests.

## Development

Code was written with the assistance of GPT (OpenAI) and independently reviewed with Claude (Anthropic), both under the author's direction. All scientific decisions, registered criteria, data choices and verifications are the author's.

## Citation and licence

Cite Youshan Heng, *scpPower*, version 1.0.0, using [CITATION.cff](CITATION.cff). scpPower is released under [GNU GPL version 3](https://www.gnu.org/licenses/gpl-3.0.html).
'''
Path('README.md').write_text(text)
m=json.loads(Path('reports/M5/worked-examples-manifest.json').read_text())
inputs=['reports/M5/worked-examples.csv','reports/M3b/cheapest-designs.csv','reports/J/native-depth-summary.csv','reports/M2-revised/fidelity-summary.csv','reports/M4/dataset-verdicts.csv','reports/M4b/dataset-verdicts.csv']
Path('reports/M5/README-generation-manifest.json').write_text(json.dumps({'seed':m['seed'],'versions':m['versions'],'input_hashes':{p:sha(Path(p)) for p in inputs},'output_sha256':sha(Path('README.md')),'wall_seconds':time.monotonic()-start},indent=2)+'\n')
