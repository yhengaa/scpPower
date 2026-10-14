# scpPower

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

### Carrier TMT: Leduc

For a twofold change, compare eight and twelve samples per condition, each with thirty-two cells.

```r
x <- scp_power("leduc2022_pSCoPE", n_per_condition = c(8, 12),
  m_cells_per_sample = 32, beta = log2(2), replicates = 200, seed = 20261013)
x$results[, c("n_per_condition", "m_cells_per_sample", "power", "testability",
  "conditional_power", "realised_fdr", "instrument_hours")]
```

| Samples/condition | Cells/sample | Power | Testability | Conditional power | Realised FDR | Hours |
| --- | --- | --- | --- | --- | --- | --- |
| 8 | 32 | 0.348 | 0.398 | 0.864 | 0.0385 | 58.58 |
| 12 | 32 | 0.364 | 0.386 | 0.949 | 0.0372 | 87.08 |

### plexDIA: timsTOF

For a twofold change, compare eight and twelve samples per condition, each with thirty-two cells.

```r
x <- scp_power("derks2022_timsTOF", n_per_condition = c(8, 12),
  m_cells_per_sample = 32, beta = log2(2), replicates = 200, seed = 20261013)
x$results[, c("n_per_condition", "m_cells_per_sample", "power", "testability",
  "conditional_power", "realised_fdr", "instrument_hours")]
```

| Samples/condition | Cells/sample | Power | Testability | Conditional power | Realised FDR | Hours |
| --- | --- | --- | --- | --- | --- | --- |
| 8 | 32 | 0.777 | 0.898 | 0.854 | 0.0365 | 128.00 |
| 12 | 32 | 0.868 | 0.891 | 0.960 | 0.0391 | 192.00 |

### Label-free: Brunner

For a twofold change, compare eight and twelve samples per condition, each with thirty-two cells.

```r
x <- scp_power("brunner2022", n_per_condition = c(8, 12),
  m_cells_per_sample = 32, beta = log2(2), replicates = 200, seed = 20261013)
x$results[, c("n_per_condition", "m_cells_per_sample", "power", "testability",
  "conditional_power", "realised_fdr", "instrument_hours")]
```

| Samples/condition | Cells/sample | Power | Testability | Conditional power | Realised FDR | Hours |
| --- | --- | --- | --- | --- | --- | --- |
| 8 | 32 | 0.647 | 0.673 | 0.905 | 0.0478 | 298.67 |
| 12 | 32 | 0.670 | 0.647 | 0.972 | 0.0481 | 448.00 |

## Cheapest registered design reaching 0.8 power

The combined [M3b grid](reports/M3b/power-grid.csv) contains updated TMT rows and unchanged non-TMT rows. The table below selects the least instrument time among designs whose mean simulated power reaches the target; it does not extrapolate beyond the registered grid. See the [report](reports/M3b/power-report.md) and [cheapest-design difference](reports/M3b/cheapest-design-diff.csv).

| Preset | Fold change | k | Samples/condition | Cells/sample | Power | Hours | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| specht2019v3 | 1.25 | 1.00 | — | — | — | — | no registered design reaches 0.8 |
| specht2019v3 | 1.25 | 0.75 | — | — | — | — | no registered design reaches 0.8 |
| specht2019v3 | 1.25 | 0.50 | — | — | — | — | no registered design reaches 0.8 |
| specht2019v3 | 1.50 | 1.00 | — | — | — | — | no registered design reaches 0.8 |
| specht2019v3 | 1.50 | 0.75 | — | — | — | — | no registered design reaches 0.8 |
| specht2019v3 | 1.50 | 0.50 | — | — | — | — | no registered design reaches 0.8 |
| specht2019v3 | 2.00 | 1.00 | — | — | — | — | no registered design reaches 0.8 |
| specht2019v3 | 2.00 | 0.75 | — | — | — | — | no registered design reaches 0.8 |
| specht2019v3 | 2.00 | 0.50 | — | — | — | — | no registered design reaches 0.8 |
| leduc2022_pSCoPE | 1.25 | 1.00 | — | — | — | — | no registered design reaches 0.8 |
| leduc2022_pSCoPE | 1.25 | 0.75 | — | — | — | — | no registered design reaches 0.8 |
| leduc2022_pSCoPE | 1.25 | 0.50 | — | — | — | — | no registered design reaches 0.8 |
| leduc2022_pSCoPE | 1.50 | 1.00 | — | — | — | — | no registered design reaches 0.8 |
| leduc2022_pSCoPE | 1.50 | 0.75 | — | — | — | — | no registered design reaches 0.8 |
| leduc2022_pSCoPE | 1.50 | 0.50 | — | — | — | — | no registered design reaches 0.8 |
| leduc2022_pSCoPE | 2.00 | 1.00 | — | — | — | — | no registered design reaches 0.8 |
| leduc2022_pSCoPE | 2.00 | 0.75 | — | — | — | — | no registered design reaches 0.8 |
| leduc2022_pSCoPE | 2.00 | 0.50 | — | — | — | — | no registered design reaches 0.8 |
| derks2022_timsTOF | 1.25 | 1.00 | — | — | — | — | no registered design reaches 0.8 |
| derks2022_timsTOF | 1.50 | 1.00 | — | — | — | — | no registered design reaches 0.8 |
| derks2022_timsTOF | 2.00 | 1.00 | 12 | 32 | 0.865 | 192.00 | Target reached in registered grid |
| derks2022_QExactive | 1.25 | 1.00 | — | — | — | — | no registered design reaches 0.8 |
| derks2022_QExactive | 1.50 | 1.00 | — | — | — | — | no registered design reaches 0.8 |
| derks2022_QExactive | 2.00 | 1.00 | — | — | — | — | no registered design reaches 0.8 |
| brunner2022 | 1.25 | 1.00 | — | — | — | — | no registered design reaches 0.8 |
| brunner2022 | 1.50 | 1.00 | — | — | — | — | no registered design reaches 0.8 |
| brunner2022 | 2.00 | 1.00 | — | — | — | — | no registered design reaches 0.8 |

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
4. Cell-depth dependence is not modelled. Real depth SD is roughly 5 to 7 times simulated SD (recorded endpoint ratios: 5.17 and 7.50); this did not affect the registered validation verdict. See [Diagnostic J](reports/J/native-depth-summary.csv).
5. Q Exactive and timsTOF could not be validated at useful power for lack of public cells. Q Exactive passed only at low power; timsTOF was not evaluable.
6. The primary endpoint tests per-sample means with limma-trend. Other tests will give other numbers.

## Validation

| Stage | Evaluated datasets | Result | Evidence |
| --- | --- | --- | --- |
| M2 (amendment G) | 5 | 5/5 PASS | [Fidelity report](reports/M2-revised/README.md) |
| M4, v0.3.0 | 4 | specht2019v3: 67/72; leduc2022_pSCoPE: 66/69; derks2022_QExactive: 6/6; brunner2022: 57/57; overall PASS | [Validation report](reports/M4/validation-report.md) |
| M4b, v0.3.1 | 2 | specht2019v3: 72/72; leduc2022_pSCoPE: 69/69; both PASS | [Validation report](reports/M4b/validation-report.md); [M2b](reports/M2b/README.md) |

M4 validates cell-level noise, missingness and the registered test. It cannot validate donor-variance transfer or k. The [vignette](vignettes/designing-scp.Rmd) explains testability, conditional power and compression sensitivity. Every reported data estimate has a public-input reproduction script in [inst/scripts](inst/scripts), including [release reproduction](inst/scripts/reproduce_release.sh); raw input hashes are retained in the stage manifests.

## Development

Code was written with the assistance of GPT (OpenAI) and independently reviewed with Claude (Anthropic), both under the author's direction. All scientific decisions, registered criteria, data choices and verifications are the author's.

## Citation and licence

Cite Youshan Heng, *scpPower*, version 1.0.0, using [CITATION.cff](CITATION.cff). scpPower is released under [GNU GPL version 3](https://www.gnu.org/licenses/gpl-3.0.html).
