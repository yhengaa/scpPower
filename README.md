# scpPower

[![R CMD check](https://github.com/yhengaa/scpPower/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yhengaa/scpPower/actions/workflows/R-CMD-check.yaml)

Calibration for mass-spectrometry single-cell proteomics experimental design.
The governing [mandate is frozen](scpPower_mandate_v1.1.md). The package provides
`scp_calibrate()` and `scp_simulate()`; calibration status and limitations
are recorded in [the final M1 evidence](reports/M1-final).

M1 closes under [frozen addendum v1.1.2](scpPower_mandate_v1.1.md).
The [final M1 evidence](reports/M1-final) calibrates every registered dataset and
both Derks instruments. TMT missingness uses raw reporter-ion PSMs; its accepted
ratio-based variance is retained. Brunner combines cell and acquisition-run
noise as the identified within-line variance. The earlier stopped and resumed
runs remain in [reports/M1](reports/M1) and [reports/M1-resumed](reports/M1-resumed).

`scp_calibrate()` accepts a protein-by-cell matrix with cell metadata or an
explicit protein assay from a `QFeatures` object. It returns variance components,
abundance quintiles, a fitted quantification curve, a registered carrier
compression assumption, mapped bulk donor variance, and a run manifest. Unavailable quantities
receive explained `NA` values. No calibration values are imputed.

```r
renv::restore()
calibration <- scp_calibrate(
    protein_log2, cell_metadata,
    sample_col = "biological_sample",
    condition_col = "condition",
    batch_col = "plex",
    bulk = bulk_log2, bulk_metadata = donor_metadata
)
calibration$tiers
calibration$manifest
```

For the registered cell-line datasets, `biological_unit="culture"` treats each
cell line as one culture sample and fits line means with shared batch effects.
These units do not identify donor variance. Independent donor datasets retain
the sample random-effect mode. One
cell per run cannot separate acquisition variance from cell noise. Bulk donor
variance transfers only under the registered cell-type assumption, which
requires a lab pilot and cannot be validated by the semi-synthetic M4 gate.

Every calibration output includes R, scp, limma, edgeR and the scplainer workflow
version. Scplainer is implemented in scp; it is not treated as a separate R
package. Version differences warn and do not block execution.

The scripts in [inst/scripts](inst/scripts) download public inputs, record
hashes, preserve seeds and wall time, and regenerate reported results. Use
`bootstrap.R` to prepare the R library on the cluster and `inspect_resumed.R` to
fetch and inspect the four registered datasets. Raw data and local libraries
are excluded from Git. The committed `renv.lock` records the environment.

To reproduce final M1 from public data after restoring the committed environment:

```sh
python inst/scripts/fetch_resumed_sources.py
python inst/scripts/fetch_bulk_metadata.py
Rscript inst/scripts/inspect_resumed.R
Rscript inst/scripts/calibrate_resumed.R
python inst/scripts/extract_paper_benchmarks.py
Rscript inst/scripts/extract_tmt_psms.R
Rscript inst/scripts/close_m1.R
python inst/scripts/extract_paper_benchmarks.py M1-final
Rscript inst/scripts/report_final_m1.R
python inst/scripts/package_m1_evidence.py
```

The raw TMT source and variance inputs are paired by their protein and cell
identifiers, retaining the accepted population without new QC thresholds.
The two measurement scales are named in the manifest. Derks metadata are
matched by physical file and label; its accepted estimates are unchanged.
No estimate is adjusted to match a paper.

Compression is a registered user assumption, never an estimate from these datasets.
The manifest records the default and sensitivity assumptions and the ratio-standard
pilot requirement. Walzik supplies whole-PBMC donor variance, including donor
composition differences; sorted-subset proteomics is unavailable.

`scp_simulate()` generates log2 measurements with missing values, cell metadata
(sample, condition, batch), injected truth and a run manifest. The YAML presets
in [inst/presets](inst/presets) hold the frozen calibration; a user can supply a
YAML file or complete preset list. For example:

```r
simulated <- scp_simulate(
    NULL, n_per_condition = 4, m_cells_per_sample = 16,
    beta = log2(1.5), preset = "derks2022_timsTOF",
    layout = "balanced", k = 1, seed = 20261013
)
simulated$matrix
simulated$cell_metadata
```

The registered [M2 fidelity report](reports/M2/fidelity-report.md) records every
PASS and FAIL, the observed/simulated numbers, and suspected causes. Failed
datasets stop at that report; the next directive decides any change. Calibration
and simulation were not tuned to obtain a pass. Brunner's separate run variance
is absorbed into its cell term. TMT baseline and missingness use raw intensities
while variance retains the accepted ratio scale.

To reproduce the original M2 evidence, use source at
[b34d9d4](https://github.com/yhengaa/scpPower/commit/b34d9d4bd89a4cfb17f8fb1af5e8c37bba6c5f6a)
after regenerating final M1. The revised commands below use the current source.
Original M2 commands:

```sh
Rscript inst/scripts/generate_presets.R
Rscript inst/scripts/fidelity_m2.R
Rscript inst/scripts/report_m2.R
```

Author and maintainer: Youshan Heng,
[ORCID 0000-0002-2758-171X](https://orcid.org/0000-0002-2758-171X),
poggar@hotmail.com (verified through the connected GitHub profile).

The carrier reference identified by PMID 33288958 is Cheung et al.,
[*Nature Methods*](https://doi.org/10.1038/s41592-020-01002-5).
The frozen mandate retains the PI's text verbatim; the bibliographic correction
is recorded here and in calibration documentation.

The approved [addendum v1.1.3](scpPower_mandate_v1.1.md) revises the simulator
under amendment G. Registered presets jointly resample calibrated protein rows
within tiers, retaining each protein's mean, variance components and detection
parameters. Tier summaries remain for reporting. TMT uses a protein's plex
identification probability followed by conditional channel detection; both
probabilities are independent of the effect because the carrier pools both
conditions. DIA and label-free presets keep protein-specific baseline
detection and use the fitted intensity slope for expected-intensity shifts.

The [revised M2 evidence](reports/M2-revised/fidelity-report.md) preserves the
registered tolerances, native designs, replicate count and metric definitions.
It also reports the fraction quantified in at least half the cells, outside
the gate. The [original M2 evidence](reports/M2/fidelity-report.md) remains
available. M3 and the v0.2.0 tag require all series to pass the revised gate.
The M4 registration is unchanged.

After reproducing final M1, regenerate amendment G and its evidence with:

```sh
Rscript inst/scripts/generate_presets_revised.R
Rscript inst/scripts/fidelity_m2_revised.R
Rscript inst/scripts/report_m2_revised.R
```

All registered series pass the revised M2 gate. `scp_power()` now evaluates
the approved M3 endpoint: per-sample log2 means with a minimum of three
quantified cells in every sample, limma-trend, BH and injected-sign agreement.
It reports power, testability, realised FDR, Monte Carlo standard errors,
actual acquisition runs, nominal instrument hours and the cheapest tested
design meeting the requested power target.

```r
power <- scp_power(
    "specht2019v3", n_per_condition = c(3, 4, 6, 8, 12),
    m_cells_per_sample = c(1, 2, 4, 8, 16, 32, 64),
    beta = log2(c(1.25, 1.5, 2)), layout = "balanced",
    k = c(1, 0.75, 0.5), replicates = 200, seed = 20261013
)
power$results
power$cheapest
```

The [M3 report](reports/M3/power-report.md) includes power and testability
curves, the TMT sensitivity assumptions, realised FDR, instrument-hour
estimates, point runtimes and seeds/hashes. A target that no registered
design reaches is reported as unavailable. Instrument hours use cited
method/gradient durations and actual occupied batches; unreported overhead
is excluded. Passing M2 does not establish power validity or bulk-donor
transfer. No M4 evaluation has been performed.

Reproduce M3 after the revised M2 gate passes:

```sh
python inst/scripts/extract_instrument_times.py
Rscript inst/scripts/add_instrument_times.R
Rscript inst/scripts/verify_revised_replay.R
Rscript inst/scripts/register_m3_grid.R
# Choose your cluster account in the submitted file.
sbatch --array=1-27%27 inst/scripts/run_m3_chunk.sbatch
# Once all registered points complete:
Rscript inst/scripts/report_m3.R
```

For an interactive compute allocation, the same registered points can run
with `srun Rscript inst/scripts/run_m3_parallel.R`; the launcher uses the
allocated CPU count and retains already-completed points. This changes only
scheduling, with the grid and replicate seeds preserved.
