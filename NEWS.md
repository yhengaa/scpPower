# scpPower 0.3.1

* Both M2b rows pass. M4b passes on all counted points: Specht 72/72 and
  Leduc 69/69. In the descriptive power >= 0.2 window, agreement is 4/4
  and 10/10; mean signed predicted-minus-observed gaps are -0.80 and +1.45 pp.
  Preserves the corrected startup error before any M4b prediction or observation.

* Frozen addendum v1.1.5 records Diagnostic J before any preset change.
  Endpoint matrices are sparser than raw PSM masks by 9.781 pp for Specht
  macrophages and 7.080 pp for Leduc monocytes in the v0.3.0 validation halves.
  Both exceed the fixed 5-pp threshold, authorising amendment I.
* Amendment I uses the endpoint aggregated assay mask for TMT p_g/q_g and
  the M2b missingness comparison. Intensity, variance, simulation draws and
  endpoint definitions are unchanged. M2b and M4b retain registered gates.
* Adds limitation K to the validation report and README; removes the old
  TMT power-grid example pending the diagnostic round. Reports actual recorded
  conditional-power gaps rather than asserting the unreproduced 1-pp bound.
  The accepted v0.3.0 M4 verdict and all non-TMT results stand. No M5 work.

# scpPower 0.3.0

* Registered M4 passes overall: Specht 67/72, Leduc 66/69, Q Exactive 6/6
  and Brunner 57/57 counted points within 10 percentage points. TimsTOF is
  not evaluable. Retains all eight discrepant points without adjustment.
  M4 does not validate bulk donor variance transfer or compression k.

* Frozen addendum v1.1.4 records H1–H4: single-line validation source cells,
  m >= 3 in the pass fraction, prediction using the first assignment's native
  batch layout, mapped per-protein donor injection and separate half-presets.
  The 200 replicates, 10 percentage point tolerance, 80% dataset pass fraction
  and requirement for two passing datasets remain unchanged.
* Adds conditional power beside unconditional power, mean rejections and the
  zero-rejection replicate fraction beside realised FDR. Existing endpoint,
  simulation and M3 results are unchanged; the added columns have tests.
* Adds the requested interpretation of the frozen M3 grid to its report.
  M4 cannot validate donor variance transfer or compression. No M5 work.

# scpPower 0.2.0

* All five calibrated series pass the unchanged M2 fidelity gate under
  amendment G, authorising M3. Preserves both earlier M2 reports.
* Introduces documented and tested `scp_power()` for the registered per-sample
  means endpoint, limma-trend, BH 0.05 and injected-sign agreement. Reports
  power, testability, realised FDR, Monte Carlo standard errors and actual
  layout-dependent run counts and cited instrument-hour estimates.
* M3 evaluates the registered grid and TMT compression sensitivity assumptions.
  No M4 evaluation is performed and its criterion remains unchanged.

# scpPower 0.1.2

* Frozen addendum v1.1.3 implements approved amendment G: jointly resample
  per-protein intensity and variance parameters within tiers; model TMT plex
  identification followed by conditional channel detection; preserve each DIA
  protein's baseline detection rate with the registered intensity slope.
* Retains the M2 gate definitions, tolerances, seeds and replicate count.
  Adds the fraction quantified in at least half the cells as a descriptive
  quantity outside the gate. Preserves the original M2 report.
* M3 and tagging v0.2.0 are conditional on all five series passing the rerun.
  The registered M4 criterion is unchanged.

# scpPower 0.1.1

* Introduces documented and tested `scp_simulate()` with biological-sample,
  batch and cell effects, expected-intensity missingness, seeded truth, and
  balanced, confounded or custom sample-by-batch layouts.
* Provides a YAML preset for every final M1 dataset/instrument series and
  the registered M2 fidelity report, preserving all failures without tuning.

# scpPower 0.1.0

* Frozen addendum v1.1.2: TMT missingness is fitted from raw reporter-ion
  PSM assays while preserving the accepted ratio-based variance components.
  Brunner's one-cell-per-run variance is reported as total within-line cell
  variance after cell median shifts, with batch variance absorbed into it.
* Closes M1 with five dataset/instrument calibrations and registers the
  M2 simulator and its fidelity gate, without changing the M4 criterion.

# scpPower 0.0.1

* Introduces `scp_calibrate()` for processed matrices and explicitly selected
  protein assays in QFeatures objects.
* Reports unavailable or aliased components with reasons, fits an unimputed
  quantification curve, and maps variance from independent bulk donors.
* Includes a frozen mandate, reproducible calibration scripts, run manifests,
  and checks in the Bioconductor release container.
* Simulation, power calculations, and M4 evaluation remain outside this release.

* Frozen addendum v1.1.1 (13 October 2026): replace processed protein assays
  with the specified raw peptide/precursor inputs and protein medians; use
  unnormalised protein means for missingness; register cell lines as culture
  samples for residual and batch components; apply implausibility stops per
  dataset; replace empirical compression presets with the default k=1
  assumption and future sensitivity assumptions 0.75 and 0.5; drop the
  compression plot. Preserve Walzik and verify whether sorted subsets exist.
  The registered M4 criterion and earlier stopped-run evidence are unchanged.
