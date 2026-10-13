# scpPower mandate version 1.1

13 October 2026

**Status: FROZEN.** The principal investigator approved decisions 1 to 5 in the M0 closure directive of 13 October 2026. They are reproduced verbatim below. This mandate supersedes the unfrozen record. Build directive 1, section 1 onward, is authorised immediately. No further approval of these decisions is required.

The package concerns two conditions, independent biological samples, and cells within samples. The current build directive covers the repository, package skeleton, and `scp_calibrate()` only. The registered M4 criterion below is preserved; M4 execution and implementation of `scp_simulate()` or `scp_power()` remain outside this directive.

## Decision 1: hierarchical model

Keep the v1.0 model: log2 y_gi = mu_g + beta_g C_i + u_g,s(i) + v_g,b(i) + e_gi, with normal random effects of variance sigma2_s, sigma2_b, sigma2_e per protein, pooled by abundance tier (quintiles of mean log2 intensity).

- **Missingness.** Logistic regression of the probability that a protein is quantified in a cell on the protein's expected log2 intensity in that cell, fitted per dataset by maximum likelihood with protein-level random intercepts omitted (one curve per dataset). If the fitted slope's 95% interval includes zero, report that and retain the fitted intercept as an MCAR rate. Shared identification failures across a plex are treated as missing at the cell level; no imputation anywhere in calibration.
- **Carrier compression.** beta_observed = k times beta_true. k is a preset parameter, not estimated from a single-carrier dataset. Default k per carrier ratio is taken from Cheung et al. 2021, J Proteome Res, PMID 33288958; record the figure or table it comes from. Report power at the default k and at k = 1 so the user sees the cost of compression.
- **Batch.** sigma2_b is estimated only where plex or run is not aliased with sample or condition. Where aliased, report NA with the reason.

## Decision 2: endpoints

Primary endpoint: per-sample means of log2 protein intensity across the cells of each biological sample, requiring at least 3 quantified cells per sample for a protein to be testable; limma-trend with design ~ condition; Benjamini-Hochberg FDR 0.05; a protein is detected when adjusted p is below 0.05 and the estimated log2 fold change has the injected sign. Power is the fraction of truly changed proteins detected. Testability is the fraction of proteins meeting the cell minimum in every sample.

Secondary endpoint: the scplainer cell-level linear model with condition and batch as fixed terms, reported for comparison only. It is not used for the M4 gate because it does not model donors as random units.

## Decision 3: the four SCP calibration datasets

specht2019v3, leduc2022_pSCoPE, derks2022, brunner2022 (label-free DIA, scpdata). For each, calibrate on the earliest protein-level assay that has no imputation and no batch correction; name the assay used. For derks2022, use the single-cell experiment only. For leduc2022, do not use the proteins_processed assay for missingness or batch estimates.

## Decision 4: multi-donor bulk cohort for between-donor variance

Selection rule: a public MS-based proteomics dataset of a primary human cell type matching an SCP use case (PBMC, monocytes or T cells), with at least 4 healthy donors, protein-level quantities available without imputation, on PRIDE or a journal supplement. Among qualifying datasets choose the one with the most donors. Default if the search finds nothing better: Rieckmann et al. 2017, Nature Immunology, PRIDE PXD004352 (immune cell types from 4 donors), using one cell type. Report the cohort chosen, the alternatives found, the donor count, and the protein overlap with each SCP dataset. Between-donor variance is estimated per protein from donor means on log2 scale and pooled by abundance tier. Register the assumption that this variance transfers from bulk to single-cell measurements of the same cell type; it is tested only by the lab pilot, not by M4.

## Decision 5: semi-synthetic M4 registration

- **Truth construction.** For each SCP dataset, take cells from one condition only. Randomly assign cells to pseudo-samples of m cells. Assign pseudo-samples to two conditions. For 10% of proteins (random, seed recorded), add a log2 effect beta to one condition. Add a per-pseudo-sample random effect per protein drawn from the calibrated between-donor variance for that protein's tier.
- **Separation.** Split proteins into calibration and validation halves by a recorded seed. Calibrate on one half; inject and evaluate on the other.
- **Eligible grid.** All combinations of m in {1, 2, 4, 8, 16, 32, 64}, n in {3, 4, 6, 8, 12}, beta in {log2 1.25, log2 1.5, log2 2} where 2 times n times m does not exceed the dataset's cell count.
- **Observed power.** Fraction of injected proteins detected by the primary endpoint, averaged over 200 replicates per grid point, so the simulation standard error is under 3.5 percentage points. Realized FDR is reported beside it.
- **Predicted power.** scp_power() for the same m, n, beta, preset and calibration, computed before observed power is seen.
- **Decision.** Per dataset: pass if predicted is within 10 percentage points of observed on at least 80% of eligible grid points. Overall: pass if at least two datasets pass. Fail is reported and work stops; the simulator is not tuned to validation data.
- **What M4 does and does not show.** M4 validates cell-level noise, missingness and the test's behaviour. It cannot validate the between-donor variance assumption, because that variance is injected, not observed. State this in the validation report.

## Standing engineering rules

These instructions are specified directly by build directive 1 and are preserved for the first authorised package commit:

1. Create the GitHub repository on the first day of package work. `DESCRIPTION` must name the principal investigator as author and maintainer using the real name, ORCID, and email. No placeholder author is permitted.
2. Run `R CMD check` on every push in GitHub Actions using the current Bioconductor release container. Put its badge in the README. Fix a red check before the next commit.
3. Only add tests. Never delete, move out, or skip tests to obtain a green check. Test each exported function before its introducing directive closes.
4. Keep only user-facing code in `R/`. Put calibration scripts, benchmarks, and evaluation code in `inst/scripts/` or a separate repository.
5. Reproduce every reported data-derived number from public data through an `inst/scripts/` script that downloads the data, records its hash, and regenerates the number.
6. Record R, scp, limma, edgeR, and scplainer versions in every output. Warn on version mismatch; never refuse to run solely because of a mismatch.
7. Document every exported function with roxygen2. Maintain `NEWS.md` from version 0.0.1. Use `renv` and commit its lockfile.
8. Record seed, input hashes, and wall time for every run.

The current instruction authorises repository creation when the mandate precondition is satisfied. The older general restriction on a public push does not require asking again for that expressly authorised action. The M0 closure directive satisfies the mandate precondition and authorises proceeding immediately with M1.

## M1 deliverables and stop rules

Build directive 1 covers the package skeleton and `scp_calibrate()` only. Its required outputs are:

1. A repository link with green CI.
2. A table with one row per dataset and abundance tier containing proteins, cells, samples, all three variance components, missingness parameters, and compression.
3. Plots of variance components by tier, missingness versus mean intensity with the fitted curve, and compression versus carrier ratio where multiple ratios exist.
4. A one-page plausibility note comparing estimates with the source papers, flagging every disagreement without smoothing, capping, or adjustment to match.
5. A run manifest containing seeds, hashes, software versions, and wall times.

Do not implement `scp_simulate()`, `scp_power()`, or execute M4 under this directive. Do not add datasets or cohorts. Stop if the specified bulk cohort is unusable or a calibration estimate is implausible. M1 closes only when its five deliverables are returned.

## Authority and source record

The governing instructions are the build directive 1 supplied in this chat and the [M0 closure directive](</home/poggar/.codex/attachments/8379d4ea-34fd-4f5b-9153-e9eedc518a0e/Pasted text.txt>). The closure directive resolves the previously absent scientific decisions. Its approved specification takes precedence over provisional suggestions in the [earlier audit](https://chatgpt.com/s/t_6acdcb163018819182788927c5f808ee) and the original concept plan. Factual calibration blockers are reported after this freeze and do not reopen the approved decisions.

## Frozen addendum v1.1.1 — M1 resumption, approved 13 October 2026

**Status: FROZEN.** The principal investigator accepted the original mandate, repository, Walzik PXD058573 selection, prior stop, and Cheung citation correction to Nature Methods. This addendum implements the approved resumption directive. The original text above and the recorded stop are preserved. Only the clauses explicitly replaced or clarified below change. Decision 5, including the registered M4 criterion, is unchanged.

### A. Quantitative calibration inputs (replaces Decision 3 assay clause)

Use the earliest unimputed, uncorrected quantitative assay and aggregate to proteins:

| Dataset | Input | Protein aggregation |
| --- | --- | --- |
| specht2019v3 | peptides | Log2 positive quantities, then per-cell median across peptides mapped to a protein; NA only when every peptide is missing |
| leduc2022_pSCoPE | peptides | Same; document why peptides_log is or is not identical to these values logged |
| derks2022 | tims_prec_extracted and qe_prec_extracted | Single-cell columns only; log2 precursor quantities, then per-cell protein median; calibrate and report instruments separately |
| brunner2022 | proteins | Raw linear intensities; log2 only |

The missingness predictor is each protein's mean observed log2 intensity across all quantified cells in the unnormalised input. Do not centre proteins before fitting. Per-cell median shifts may be used only for variance estimation and must be named in the manifest.

### B. Culture biological units (clarifies Decision 1)

These cell-line datasets have no donors. Their biological sample is the cell line: SampleType for Specht and Leduc; Celltype for Derks; HeLa for Brunner. There is one culture sample per line, not an independent donor. The table's samples count is the number of cell lines. SCP sigma2_s remains unavailable; the mapped donor component comes from Walzik.

Sigma2_e is residual cell-level log2 variance around the cell line's mean within its batch. Sigma2_b is variance of batch means after removing the cell line mean where batch contains cells of the same line. Use Set for Specht/Leduc and run for Derks/Brunner. Where a batch holds only one cell line, batch and line are aliased for that line: report NA with that reason per line. Identification failures caused by insufficient replication must also be stated explicitly.

### C. Dataset-specific stop (replaces implausibility run stop)

An implausible estimate stops that dataset's row, not other datasets. Calibrate every registered dataset; flag every implausible estimate in the plausibility note without smoothing, capping, adjusting, or changing the model. Stop the entire run only when no dataset calibrates. The unusable-bulk-cohort stop remains.

### D. Compression assumption (replaces Decision 1 compression clause)

Cheung et al. 2021, Nature Methods, PMID 33288958 reports instrument- and ion-sampling-dependent carrier effects; no published per-ratio compression table is retrievable. k is a registered user assumption, not a calibrated quantity: default k=1 (no compression). Future power reporting must include the sensitivity assumptions k in {0.75, 0.5}. This directive does not implement that power function. Record that k was not estimated and recommend a two-proteome ratio standard in the lab pilot. The compression plot is dropped from M1.

### Walzik clarification and resumed deliverables

Keep PXD058573. Establish whether it contains proteomics of sorted immune subsets in addition to whole PBMC. If present, report whole-PBMC donor variance and variance of the subset closest to each SCP line's lineage: monocyte-like for U-937 and THP-1 style lines; lymphoid for other lines. State that whole-PBMC variance includes between-donor composition differences. Do not substitute another cohort or invent sorted subsets.

Return together: (1) repository with green CI at the resumed commit; (2) calibration table for every dataset and abundance tier, with both Derks instruments separate; (3) variance-by-tier and fitted missingness panels; (4) a one-page note comparing missing fraction, proteins per cell, and cell CV with source papers, flagging disagreements without adjustment; (5) seeds, input hashes, versions, and wall-time manifest. Preserve the earlier stop evidence. No scp_simulate(), scp_power(), or M4 work is authorised.

## Frozen addendum v1.1.2 — M1 closure and M2, approved 13 October 2026

**Status: FROZEN.** The principal investigator accepts the resumed M1 evidence and approves amendments E and F and the M2 gate below. This addendum replaces only the specified TMT missingness and one-cell-per-run variance clauses. M4 registration in Decision 5 is unchanged. The approved directive is reproduced verbatim below.

# Directive: M1 closure and M2 simulator

13 October 2026

The resumed M1 bundle is accepted. Derks timsTOF, Derks Q Exactive and Brunner are calibrated. The CI run, the 77 tests and the plausibility comparisons are accepted as delivered. Two amendments close the remaining rows; they are decisions, approved now. Add them to the frozen mandate as addendum v1.1.2, record the diff in NEWS.md, then proceed to M2 in the same directive.

## Amendment E: abundance predictor for TMT datasets (Specht, Leduc)

The scpdata peptide assays for both TMT datasets are reference-channel ratios or row-normalised values, so they carry no abundance scale and the missingness curve cannot be fitted from them. The raw scale exists in the per-run PSM assays (the 177 Specht and 134 Leduc run assays), whose columns are the raw reporter-ion intensities of each channel.

- Missingness indicator per protein per single-cell channel: quantified if the raw reporter-ion intensity is finite and above zero in at least one PSM of that protein in that run; otherwise missing. Carrier, reference, empty and blank channels are excluded by `SampleType`.
- Abundance predictor per protein: the mean over single-cell channels of the log2 raw reporter-ion intensity, aggregated from PSMs to protein by median, computed only where quantified. No normalisation before this step.
- Fit the registered logistic curve on these two quantities. Report the curve and the slope interval as for the other datasets.
- Variance components for Specht and Leduc stay as delivered (from the ratio-based peptide assays, which is the scale labs test on). Mark both datasets CALIBRATED once the curve is fitted, with a note that the missingness predictor and the variance scale come from different assays.

## Amendment F: batch term in one-cell-per-run data (Brunner)

When every cell is its own acquisition run, batch and cell noise are the same unit. For such datasets report sigma2_e as the total within-line variance of a cell's log2 value after the per-cell median shift, set sigma2_b to 0 with the reason "absorbed into sigma2_e: one cell per run", and record that this is the per-cell variance a label-free design will actually see. This is not an imputation; it is the identifiable quantity.

## M1 closure

With E and F applied, M1 closes with five calibrated rows (Specht, Leduc, Derks timsTOF, Derks Q Exactive, Brunner). Regenerate the table, the two plots, the plausibility note and the manifest; tag the repository v0.1.0 and publish a GitHub release with the M1 evidence attached.

## M2: scp_simulate()

Implement `scp_simulate(calibration, n_per_condition, m_cells_per_sample, beta, preset, layout, k = 1, seed)` returning a protein-by-cell log2 matrix with NA for missing values and a cell metadata table (sample, condition, batch).

Generation, per protein g in tier t:
1. mu_g drawn from the calibrated mean log2 of tier t.
2. Sample effect u_gs ~ N(0, sigma2_s,t) from the Walzik tier estimate, one per biological sample.
3. Condition effect k * beta_g added to one condition for the proteins flagged as changed (10% by default; argument).
4. Batch effect v_gb ~ N(0, sigma2_b,t), one per batch, with batches assigned by `layout`: `balanced` (each plex or run holds equal cells from both conditions), `confounded` (each batch holds one condition), or a user matrix. Label-free presets use one batch per cell and sigma2_b = 0 by amendment F.
5. Cell residual e_gi ~ N(0, sigma2_e,t).
6. Missingness: each value is set to NA with probability from the dataset's logistic curve evaluated at the expected log2 intensity mu_g + u_gs + v_gb + condition effect (not the realised value).

Presets are YAML files in `inst/presets/`, one per calibrated dataset row, holding the tier table, the curve, the carrier ratio and k. A user can supply their own.

## M2 gate, registered now

For each calibrated dataset, simulate the dataset's own design (its cell lines as samples, its real batch sizes, beta = 0) with the same cell count, 20 replicates, and compare with the observed data:

| Quantity | Tolerance |
| --- | --- |
| Missing fraction per tier | within 5 percentage points |
| Mean proteins per cell | within 10% relative |
| Median within-line cell SD per tier (log2) | within 20% relative |
| Between-batch SD per tier, where identified | within 25% relative |

A dataset passes when all four hold. Report every dataset, pass or fail, with the numbers. Do not change the calibration or the simulator to pass; a failure is reported with the suspected cause and the next directive decides.

## Deliverables

1. Repository link at the tagged v0.1.0 release and at the M2 commit, CI green.
2. Final M1 table, plots, plausibility note and manifest.
3. `scp_simulate()` documented and tested (at least: reproducibility under a seed, correct dimensions, missing fraction responds to the curve, batch effect responds to layout).
4. M2 fidelity report: the table above per dataset, with plots of simulated versus observed missingness by intensity and cell SD by tier.
5. Run manifest with seeds, hashes, versions and wall times.

No work on scp_power() or M4. Stop rules remain per dataset.

## Frozen addendum v1.1.3 — M2 revision and conditional M3, approved 13 October 2026

**Status: FROZEN.** Amendment G changes simulator structure only. The registered M2 tolerances, replicate count and metric definitions and the M4 criterion remain unchanged. M3 and the v0.2.0 tag require all five series to pass. The approved directive is reproduced verbatim below.

# Directive: M2 revision and M3

13 October 2026

The v0.1.0 release, the final M1 calibration and the M2 fidelity report are accepted as delivered. Four datasets failed the registered gate and the report says why without tuning; that is what M2 is for. The failure pattern points to two structural errors in the simulator, not to the calibration or the tolerances. Amendment G below fixes them. It is a decision, approved now; add it to the mandate as addendum v1.1.3 and record it in NEWS.md. The M2 tolerances, replicate count and metric definitions do not change. The M4 criterion does not change.

## Why the gate failed

1. **Cell and batch SD too high everywhere.** The simulator gives every protein in a tier the tier's mean variance. The gate compares the median per-protein SD. The mean of a right-skewed variance distribution is well above the median, so every simulated SD is inflated. Brunner passed only because its variances are less skewed.
2. **TMT missingness flat across tiers.** In Specht the fitted intensity slope is about zero, so one curve predicts the global missing fraction for every tier while the observed tiers range from 47% to 85%. The per-protein missingness plot is bimodal: a protein is quantified in most cells or in almost none. That is the signature of plex-level identification driven by the carrier, with protein-specific detectability, not of cell-level intensity. One global curve cannot represent it.

## Amendment G: simulator structure

**G1. Per-protein parameters by resampling.** A simulated protein is a resampled real protein. Draw, jointly and with replacement, one calibrated protein from the preset's tier and take its mean log2 intensity, its sigma2_e, its sigma2_b (0 where unidentified, with the tier value substituted only when the protein's own estimate is NA), its mapped sigma2_s (tier value when unmapped), and its detection parameters below. Tier summaries remain in the table for reporting; they are no longer the generating parameters. The preset YAML therefore carries the per-protein parameter table (one row per calibrated protein), not only tier rows.

**G2. Two-level missingness for TMT presets.** For each protein g: p_g = fraction of plexes (runs) in which g has at least one PSM; q_g = fraction of single-cell channels in which g is quantified, among plexes where it is identified. Both come from the raw PSM mask already extracted for M1. In simulation, a protein is identified in a plex with probability p_g; if identified, each cell in that plex is quantified with probability q_g; otherwise the whole plex is NA for that protein. The condition effect does not change p_g or q_g, because the carrier that drives identification pools both conditions. Record this as a registered modelling consequence: in carrier-TMT designs, testability depends on the number of plexes and barely on the effect size.

**G3. Per-protein detection for DIA and label-free presets.** Each protein keeps its own observed quantification rate r_g as the baseline. The registered intensity curve supplies only the slope: the logit of detection for a cell is logit(r_g) + slope times (expected log2 intensity minus the protein's calibrated mean). The condition effect therefore shifts detection in the changed condition, which is the correct behaviour for DIA.

**G4. Everything else unchanged.** Sample, batch and residual draws, the layout argument, k, seeds, hashes and manifests stay as implemented.

## Gate rerun

Rerun the registered M2 gate exactly as before: native designs, beta = 0, 20 replicates, the same four tolerances and metric definitions. Add one reported quantity, not a gate: the fraction of proteins quantified in at least half of the cells, observed versus simulated. Report every dataset pass or fail. If any dataset fails, stop after the report with the suspected cause.

## M3, only if all five pass

Implement `scp_power(preset, n_per_condition, m_cells_per_sample, beta, layout, k, replicates = 200, seed)` over the registered grid: m in {1, 2, 4, 8, 16, 32, 64}, n in {3, 4, 6, 8, 12}, beta in {log2 1.25, log2 1.5, log2 2}, for every preset. For each grid point it simulates, applies the primary endpoint (per-sample means over cells with at least 3 quantified cells, limma-trend, BH 0.05, sign agreement), and returns power, testability, realised FDR, and the simulation standard error.

Instrument time: each preset YAML gains `cells_per_run` (TMT: single-cell channels per plex; plexDIA: 3; label-free: 1) and `minutes_per_run` taken from the source paper's gradient, with the citation in the YAML. `scp_power()` reports runs and instrument hours per design and the cheapest design in hours reaching a target power (default 0.8) for each beta.

Run the grid on the cluster. Deliver power curves (power against m for each n, one panel per beta and preset), the testability curves, the realised FDR table, run time per grid point, and the manifest. Report power at k = 1 and at k = 0.75 and 0.5 for TMT presets.

## Deliverables

1. Repository link, CI green, tagged v0.2.0 after the gate passes.
2. Revised M2 fidelity report with the rerun table and the two plots.
3. If authorised: `scp_power()` documented and tested, the M3 power and testability curves, the FDR table, and instrument-hour tables.
4. Run manifest.

No M4 work. Stop rules remain per dataset.

## Frozen addendum v1.1.4 — M4 validation, approved 13 October 2026

**Status: FROZEN.** The principal investigator accepts v0.2.0, the five M2
passes and the completed M3 grid. Decision 5 is now authorised. Its criterion
is unchanged: 200 replicates per grid point; absolute predicted–observed power
gap at most 10 percentage points on at least 80% of counted points per dataset;
at least two datasets must pass overall. No tuning after validation outcomes.

**H1. Source cells.** Use the largest single cell line: Specht macrophage,
Leduc monocyte, Brunner HeLa, and the largest Derks Q Exactive line. Derks
timsTOF has ten cells per line and no eligible grid point with three cells per
sample; report it as not evaluable and do not run degenerate points.

**H2. Pass fraction.** List eligible m < 3 points with their constructional
zero predicted and observed power, but exclude them from the 80% pass fraction.
Only eligible points with m >= 3 count. This tightens the registered gate.

**H3. Prediction layout.** Record the pseudo-sample-by-native-batch count
matrix from the first replicate's source-cell assignment for every point.
Pass that matrix as the custom layout to scp_power(). For Brunner there is one
cell per batch. Finish predicted power before observing observed power.

**H4. Donor effect and halves.** Split proteins into calibration and validation
halves with a recorded seed. Fit the half-preset from calibration proteins only
and store it under inst/presets/m4/. Released presets stay untouched. Inject
mapped per-protein sigma2_s under G1 into validation pseudo-samples, using the
tier value where unmapped. Between-donor variance is injected, not observed.

**R1 and R2. Reporting only.** Add conditional power (detected changed proteins
divided by testable changed proteins) beside registered unconditional power.
Add mean rejections and the fraction of replicates with no rejections beside
realised FDR. Existing definitions and M3 numerical results do not change.

Report eligible points, predicted and observed power, gap, within-10-pp flag,
observed testability and FDR, simulation SE and one comparison plot per dataset.
Report each verdict and the overall verdict; on failure stop after the report
with suspected causes and no tag. The validation report's first paragraph
states that M4 validates cell-level noise, missingness and the test and cannot
validate bulk donor variance transfer or k. Add a one-page Reading the grid
section to the M3 report on testability, 1.25-fold changes and k sensitivity.
Record split and replicate seeds, hashes, versions, wall times and cluster
accounting. On overall pass, tag v0.3.0 and create a GitHub release with M4
evidence. No M5 vignette, README design examples or Zenodo work is authorised.
