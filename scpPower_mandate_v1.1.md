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
