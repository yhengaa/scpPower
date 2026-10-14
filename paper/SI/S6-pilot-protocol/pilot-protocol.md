# Lab pilot: calibrate your own cell type

Use **4 donors × 2 conditions × 20 cells per donor and condition** as a practical default, not a derived minimum or a powered design. Record donor, condition, acquisition batch, single-cell channel and blank status. Distribute every donor and both conditions across mixed acquisition batches; include replicate cells from each sample within batches. Use a unique donor–condition sample ID. The independent-sample endpoint does not analyse donor pairing; the pilot can measure donor variance without establishing power for a paired analysis.

Measure three quantities. First, estimate per-protein between-donor variance in your cell type from log2 donor means in a single reference condition, using at least four quantified donors for each protein. Supply this protein-by-donor matrix as `donor_means`; `donor_info$donor` identifies its columns. This uses the registered unbiased variance estimator and exact protein-ID mapping. Sampling noise in pilot donor means remains part of that estimate; record its uncertainty and overlap rather than treating this small pilot as precise. The SCP matrix supplies residual cell and identifiable batch variance. Second, in carrier TMT include a two-proteome standard with known ratios, processed through the intended carrier ratio, instrument and endpoint pipeline. Estimate the log2 ratio attenuation externally and record it as `pilot_k`, with its uncertainty; scp_calibrate() records this supplied assumption, rather than estimating compression itself. Third, include blank wells in each acquisition batch to measure background and contamination. Blanks are controls, excluded from biological cells and from the calibration matrix; report their signal separately without imposing a new background correction.

Use an unimputed protein assay on the endpoint scale, before protein centring or batch correction. For TMT, retain the endpoint protein matrix after the same peptide/cell filters and aggregation that the test will see; use raw abundance measurements through `missingness_input` if the endpoint contains ratios. Median-shift cell normalisation applies only to variance estimation. In the call below, `pilot` is QFeatures, `proteins` is an already-log2 protein assay of biological cells, its colData contains `sample`, `condition` and `plex`, `raw_proteins` is the aligned unimputed log2 abundance matrix, and bulk identifiers exactly match assay row names.

```r
calibration <- scpPower::scp_calibrate(
  pilot, assay = "proteins", log2_input = TRUE,
  sample_col = "sample", condition_col = "condition", batch_col = "plex",
  biological_unit = "sample", normalize_cells = TRUE,
  bulk = donor_means, bulk_metadata = donor_info, bulk_log2_input = TRUE,
  missingness_input = raw_proteins, missingness_log2_input = TRUE,
  carrier_ratio = pilot_carrier_ratio, compression_k = pilot_k, seed = 20261013)
source(system.file("scripts", "pilot_to_preset.R", package = "scpPower"))
preset <- pilot_to_preset(calibration,
  endpoint = SummarizedExperiment::assay(pilot[["proteins"]]),
  metadata = as.data.frame(SummarizedExperiment::colData(pilot[["proteins"]])),
  batch_col = "plex", cells_per_run = pilot_single_cell_slots,
  minutes_per_run = pilot_run_minutes, carrier = TRUE)
```

Review `calibration$proteins`, `$tiers`, component reason columns and `$manifest` before using the preset. Unavailable donor or residual tier variances must be reported and prevent meaningful simulation; no value is invented. The conversion utility uses the existing G1 fallbacks and amendment-I endpoint mask, without changing the model. Use `scp_simulate(NULL, 4, 20, beta = log2(2), preset = preset)` and `scp_power(preset, 4, 20, beta = log2(2))` with recorded seeds to explore the pilot assumptions. For label-free runs, supply `carrier_ratio = NULL`, `batch_col = "run"`, `absorb_cell_batches = TRUE`, use one cell per run, and set `carrier = FALSE`; separate cell and run noise cannot be identified. Record hashes, software versions, run time, standard ratios and blank measurements alongside the calibration manifest.
