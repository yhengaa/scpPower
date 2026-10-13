# scpPower

[![R CMD check](https://github.com/yhengaa/scpPower/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yhengaa/scpPower/actions/workflows/R-CMD-check.yaml)

Calibration for mass-spectrometry single-cell proteomics experimental design.
The governing [mandate is frozen](scpPower_mandate_v1.1.md). This release covers
the package skeleton and `scp_calibrate()`; calibration status and limitations
are recorded in [the M1 evidence](reports/M1).

M1 has resumed under the [frozen addendum v1.1.1](scpPower_mandate_v1.1.md#frozen-addendum-v111--m1-resumption-approved-13-october-2026). All registered datasets and both Derks instruments have been evaluated. The [resumed evidence](reports/M1-resumed) contains the calibration table, plots, one-page plausibility note and manifest. Dataset-specific failures and unidentifiable components are explicitly flagged. The [previous stopped run](reports/M1) is preserved.

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

To regenerate the resumed deliverables from public inputs, run from the package root
after restoring the committed environment (or preparing the cluster library with
`bootstrap.R`):

```sh
python inst/scripts/fetch_resumed_sources.py
python inst/scripts/fetch_bulk_metadata.py
Rscript inst/scripts/inspect_resumed.R
Rscript inst/scripts/calibrate_resumed.R
python inst/scripts/extract_paper_benchmarks.py
Rscript inst/scripts/report_resumed.R
python inst/scripts/hash_hub_inputs.py M1-resumed
python inst/scripts/finalize_resumed_manifest.py
```

The resumed calibration uses the specified quantitative assays, aggregates measured
log2 features by protein median, preserves missingness, and applies cell median
shifts only for variance estimation. Specht's supplied peptide input is already
normalised and logged; its abundance calibration is flagged and stopped.
Leduc also retains a significantly negative quantification slope on its
reference-channel ratios and is stopped separately.
Derks metadata are matched by physical file and label after detecting permuted row
names. Brunner has one cell per run, leaving separate cell and run variance
unidentifiable. No dataset is replaced and no fit is adjusted to match a paper.

Compression is a registered user assumption, never an estimate from these datasets.
The manifest records the default and sensitivity assumptions and the ratio-standard
pilot requirement. Walzik supplies whole-PBMC donor variance, including donor
composition differences; sorted-subset proteomics is unavailable.

Author and maintainer: Youshan Heng,
[ORCID 0000-0002-2758-171X](https://orcid.org/0000-0002-2758-171X),
poggar@hotmail.com (verified through the connected GitHub profile).

The carrier reference identified by PMID 33288958 is Cheung et al.,
[*Nature Methods*](https://doi.org/10.1038/s41592-020-01002-5).
The frozen mandate retains the PI's text verbatim; the bibliographic correction
is recorded here and in calibration documentation.
