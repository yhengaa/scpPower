# scpPower

[![R CMD check](https://github.com/yhengaa/scpPower/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yhengaa/scpPower/actions/workflows/R-CMD-check.yaml)

Calibration for mass-spectrometry single-cell proteomics experimental design.
The governing [mandate is frozen](scpPower_mandate_v1.1.md). This release covers
the package skeleton and `scp_calibrate()`; calibration status and limitations
are recorded in [the M1 evidence](reports/M1).

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

Author and maintainer: Youshan Heng,
[ORCID 0000-0002-2758-171X](https://orcid.org/0000-0002-2758-171X),
poggar@hotmail.com (verified through the connected GitHub profile).

The carrier reference identified by PMID 33288958 is Cheung et al.,
[*Nature Methods*](https://doi.org/10.1038/s41592-020-01002-5).
The frozen mandate retains the PI's text verbatim; the bibliographic correction
is recorded here and in calibration documentation.
