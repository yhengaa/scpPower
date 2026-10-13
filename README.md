# scpPower

[![R CMD check](https://github.com/yhengaa/scpPower/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yhengaa/scpPower/actions/workflows/R-CMD-check.yaml)

Calibration for mass-spectrometry single-cell proteomics experimental design.
The governing [mandate is frozen](scpPower_mandate_v1.1.md). This release covers
the package skeleton and `scp_calibrate()`; calibration status and limitations
are recorded in [the M1 evidence](reports/M1).

**M1 is stopped, not closed.** The Leduc fit has an implausible negative
quantification slope in the normalized protein matrix. Its result is preserved;
Derks and Brunner were not calibrated after the stop. The source matrices and
unavailable variance and compression components are explained in the
[plausibility note](reports/M1/plausibility-note.md).

`scp_calibrate()` accepts a protein-by-cell matrix with cell metadata or an
explicit protein assay from a `QFeatures` object. It returns variance components,
abundance quintiles, a fitted quantification curve, externally sourced carrier
presets, mapped bulk donor variance, and a run manifest. Unavailable quantities
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

Biological sample IDs must describe genuine donors or cultures. Acquisition
runs, glass slides, and cell types are not substituted for those units. One
cell per run cannot separate acquisition variance from cell noise. Bulk donor
variance transfers only under the registered cell-type assumption, which
requires a lab pilot and cannot be validated by the semi-synthetic M4 gate.

Every calibration output includes R, scp, limma, edgeR and the scplainer workflow
version. Scplainer is implemented in scp; it is not treated as a separate R
package. Version differences warn and do not block execution.

The scripts in [inst/scripts](inst/scripts) download public inputs, record
hashes, preserve seeds and wall time, and regenerate reported results. Use
`bootstrap.R` to prepare the R library on the cluster and `inspect_scp.R` to
fetch and inspect the four registered datasets. Raw data and local libraries
are excluded from Git. The committed `renv.lock` records the environment.

To reproduce the preserved run from public inputs, run these from the package root:

```sh
python inst/scripts/fetch_sources.py
python inst/scripts/fetch_bulk_metadata.py
RENV_CONFIG_AUTOLOADER_ENABLED=FALSE Rscript inst/scripts/inspect_scp.R
RENV_CONFIG_AUTOLOADER_ENABLED=FALSE Rscript inst/scripts/calibrate_m1.R
RENV_CONFIG_AUTOLOADER_ENABLED=FALSE Rscript inst/scripts/report_m1.R
RENV_CONFIG_AUTOLOADER_ENABLED=FALSE Rscript inst/scripts/render_note_pdf.R
python inst/scripts/hash_hub_inputs.py
```

The renderer reports the recorded stop and does not fit further models.

Author and maintainer: Youshan Heng,
[ORCID 0000-0002-2758-171X](https://orcid.org/0000-0002-2758-171X),
poggar@hotmail.com (verified through the connected GitHub profile).

The carrier reference identified by PMID 33288958 is Cheung et al.,
[*Nature Methods*](https://doi.org/10.1038/s41592-020-01002-5).
The frozen mandate retains the PI's text verbatim; the bibliographic correction
is recorded here and in calibration documentation.
