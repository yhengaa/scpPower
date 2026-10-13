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
