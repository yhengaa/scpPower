# Contributing to scpPower

Report bugs and proposed changes through [GitHub issues](https://github.com/yhengaa/scpPower/issues). Include a small reproducible example, the dataset or preset, the seed and your software versions. Do not post private data or credentials.

Open a pull request against `main`. Explain the problem, the resulting behaviour and how you checked it. User-facing functions belong in `R/`; calibration, benchmark and evaluation scripts belong in `inst/scripts/`. Document exported functions with roxygen2, add meaningful tests and record changes in `NEWS.md`. Existing tests must remain in place and run without skips introduced to hide failures.

Use the committed `renv.lock` for the registered environment. Run the package tests and `R CMD check` before submitting. The Bioconductor release-container check must pass before merging. Version differences should warn and continue.

The frozen mandate, presets and reports are scientific records. Do not overwrite a registered result or change a gate to make it pass. Propose a separate, explicit amendment and retain the original evidence. A new reported numerical result needs a public-data reproduction script, seed, input hashes, software versions and wall time. Manuscript corrections belong in the manuscript; audit discrepancies never authorise changes to registered values.

Contributions are distributed under the package’s GNU GPL version 3 licence. Follow [the Code of Conduct](CODE_OF_CONDUCT.md).
