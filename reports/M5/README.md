# scpPower 1.0.0 release evidence

The build closes with the [M3b report](../M3b/power-report.md), [combined grid](../M3b/power-grid.csv), refreshed plots, [run manifest](../M3b/run-manifest.json), and [cheapest-design difference](../M3b/cheapest-design-diff.csv). The old TMT rows are archived alongside them. No scientific gate or registered quantity changed.

The [README](../../README.md), [vignette source](../../vignettes/designing-scp.Rmd), [pilot protocol](../../inst/doc/pilot-protocol.md), and [CITATION.cff](../../CITATION.cff) form the user documentation. The vignette executes its simulation, power and compression calls when the package is built and is rebuilt by R CMD check.

[Worked examples](worked-examples.csv) have their own [manifest](worked-examples-manifest.json), distinct from the registered grid. [Release reproduction](../../inst/scripts/reproduce_release.sh) downloads the public inputs and records hashes through the existing acquisition scripts. Earlier stages have their registered scripts and immutable reports. The refreshed grid preserves the 315 non-TMT results without rerunning them.

The release workflow records the exact green commit in CI-receipt.json inside the unified evidence attachment before publication. The [release notes](release-notes.md) provide the three author steps for Zenodo. No DOI is available until the author enables that integration and Zenodo processes a release; the requested DOI badge and citation update therefore await v1.0.1.
