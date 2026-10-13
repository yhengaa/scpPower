# M1 evidence from a stopped calibration run

The first package CI run is green. The scientific calibration is stopped at Leduc; M1 is not closed.

- [Calibration table](calibration-table.csv): all registered datasets and quintiles; unrun or unidentifiable entries remain NA.
- [Variance components](variance-components.png)
- [Missingness fit](missingness-fit.png)
- [Compression status](compression-status.png): no numerical compression series available.
- [Plausibility note](plausibility-note.md)
- [Run manifest](run-manifest.json)
- [Bulk cohort and alternatives](bulk-cohort-selection.json)
- [Dataset versions, assays, and hashes](dataset-inspection.json)

No simulations, power calculations, or M4 evaluation were performed.

Software: R=4.6.1; scp=1.22.0; limma=3.68.5; edgeR=4.10.5; scplainer=workflow in scp 1.22.0; no separate package
