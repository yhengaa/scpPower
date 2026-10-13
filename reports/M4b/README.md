# M4b evidence — amendment I

[Validation report](validation-report.md), [eligible-point table](validation-table.csv),
[dataset verdicts](dataset-verdicts.csv), [descriptive power window](high-power-window.csv),
[run manifest](run-manifest.json), [registration](registration.json),
[prediction freeze](predictions-frozen.json), [replicate metrics](replicate-metrics.csv.gz),
[first assignments and native count matrices](first-assignments-and-layouts.tar.gz),
[cluster accounting](cluster-accounting.txt).

| Dataset | Table | Predicted versus observed |
| --- | --- | --- |
| Specht | [Table](specht2019v3-validation-table.csv) | [Plot](specht2019v3-predicted-observed.png) |
| Leduc | [Table](leduc2022_pSCoPE-validation-table.csv) | [Plot](leduc2022_pSCoPE-predicted-observed.png) |

[Diagnostic J](../J/diagnostic-report.md) was reported before changing masks;
[M2b](../M2b/fidelity-report.md) passed before M4b. The failed startup and its
correction before any prediction or observation are retained under startup-aborted/.
The accepted original M4 numerical results and non-TMT results are unchanged.
The correction to K distinguishes dataset-average conditional agreement from
per-point bounds; the original verdict and criterion stand.

[Public-input reproduction](../../inst/scripts/reproduce_tmt_round.sh).
The evidence ZIP expands recorded cell assignments and native count matrices as CSVs.
No M5 work.
