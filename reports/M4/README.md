# M4 evidence

[Validation report](validation-report.md), [eligible-point table](validation-table.csv),
[dataset verdicts](dataset-verdicts.csv), [manifest](run-manifest.json),
[registration](registration.json), [prediction freeze](predictions-frozen.json),
[replicate metrics](replicate-metrics.csv.gz),
[first cell assignments and native count matrices](first-assignments-and-layouts.tar.gz),
[cluster accounting](cluster-accounting.txt).

| Dataset | Eligible points | Predicted versus observed |
| --- | --- | --- |
| Specht | [Table](specht2019v3-validation-table.csv) | [Plot](specht2019v3-predicted-observed.png) |
| Leduc | [Table](leduc2022_pSCoPE-validation-table.csv) | [Plot](leduc2022_pSCoPE-predicted-observed.png) |
| Derks timsTOF | [Table](derks2022_timsTOF-validation-table.csv) | [Not evaluable](derks2022_timsTOF-predicted-observed.png) |
| Derks Q Exactive | [Table](derks2022_QExactive-validation-table.csv) | [Plot](derks2022_QExactive-predicted-observed.png) |
| Brunner | [Table](brunner2022-validation-table.csv) | [Plot](brunner2022-predicted-observed.png) |

[Reading the grid](../M3/reading-grid.pdf) is the requested interpretation of
the frozen M3 results. The original M3 numeric tables and plots are unchanged.
New presets use calibration proteins only and live in
[inst/presets/m4](../../inst/presets/m4/).

Reproduce using [the public-input workflow](../../inst/scripts/reproduce_m4.sh).
The ZIP expands the recorded assignments and layouts into directly readable CSVs.
