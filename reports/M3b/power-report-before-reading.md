# M3b updated carrier-TMT simulation power report
13 October 2026. Registered M3 grid v1.1.3, with amendment-I masks under v1.1.5. All five M2 series passed before the original M3. Only the two TMT presets use updated endpoint masks; variances, parameter draws, layouts, k assumptions and registered seeds remain fixed. This is simulation-based planning; the accepted M4 and M4b verdicts are unchanged.

The combined table contains every registered grid point: 945 points with 200 independent replicates each. TMT uses k=1, 0.75 and 0.5; DIA and label-free use k=1.
Maximum power simulation standard error: 0.63 percentage points. Grid-point wall times range from 1.9 to 193.5 seconds (sum 7.52 worker wall hours; excludes scheduler waiting and source preparation).

## Endpoint and interpretation

Sample means require at least three quantified cells in every sample. Limma-trend fits design ~ condition; BH uses 0.05 among testable proteins. Changed proteins require sign agreement for detection. Power has all truly changed proteins in the denominator, including untestable proteins. Testability counts all proteins satisfying the cell minimum. FDR is mean replicate FDP among BH rejections, null rejections divided by all rejections, and zero with no rejections; wrong-sign changed rejections remain separately recorded. Standard errors are replicate SD divided by sqrt(replicates).

Designs with fewer than three cells per sample cannot meet the endpoint and have zero testability and power. Every registered point is still simulated; no points are removed. Increasing cells can reduce cell noise and improve testability, but biological-sample variance remains. In TMT, identification/detection does not respond to effect size because the pooled carrier governs p_g and q_g.

## Realised FDR

| Preset | k | Maximum realised FDR over grid |
| --- | ---: | ---: |
| specht2019v3 | 1 | 0.0725 |
| specht2019v3 | 0.75 | 0.0620 |
| specht2019v3 | 0.5 | 0.0500 |
| leduc2022_pSCoPE | 1 | 0.0542 |
| leduc2022_pSCoPE | 0.75 | 0.0550 |
| leduc2022_pSCoPE | 0.5 | 0.0558 |
| derks2022_timsTOF | 1 | 0.0667 |
| derks2022_QExactive | 1 | 0.0675 |
| brunner2022 | 1 | 0.0575 |

Some simulated designs exceed nominal BH 0.05 in realised FDR. The full table shows the extent by design; no variance, detection, test, or threshold was adjusted. Monte Carlo standard errors are reported alongside estimates; these are estimates under the registered model, not guarantees of experimental FDR.

## Cheapest tested design reaching target power

| Preset | Fold change | k | n per condition | m cells per sample | Power | Instrument hours | Status |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| specht2019v3 | 1.25 | 1 | NA | NA | NA | NA | No registered design reaches target |
| specht2019v3 | 1.25 | 0.75 | NA | NA | NA | NA | No registered design reaches target |
| specht2019v3 | 1.25 | 0.5 | NA | NA | NA | NA | No registered design reaches target |
| specht2019v3 | 1.50 | 1 | NA | NA | NA | NA | No registered design reaches target |
| specht2019v3 | 1.50 | 0.75 | NA | NA | NA | NA | No registered design reaches target |
| specht2019v3 | 1.50 | 0.5 | NA | NA | NA | NA | No registered design reaches target |
| specht2019v3 | 2.00 | 1 | NA | NA | NA | NA | No registered design reaches target |
| specht2019v3 | 2.00 | 0.75 | NA | NA | NA | NA | No registered design reaches target |
| specht2019v3 | 2.00 | 0.5 | NA | NA | NA | NA | No registered design reaches target |
| leduc2022_pSCoPE | 1.25 | 1 | NA | NA | NA | NA | No registered design reaches target |
| leduc2022_pSCoPE | 1.25 | 0.75 | NA | NA | NA | NA | No registered design reaches target |
| leduc2022_pSCoPE | 1.25 | 0.5 | NA | NA | NA | NA | No registered design reaches target |
| leduc2022_pSCoPE | 1.50 | 1 | NA | NA | NA | NA | No registered design reaches target |
| leduc2022_pSCoPE | 1.50 | 0.75 | NA | NA | NA | NA | No registered design reaches target |
| leduc2022_pSCoPE | 1.50 | 0.5 | NA | NA | NA | NA | No registered design reaches target |
| leduc2022_pSCoPE | 2.00 | 1 | NA | NA | NA | NA | No registered design reaches target |
| leduc2022_pSCoPE | 2.00 | 0.75 | NA | NA | NA | NA | No registered design reaches target |
| leduc2022_pSCoPE | 2.00 | 0.5 | NA | NA | NA | NA | No registered design reaches target |
| derks2022_timsTOF | 1.25 | 1 | NA | NA | NA | NA | No registered design reaches target |
| derks2022_timsTOF | 1.50 | 1 | NA | NA | NA | NA | No registered design reaches target |
| derks2022_timsTOF | 2.00 | 1 | 12 | 32 | 0.865 | 192.00 | Target reached in registered grid |
| derks2022_QExactive | 1.25 | 1 | NA | NA | NA | NA | No registered design reaches target |
| derks2022_QExactive | 1.50 | 1 | NA | NA | NA | NA | No registered design reaches target |
| derks2022_QExactive | 2.00 | 1 | NA | NA | NA | NA | No registered design reaches target |
| brunner2022 | 1.25 | 1 | NA | NA | NA | NA | No registered design reaches target |
| brunner2022 | 1.50 | 1 | NA | NA | NA | NA | No registered design reaches target |
| brunner2022 | 2.00 | 1 | NA | NA | NA | NA | No registered design reaches target |

The target is mean simulated power 0.8 within this finite grid; it is not a lower confidence bound. An unavailable target is reported explicitly rather than extrapolated.

## Instrument time

Runs count actual occupied batches in the balanced simulation layout. A three-channel plexDIA run contains two cells under the registered equal-A/B rule, leaving one slot empty; hours therefore use actual runs rather than dividing cells by nominal channel capacity. Preset durations are cited in YAML and parsed from the same registered source papers. Specht and Leduc use total published TMT method durations; Derks uses the full single-cell instrument-specific gradient; Brunner uses its published single-cell EvoSep gradient. Unreported handling and turnover overhead are excluded. These are nominal hours, not measured lab scheduling times.

Whole-PBMC donor variance transfer remains an assumption, including donor composition differences. k values are registered compression assumptions. Passing M2 does not validate power or donor transfer, and the registered M4 criterion remains unchanged.

Outputs: [full power grid](power-grid.csv), [FDR table](realised-fdr.csv), [instrument-hour table](instrument-hours.csv), [cheapest designs](cheapest-designs.csv), [compressed replicate metrics](replicate-metrics.csv.gz), [manifest](run-manifest.json).

Software: R=4.6.1; scp=1.22.0; limma=3.68.5; edgeR=4.10.5; scplainer=workflow in scp 1.22.0; no separate package
