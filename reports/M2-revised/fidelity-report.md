# Revised M2 simulator fidelity report — amendment G
13 October 2026. Frozen addendum v1.1.3. Every dataset uses its exact cell-line-by-acquisition-batch counts and beta=0. Each comparison averages the registered independent replicates. Every applicable tier must pass; no failing tier is averaged away.

| Dataset | Gate | Missingness maximum error (pp) | Coverage relative error | Cell SD maximum relative error | Batch SD maximum relative error |
| --- | --- | ---: | ---: | ---: | ---: |
| specht2019v3 | PASS | 0.75 | 1.0% | 2.8% | 4.3% |
| leduc2022_pSCoPE | PASS | 0.60 | 0.8% | 2.5% | 9.2% |
| derks2022_timsTOF | PASS | 1.24 | 0.2% | 4.8% | 7.6% |
| derks2022_QExactive | PASS | 0.59 | 0.5% | 4.0% | 6.3% |
| brunner2022 | PASS | 1.04 | 0.5% | 0.3% | Not applicable (F) |

Registered tolerances: missing fraction per tier within 5 percentage points; proteins per cell within 10% relative; median within-line cell SD per tier within 20% relative; between-batch SD per tier within 25% relative where separately identified. Brunner has no separate batch component under F.

The [full fidelity table](fidelity-table.csv) gives observed values, simulated means, errors, tolerances and PASS/FAIL for every quantity and tier. [Replicate metrics](replicate-metrics.csv) preserve all draws and seeds.

Suspected causes are interpretations of these fixed outcomes, not changes to the implementation:

- specht2019v3: PASS on all applicable quantities.
- leduc2022_pSCoPE: PASS on all applicable quantities.
- derks2022_timsTOF: PASS on all applicable quantities.
- derks2022_QExactive: PASS on all applicable quantities.
- brunner2022: PASS on all applicable quantities.

The accepted M1 estimates, including Specht's weak negative raw-intensity slope, are unchanged. Amendment G replaces TMT curve masking with plex p_g then conditional channel q_g draws, independent of effect size. For DIA and label-free data, protein r_g supplies baseline detection and the fitted curve supplies only the slope. Parameters were fixed before this rerun; no adjustment follows its outcomes. If any series fails, M3 and the v0.2.0 tag remain unauthorised and work stops after this report.

Missingness and coverage use the paired raw PSM mask for TMT and the selected quantitative protein mask for DIA. Cell SD uses the accepted variance assay after per-cell median shifts; the same shift is applied to simulations. It is the median over measured protein-by-line SDs within each tier. Batch SD is the median of protein-level SDs across mixed-batch means after line-mean subtraction. The figures bin intensities for display only; they do not alter the curves or gate tests.

Cell lines are culture samples. The injected Walzik variance remains a whole-PBMC proxy, including donor composition differences; M2 cannot validate bulk-to-cell-line donor transfer. The registered M4 criterion is unchanged.

Software: R=4.6.1; scp=1.22.0; limma=3.68.5; edgeR=4.10.5; scplainer=workflow in scp 1.22.0; no separate package

## Additional descriptive quantity (outside the gate)

| Dataset | Observed fraction quantified in at least half the cells | Simulated mean | Replicate SD |
| --- | ---: | ---: | ---: |
| specht2019v3 | 0.2249 | 0.2159 | 0.0063 |
| leduc2022_pSCoPE | 0.3150 | 0.3143 | 0.0068 |
| derks2022_timsTOF | 0.6939 | 0.7788 | 0.0074 |
| derks2022_QExactive | 0.3680 | 0.3977 | 0.0131 |
| brunner2022 | 0.4867 | 0.4952 | 0.0053 |

The timsTOF fraction quantified in half the cells differs visibly despite passing the registered tier-level metrics; the fixed model does not reproduce all protein-level detection structure. This quantity has no tolerance and cannot change PASS/FAIL. Registered consequence: TMT carrier pools both conditions, so p_g and q_g do not respond to effect size; testability depends on plex count and barely on effect size.

Per-protein parameter tables and preset hashes are saved alongside the report. Undefined means for never-detected proteins use a documented tier fallback while retaining zero detection; unmapped donor variances and missing protein component estimates use the registered tier fallback. The calibration itself is not changed.
