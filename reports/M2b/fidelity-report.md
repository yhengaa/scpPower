# M2b simulator fidelity report — amendment I
13 October 2026. Frozen addendum v1.1.5. Every dataset uses its exact cell-line-by-acquisition-batch counts and beta=0. Each comparison averages the registered independent replicates. Every applicable tier must pass; no failing tier is averaged away.

| Dataset | Gate | Missingness maximum error (pp) | Coverage relative error | Cell SD maximum relative error | Batch SD maximum relative error |
| --- | --- | ---: | ---: | ---: | ---: |
| specht2019v3 | PASS | 0.90 | 1.7% | 3.7% | 7.2% |
| leduc2022_pSCoPE | PASS | 0.53 | 0.8% | 0.7% | 11.3% |

Registered tolerances: missing fraction per tier within 5 percentage points; proteins per cell within 10% relative; median within-line cell SD per tier within 20% relative; between-batch SD per tier within 25% relative where separately identified. Only the two TMT datasets are rerun.

The [full fidelity table](fidelity-table.csv) gives observed values, simulated means, errors, tolerances and PASS/FAIL for every quantity and tier. [Replicate metrics](replicate-metrics.csv) preserve all draws and seeds.

Suspected causes are interpretations of these fixed outcomes, not changes to the implementation:

- specht2019v3: PASS on all applicable quantities.
- leduc2022_pSCoPE: PASS on all applicable quantities.

The accepted M1 estimates, including Specht's weak negative raw-intensity slope, are unchanged. Amendment G replaces TMT curve masking with plex p_g then conditional channel q_g draws; I changes only their mask source to the endpoint assay, independent of effect size. For DIA and label-free data, protein r_g supplies baseline detection and the fitted curve supplies only the slope. Parameters were fixed before this rerun; no adjustment follows its outcomes. If any series fails, M4b and the v0.3.1 tag remain unauthorised and work stops after this report.

Missingness and coverage use the endpoint aggregated assay mask for TMT under I ; non-TMT results are unchanged and not rerun. Cell SD uses the accepted variance assay after per-cell median shifts; the same shift is applied to simulations. It is the median over measured protein-by-line SDs within each tier. Batch SD is the median of protein-level SDs across mixed-batch means after line-mean subtraction. The missingness figure shows fixed-tier gate statistics and the registered 5-pp band; no smoothing or calibration adjustment is applied.

Cell lines are culture samples. The injected Walzik variance remains a whole-PBMC proxy, including donor composition differences; M2 cannot validate bulk-to-cell-line donor transfer. The registered M4 criterion is unchanged.

Software: R=4.6.1; scp=1.22.0; limma=3.68.5; edgeR=4.10.5; scplainer=workflow in scp 1.22.0; no separate package

## Additional descriptive quantity (outside the gate)

| Dataset | Observed fraction quantified in at least half the cells | Simulated mean | Replicate SD |
| --- | ---: | ---: | ---: |
| specht2019v3 | 0.1614 | 0.1553 | 0.0064 |
| leduc2022_pSCoPE | 0.2606 | 0.2589 | 0.0070 |

This descriptive quantity does not alter the gate. TimsTOF, Q Exactive and Brunner are not rerun. This quantity has no tolerance and cannot change PASS/FAIL. Registered consequence: TMT carrier pools both conditions, so p_g and q_g do not respond to effect size; testability depends on plex count and barely on effect size.

Per-protein parameter tables and preset hashes are saved alongside the report. Undefined means for never-detected proteins use a documented tier fallback while retaining zero detection; unmapped donor variances and missing protein component estimates use the registered tier fallback. The calibration itself is not changed.

The retained legacy intensity-bin CSV uses observed endpoint ratio means and
simulated raw-intensity means, which are on different scales. It is not a gate
input and is not used in the M2b plots. The fixed abundance-tier statistics and
plots are the relevant missingness comparisons under I.
