# M2 simulator fidelity report
13 October 2026. Frozen addendum v1.1.2. Every dataset uses its exact cell-line-by-acquisition-batch counts and beta=0. Each comparison averages the registered independent replicates. Every applicable tier must pass; no failing tier is averaged away.

| Dataset | Gate | Missingness maximum error (pp) | Coverage relative error | Cell SD maximum relative error | Batch SD maximum relative error |
| --- | --- | ---: | ---: | ---: | ---: |
| specht2019v3 | FAIL | 22.73 | 0.0% | 25.0% | 39.8% |
| leduc2022_pSCoPE | FAIL | 15.48 | 2.5% | 29.3% | 61.0% |
| derks2022_timsTOF | FAIL | 5.87 | 1.3% | 28.2% | 33.9% |
| derks2022_QExactive | FAIL | 5.84 | 0.7% | 11.4% | 25.5% |
| brunner2022 | PASS | 4.97 | 0.1% | 13.3% | Not applicable (F) |

Registered tolerances: missing fraction per tier within 5 percentage points; proteins per cell within 10% relative; median within-line cell SD per tier within 20% relative; between-batch SD per tier within 25% relative where separately identified. Brunner has no separate batch component under F.

The [full fidelity table](fidelity-table.csv) gives observed values, simulated means, errors, tolerances and PASS/FAIL for every quantity and tier. [Replicate metrics](replicate-metrics.csv) preserve all draws and seeds.

Suspected causes are interpretations of these fixed outcomes, not changes to the implementation:

- specht2019v3: FAIL; one global intensity curve cannot reproduce all retained tier masks; TMT tiers retain the accepted ratio scale while their predictor is raw intensity, and identification availability is heterogeneous; the unweighted mean of heterogeneous protein variances generates an RMS scale and need not reconstruct the median measured within-line SD; observed selection and normalization also differ; pooled latent batch variance need not reproduce median noisy batch-mean SD, especially with sparse observation and differing cell-line means.
- leduc2022_pSCoPE: FAIL; one global intensity curve cannot reproduce all retained tier masks; TMT tiers retain the accepted ratio scale while their predictor is raw intensity, and identification availability is heterogeneous; the unweighted mean of heterogeneous protein variances generates an RMS scale and need not reconstruct the median measured within-line SD; observed selection and normalization also differ; pooled latent batch variance need not reproduce median noisy batch-mean SD, especially with sparse observation and differing cell-line means.
- derks2022_timsTOF: FAIL; resampling observed means and adding sample/batch effects can alter tier-specific coverage; an intensity-only curve also omits protein-specific identification availability; the unweighted mean of heterogeneous protein variances generates an RMS scale and need not reconstruct the median measured within-line SD; observed selection and normalization also differ; pooled latent batch variance need not reproduce median noisy batch-mean SD, especially with sparse observation and differing cell-line means.
- derks2022_QExactive: FAIL; resampling observed means and adding sample/batch effects can alter tier-specific coverage; an intensity-only curve also omits protein-specific identification availability; pooled latent batch variance need not reproduce median noisy batch-mean SD, especially with sparse observation and differing cell-line means.
- brunner2022: PASS on all applicable quantities.

Specht retains a weak negative fitted raw-intensity slope, flagged in final M1. It was not replaced by MCAR or forced positive for this gate. Failed datasets stop at this report. No parameter, preset, missingness definition, or tolerance was altered to obtain a pass.

Missingness and coverage use the paired raw PSM mask for TMT and the selected quantitative protein mask for DIA. Cell SD uses the accepted variance assay after per-cell median shifts; the same shift is applied to simulations. It is the median over measured protein-by-line SDs within each tier. Batch SD is the median of protein-level SDs across mixed-batch means after line-mean subtraction. The figures bin intensities for display only; they do not alter the curves or gate tests.

Cell lines are culture samples. The injected Walzik variance remains a whole-PBMC proxy, including donor composition differences; M2 cannot validate bulk-to-cell-line donor transfer. The registered M4 criterion is unchanged.

Software: R=4.6.1; scp=1.22.0; limma=3.68.5; edgeR=4.10.5; scplainer=workflow in scp 1.22.0; no separate package
