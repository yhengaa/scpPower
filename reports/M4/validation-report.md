# M4 semi-synthetic validation

M4 validates cell-level noise, missingness and the primary test under the registered semi-synthetic construction. It cannot validate the transfer of between-donor variance from Walzik whole PBMC to these cell lines, because that variance is injected, or carrier compression k, which remains a user assumption. All predictions were completed and hashed before any observed endpoint was evaluated.

13 October 2026. Frozen addendum v1.1.4. Overall verdict: **PASS** (4 datasets pass; at least two required).

## Registered gate
A dataset passes when the absolute predicted–observed unconditional power gap is at most 10 percentage points at at least 80% of eligible points with m >= 3. Every counted point uses 200 prediction and 200 observation replicates. Points with m=1 or 2 are listed as deterministic zero-power exclusions; they are not simulated and do not enter the pass fraction. H1 makes timsTOF not evaluable, not a pass. No model, calibration, tolerance or criterion was changed in response to the results.

| Dataset | Source line | Cells | Counted points | Within 10 pp | Pass fraction | Verdict |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| specht2019v3 | Macrophage | 1096 | 72 | 67 | 93.1% | PASS |
| leduc2022_pSCoPE | Monocyte | 783 | 69 | 66 | 95.7% | PASS |
| derks2022_timsTOF | Melanoma_t | 10 | 0 | 0 | NA | NOT EVALUABLE |
| derks2022_QExactive | PDAC | 45 | 6 | 6 | 100.0% | PASS |
| brunner2022 | HeLa | 434 | 57 | 57 | 100.0% | PASS |

## Construction and separation
Protein identities from the accepted G presets were split once per dataset, with the split seed and both halves recorded. Half-presets were refitted from the calibration proteins across the native accepted SCP calibration cells. H1 restricts validation source cells to the largest single line; it does not remove other native calibration cells from the training protein half. The released presets remain byte-identical. All tier summaries, missingness curves, p_g/q_g or r_g, variance estimates and cell median shifts in each new preset use calibration proteins only.
The validation half retains native quantitative missingness and plex/run labels. Calibration-half cell median shifts are applied before injection; neither injected condition effects nor donor effects are normalised out. No missing value is filled and no cell-level noise or batch effect is newly drawn in the observed data: those components are already present in real cells. For each replicate, cells are randomly sampled without replacement and assigned to 2n pseudo-samples of m cells; pseudo-samples are assigned to A/B. Ten percent of validation proteins are randomly changed and receive beta in B. Independent normal donor effects use each protein’s external mapped Walzik variance, with the half-preset tier value for unmapped IDs. Frozen M1 tier labels select that fallback for validation proteins; no validation intensities are used to fit the half-preset.
The first assignment at every eligible point is recorded as cell IDs and a pseudo-sample-by-native-batch count matrix. The prediction uses that exact custom matrix. Later observed replicates use independent source assignments; their seeds and assignment hashes are recorded. Brunner has one cell per batch. k=1 is used throughout M4. Predictions jointly resample training-half protein rows under G1; observations evaluate held-out real protein rows, rather than copying or fitting their cell-noise estimates.

## Endpoint and uncertainty
The primary endpoint remains per-sample means of log2 values, at least three quantified cells in every sample, limma-trend with design ~ condition, BH 0.05 and injected-sign agreement. Unconditional power divides detected changed proteins by all changed proteins. Conditional power divides by testable changed proteins, is undefined when none are testable, and is averaged over defined replicates. FDR remains mean null-protein FDP among BH rejections, zero when nothing is rejected. Mean rejections and zero-rejection fractions are reported beside it.
Tables provide prediction and observation power SE separately, plus their quadrature sum for the gap. The gate compares mean powers with the fixed 10-pp tolerance, not confidence intervals. The plot band is ±10 pp. Finite Monte Carlo error and repeated reuse of the finite source cell pool limit interpretation; this is not an independent-donor experiment.

## Disagreements and suspected causes

### specht2019v3
PASS: 67/72 counted points within 10 pp. Maximum gap 19.5 pp. Maximum observed power SE 0.48 pp.
Observed mean rejections range from 0.030 to 49.110 per replicate; the no-rejection fraction ranges from 0.0% to 97.0%. Realised mean FDP ranges from 0.0175 to 0.0708; 18 counted points exceed nominal 0.05, without adjustment.
0 discrepant points have observed power above prediction; 5 have prediction above observation. Mean testability in the discrepant points is 0.519 predicted versus 0.319 observed.
Suspected causes: real source-line plex and channel masks can have dependencies beyond independent p_g/q_g draws; training across the native lines can differ from the selected validation line. TMT raw-PSM detection and the endpoint’s ratio-based quantitative assay can also have different finite-value masks. These are interpretations of held-out discrepancies, not fitted corrections.
[Eligible-point table](specht2019v3-validation-table.csv) · [Comparison plot](specht2019v3-predicted-observed.png)

### leduc2022_pSCoPE
PASS: 66/69 counted points within 10 pp. Maximum gap 12.8 pp. Maximum observed power SE 0.35 pp.
Observed mean rejections range from 0.020 to 63.895 per replicate; the no-rejection fraction ranges from 0.0% to 98.0%. Realised mean FDP ranges from 0.0067 to 0.0800; 21 counted points exceed nominal 0.05, without adjustment.
0 discrepant points have observed power above prediction; 3 have prediction above observation. Mean testability in the discrepant points is 0.553 predicted versus 0.411 observed.
Suspected causes: real source-line plex and channel masks can have dependencies beyond independent p_g/q_g draws; training across the native lines can differ from the selected validation line. TMT raw-PSM detection and the endpoint’s ratio-based quantitative assay can also have different finite-value masks. These are interpretations of held-out discrepancies, not fitted corrections.
[Eligible-point table](leduc2022_pSCoPE-validation-table.csv) · [Comparison plot](leduc2022_pSCoPE-predicted-observed.png)

### derks2022_timsTOF
Not evaluable: ten cells per source line cannot supply the smallest counted design (2 × 3 × 4 = 24 cells). No degenerate timsTOF simulations were run.

### derks2022_QExactive
PASS: 6/6 counted points within 10 pp. Maximum gap 0.4 pp. Maximum observed power SE 0.18 pp.
Observed mean rejections range from 0.110 to 1.115 per replicate; the no-rejection fraction ranges from 56.0% to 94.5%. Realised mean FDP ranges from 0.0380 to 0.0912; 5 counted points exceed nominal 0.05, without adjustment.
Q Exactive passes at very low power: every observed mean is below 2%, with at most 1.115 mean rejections per replicate. The absolute 10-pp criterion provides limited evidence in this region and does not establish agreement at useful high power. This limitation does not change its registered verdict.
All counted points meet the registered 10-pp tolerance. This does not validate the donor-transfer or compression assumptions, or every feature of cell-level dependence.
[Eligible-point table](derks2022_QExactive-validation-table.csv) · [Comparison plot](derks2022_QExactive-predicted-observed.png)

### brunner2022
PASS: 57/57 counted points within 10 pp. Maximum gap 5.4 pp. Maximum observed power SE 0.69 pp.
Observed mean rejections range from 0.020 to 78.055 per replicate; the no-rejection fraction ranges from 0.0% to 98.5%. Realised mean FDP ranges from 0.0125 to 0.0950; 12 counted points exceed nominal 0.05, without adjustment.
All counted points meet the registered 10-pp tolerance. This does not validate the donor-transfer or compression assumptions, or every feature of cell-level dependence.
[Eligible-point table](brunner2022-validation-table.csv) · [Comparison plot](brunner2022-predicted-observed.png)

## Rejection counts and stop decision
Read the FDR column with mean rejections and the zero-rejection fraction. Constructional zero-power points have zero rejections and zero FDP by definition; their zero-rejection fraction is one. Low FDR alone does not demonstrate useful detection.
The overall registered criterion passes. The release is authorised after engineering checks; dataset failures and all discrepant points remain in the evidence.

No M5 work was performed. [Full eligible-point table](validation-table.csv), [dataset verdicts](dataset-verdicts.csv), [replicate metrics](replicate-metrics.csv.gz), [manifest](run-manifest.json).

Software: R=4.6.1; scp=1.22.0; limma=3.68.5; edgeR=4.10.5; scplainer=workflow in scp 1.22.0; no separate package


## M4 closure and limitation K — addendum v1.1.5

In carrier-TMT presets, unconditional power is validated only up to observed power of about 0.3. At twofold, v0.3.0 predictions ran above observation by up to 20 percentage points, associated with over-predicted testability. Q Exactive and timsTOF could not be validated at useful power for lack of cells. No TMT design numbers or new design examples are presented in the README or a vignette during this round.

The 1-pp conditional-power description holds for mean signed TMT gaps: -0.82 pp for Specht and +0.94 pp for Leduc across counted points. Brunner averages +1.30 pp and Q Exactive -0.46 pp. It is not a bound on every point: among points with maximum power at least 0.2, maximum conditional gaps are 1.53 pp for Specht, 4.56 pp for Leduc and 6.05 pp for Brunner. Across all counted points the maxima are 7.77, 11.13, 7.57 and 2.32 pp respectively. The high-power agreement counts do reproduce: 3/8 Specht, 11/14 Leduc and 12/12 Brunner; Q Exactive has none. Testability is the main TMT discrepancy, while the statement that the entire gap is testability with universal agreement within 1 pp is stronger than the recorded columns support. This clarification does not change any data, criterion or accepted verdict.

[Reproducible power-window and conditional-gap table](../J/v0.3.0-power-window-diagnostics.csv). The original v0.3.0 report is preserved in validation-report-v0.3.0.md; its original manifest remains the receipt for that snapshot. No M4 numerical result or verdict is changed.
