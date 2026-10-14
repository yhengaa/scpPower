# Diagnostic J — TMT endpoint mask consistency
13 October 2026. Uses the v0.3.0 validation halves and H1 source cells; no preset changes preceded this diagnostic.

| Dataset | PSM finite | Endpoint finite | Endpoint sparsity gap | PSM finite / endpoint NA pairs | Median p difference | Median q difference |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| specht2019v3 | 30.193% | 20.412% | 9.781 pp | 10.037% | -0.06977 | 0.00000 |
| leduc2022_pSCoPE | 34.844% | 27.764% | 7.080 pp | 7.744% | -0.05263 | -0.01272 |

## Decision before any change
Amendment I applies to both TMT datasets under the fixed rule: the endpoint aggregated matrix is sparser than the PSM mask by 9.781 pp in Specht macrophages and 7.080 pp in Leduc monocytes, exceeding the 5-pp threshold in each dataset. Detection parameters and the M2 missing-fraction comparison must therefore use the mask of the quantitative assay evaluated by the endpoint. This is an assay-consistency change, not a refit of the noise model or a change to any tolerance. The accepted v0.3.0 M4 verdict stands. M2b must pass on both TMT datasets before M4b proceeds.

## Full comparisons
[Overall and per-tier finite fractions and paired p/q differences](mask-summary-by-tier.csv); [Specht per-protein p/q](specht2019v3-protein-mask-differences.csv); [Leduc per-protein p/q](leduc2022_pSCoPE-protein-mask-differences.csv). Differences are endpoint minus PSM. Rates are estimated on the same validation proteins, source cells and native runs. p is the unweighted fraction of occupied native plexes with any finite protein value; q is quantified channels divided by channels in those identified plexes. The two masks are not nested; reverse discordant pairs are reported as well.

## Pipeline attribution
The endpoint finite mask is reproduced exactly from the deposited peptide assay by the registered log transformation and protein median. A PSM-positive / endpoint-NA pair has no finite retained peptide quantity in that assay. The local protein median does not remove a finite included peptide, and no compared cell is removed between the two paired masks. The identifiable loss is upstream of the local protein aggregation: peptide retention/quantification and the producer ratio processing, not a newly applied local cell filter.
Specht’s producer exports Peptides-raw.csv from cr_norm_log(t3), after reference-channel normalisation, cell QC and peptide missing-data filtering. Leduc’s scpdata builder imports a processed peptide table and builds its peptide/protein-group mapping; the two deposited objects do not identify which particular upstream filter removed every discordant pair. We do not assign a filter to a pair without a retention trace. [Pipeline audit](pipeline-audit.csv).
Sources: [Specht producer code](https://github.com/cvanderaa/SCoPE2/blob/master/code/SCoPE2_analysis.R); [Specht scpdata builder](https://github.com/UCLouvain-CBIO/scpdata/blob/master/inst/scripts/make-data_specht2019v3.R); [Leduc scpdata builder](https://github.com/UCLouvain-CBIO/scpdata/blob/master/inst/scripts/make-data_leduc2022_pSCoPE.R). Downloaded source hashes are retained in the existing source manifests.

## Native cell-depth variation
| Dataset | Mask | Observed depth SD | Simulated mean depth SD |
| --- | --- | ---: | ---: |
| specht2019v3 | PSM | 69.01 | 13.39 |
| specht2019v3 | aggregated | 54.94 | 10.62 |
| leduc2022_pSCoPE | PSM | 95.16 | 13.36 |
| leduc2022_pSCoPE | aggregated | 91.29 | 12.17 |

Depth is quantified validation proteins per cell. Each mask is compared with twenty native-layout, detection-only G2 replicates using that mask’s own per-protein p/q, fixed protein identities and independent plex/channel Bernoulli draws. TMT masking is intensity-independent, so variance draws are unnecessary for this mask diagnostic. This is a descriptive diagnostic, not a new gate or a calibration half-preset. [Cell depths](specht2019v3-cell-depth.csv), [Leduc depths](leduc2022_pSCoPE-cell-depth.csv), [depth summary](native-depth-summary.csv), [replicate records](depth-replicate-metrics.csv).
Both masks show much greater native cell-depth variation than independent G2 draws. Amendment I addresses the demonstrated assay mismatch; it does not establish that cell/channel dependence is resolved. No dependence term is introduced under this directive.

Software: R=4.6.1; scp=1.22.0; limma=3.68.5; edgeR=4.10.5; scplainer=workflow in scp 1.22.0; no separate package

[Run manifest](run-manifest.json). No M5 work.
