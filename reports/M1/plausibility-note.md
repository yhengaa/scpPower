# scpPower M1 plausibility note

13 October 2026

**Disposition: STOPPED. M1 is not closed. The model and registered M4 criterion were not changed.**

**Leduc pSCoPE.** `proteins_norm2` retains missing values. After removing 1042 rows annotated as reverse/contaminant groups, 1802 proteins across 1543 cells were examined. The missing fraction is 61.22%. The quantification slope is -0.1599 (95% interval -0.1727 to -0.1470), so the fitted missing fraction rises with expected intensity. This is flagged as implausible for an abundance-dependent detection preset and triggers the stop rule. Protein-wise normalization has removed the original between-protein abundance scale; the quintiles describe this supplied normalized matrix. The exclusion rule and result are preserved without another fit. Scpdata documents that `proteins_processed` is imputed and batch corrected; it was not used.

**Specht SCoPE2.** The downloaded `proteins` assay is complete and contains the dimensions used for the paper's imputed analysis. The paper explicitly describes missing quantifications and nearest-neighbor imputation. Treating the complete matrix as an original missingness estimate would disagree with that evidence. This assay is reported as uncalibratable for original missingness and variance; no substitute dataset or new protein assay was introduced.

**Derks plexDIA and Brunner label-free DIA.** Both inputs were downloaded and inspected before the stop, but neither was calibrated afterward. Derks' supplied protein table is complete; its original quantification mask cannot be certified. The source paper's high within-plex completeness is not global matrix completeness. Brunner's protein-group matrix retains missing values and the source paper describes abundance-related detection, but no fitted estimate is claimed here. No cross-instrument estimates were pooled and no bulk Derks columns were calibrated.

**Variance and biological units.** The selected cohort is Walzik et al., PXD058573: 23 healthy PBMC donors with 46 baseline observations, mapped using the checksum-verified author supplement. Leduc has 1468 exact identifier overlaps and 1455 proteins with variance estimable from at least four quantified donors. Its pooled mapped log2 donor variances range from 0.1670 to 0.2074. These are PBMC proxies, not established donor variances for the SCP cell lines. The bulk paper's reported CVs use different processing and are not directly comparable to these unimputed log2 variances. Transfer requires a matching-cell-type lab pilot and cannot be validated by M4. Verified donor/culture IDs are absent from the SCP metadata, so cell and batch components and biological sample counts are NA; technical slides and runs were not promoted to biological samples.

**Carrier effect.** Cheung's PMID resolves to Nature Methods, not J Proteome Res. Its paper demonstrates dependence on carrier level and ion sampling. Publisher spreadsheet requests returned browser challenges or failed connections, so no defensible per-ratio default k was obtained. Compression is NA with a reason; no single-carrier coefficient, correlation-as-compression proxy, or uncompressed default was substituted. Power comparisons belong to a later directive and were not run.

Sources: [Specht](https://doi.org/10.1186/s13059-021-02267-5), [Leduc assay documentation](https://bioconductor.org/packages/release/data/experiment/html/scpdata.html), [Derks](https://doi.org/10.1038/s41587-022-01389-w), [Brunner](https://doi.org/10.15252/msb.202110798), [Walzik](https://doi.org/10.1038/s41467-025-68101-9), [Cheung](https://doi.org/10.1038/s41592-020-01002-5).

Software: R=4.6.1; scp=1.22.0; limma=3.68.5; edgeR=4.10.5; scplainer=workflow in scp 1.22.0; no separate package

Reproduction: public downloads and hashes are handled by `fetch_sources.py`, `fetch_bulk_metadata.py`, and `inspect_scp.R`; `calibrate_m1.R` regenerates the preserved fit and stop; this script renders the note and plots without fitting any further models.
