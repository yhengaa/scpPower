---
title: "scpPower: calibrated power and sample-size calculations for two-condition single-cell proteomics experiments"
author: "Youshan Heng^1^*"
---

^1^ [AFFILIATION, CITY, COUNTRY]. ORCID 0000-0002-2758-171X.
\* Corresponding author: poggar@hotmail.com

**Article type:** Technical Note

**Keywords:** single-cell proteomics; experimental design; power analysis; sample size; missing values; isobaric carrier; plexDIA; label-free DIA; simulation

# Abstract

Single-cell proteomics (SCP) by mass spectrometry now supports studies that compare conditions across independent biological samples, but no tool exists to decide how many samples and cells such a study needs. Generic calculators fail because SCP missingness is structured by plex and carrier rather than by cell, because a protein can only be tested where enough cells are quantified in every sample, and because between-sample variance dominates once cells are averaged. We present scpPower, an R package that calibrates a hierarchical model of protein intensity, batch noise, cell noise and quantification probability from any unimputed SCP dataset, maps between-donor variance from bulk proteomics, and estimates power, testability, false discovery rate and instrument hours for a design grid by simulation under a per-sample limma endpoint. Five presets calibrated from four public datasets cover carrier TMT, plexDIA and label-free DIA. Acceptance criteria were fixed before each development stage. The simulator reproduces the source data within registered tolerances and predicts power on held-out real proteins with injected effects within 3.3 percentage points on two TMT datasets and 5.4 on a label-free dataset. The registered grid shows that 1.25-fold changes are undetectable with twelve samples per condition, that plex count rather than cell count limits carrier-TMT designs, and that carrier compression halves TMT power. scpPower is available at https://github.com/yhengaa/scpPower under GPL-3.

# Introduction

Mass spectrometry-based single-cell proteomics has moved from method papers on cultured cell lines^1–4^ to primary tissue, where hepatocyte zonation has been mapped in liver^5,7^ and spinal motor neurons profiled from patients with amyotrophic lateral sclerosis,^6^ and such studies compare biological conditions across donors, animals or patients. Such studies are expensive: a single-cell run occupies the instrument for 30 to 95 minutes and a design with twelve samples per condition and a few dozen cells per sample consumes hundreds of instrument hours. The first question such a study raises, how many samples and how many cells per sample are needed to detect a given change, currently has no tool. Investigators design by precedent.

Power calculators exist for neighbouring fields. scPower^8^ answers the question for multi-sample single-cell RNA sequencing with a calibrated model of count noise and detection, and is the model for the present work. Boekweg et al.^9^ computed sample sizes for temporal dynamics in SCP with a set of scripts tied to a time-course design. peppwR^10^ fits abundance distributions and runs power analysis for bulk phosphopeptide data. The scp and scpdata Bioconductor packages^11,12^ standardise SCP data handling and curate public datasets, and scplainer^13^ models cell-level SCP data with linear models, but neither addresses design. None of these transfers to a two-condition SCP study, for three reasons.

First, SCP missingness is not a cell-level coin toss. In carrier-multiplexed TMT, peptide identification is driven by the pooled carrier channel,^14^ so a protein is either identified in a whole plex or absent from every cell in it; the single-cell channels then decide only whether an identified protein is quantified in a given cell. Second, a protein can be tested only where enough of its cells are quantified in every biological sample. The fraction of testable proteins therefore depends on how cells are laid out across plexes or runs, not just on how many cells there are. Third, the replication unit for a biological claim is the sample, not the cell, and once cells are averaged within a sample the between-sample variance usually exceeds the cell-level noise.

scpPower addresses these three points directly. It calibrates per-protein noise and detection parameters from real SCP data, maps between-donor variance from an independent bulk cohort, simulates experiments that reproduce the missingness structure of each labelling scheme, and reports power, testability, false discovery rate and instrument hours for any design. The package was developed under a frozen mandate with acceptance criteria fixed before each stage, and every failure is reported in the repository rather than tuned away.

# Methods

## Model

Each log2 intensity of protein *g* in cell *c* of biological sample *s* acquired in batch *b* is modelled as

*y*~gsbc~ = *μ*~g~ + *k* *β*~g~ *C*~s~ + *u*~gs~ + *v*~gb~ + *ε*~gsbc~,

where *μ*~g~ is the protein mean, *β*~g~ the log2 condition effect for the changed proteins, *C*~s~ the condition indicator, *k* the carrier-compression factor, and *u*, *v* and *ε* independent normal sample, batch and cell effects with per-protein variances *σ*^2^~s,g~, *σ*^2^~b,g~ and *σ*^2^~e,g~ (Figure 1a). Proteins are grouped into quintiles of mean log2 intensity (tiers) for reporting and as fallbacks, but the generating unit is a resampled real protein: each simulated protein draws one calibrated protein from the preset together with its mean, its three variances and its detection parameters, so the joint distribution of these quantities is preserved.

## Calibration data

Five series were calibrated from four public datasets curated in scpdata 1.20.0^12^ (Table 1): the SCoPE2 macrophage and monocyte data of Specht et al.^1^ (carrier TMT, Q Exactive), the pSCoPE monocyte and melanoma data of Leduc et al.^2^ (carrier TMT, Q Exactive, TMTpro 18-plex), the plexDIA data of Derks et al.^3^ on a timsTOF SCP and on a Q Exactive, calibrated separately, and the label-free DIA HeLa data of Brunner et al.^4^ acquired one cell per run. Variance components were estimated from the earliest unimputed, uncorrected quantitative assay of each dataset: peptide-level assays aggregated to proteins by the per-cell median for the TMT and plexDIA data, and the raw protein table for Brunner. Cell lines are the biological units of these datasets; *σ*^2^~e,g~ is the residual cell variance within line and batch after a per-cell median shift, and *σ*^2^~b,g~ is the variance of batch means after removing line means, estimated only where a batch mixes lines (plexes for TMT, runs for plexDIA). Where every cell is its own run, batch and cell noise are the same unit and *σ*^2^~b,g~ is absorbed into *σ*^2^~e,g~.

## Missingness

Two detection models are used, chosen by labelling scheme (Figure 1a). For carrier TMT, a protein is identified in a plex with probability *p*~g~ and, if identified, quantified in each single-cell channel with probability *q*~g~; otherwise the whole plex is missing for that protein. Both probabilities are estimated per protein from the quantified mask of the assay the statistical endpoint uses, and they do not respond to the condition effect, because the carrier that drives identification pools both conditions. For plexDIA and label-free DIA, each protein keeps its observed quantification rate *r*~g~ as a baseline and the logit of detection in a cell shifts by a dataset-level slope times the difference between the cell's expected log2 intensity and the protein's calibrated mean; the condition effect therefore changes detection in the changed condition. No imputation is performed anywhere in calibration or simulation.

## Between-donor variance and carrier compression

Public SCP datasets are cell lines with one culture per line, so between-donor variance cannot be estimated from them. It was mapped per protein from the bulk proteomics of Walzik et al.,^15^ 23 healthy donors, whole peripheral blood mononuclear cells, as the unbiased variance of log2 donor means at baseline, matched by protein identifier (1,539 to 2,545 proteins per series) with the tier median as fallback for unmapped proteins. The transfer of this variance to single cells of a given cell type, and the inclusion of between-donor composition differences in whole-PBMC data, are declared assumptions tested only by a pilot in the intended cell type. The carrier-compression factor *k* is likewise a declared assumption: no public dataset measures it, so the default is *k* = 1 with sensitivity values 0.75 and 0.5 reported for the TMT presets.

## Endpoint, design grid and instrument time

The primary endpoint is the per-sample mean of log2 intensity across the cells of each biological sample. A protein is testable when at least three of its cells are quantified in every sample. Testable proteins are tested with limma-trend^16^ under the design ~ condition at a Benjamini–Hochberg false discovery rate of 0.05, and a changed protein is detected when its adjusted *p* is below 0.05 and the estimated sign matches the injected sign. Power is the fraction of all changed proteins detected; conditional power restricts the denominator to testable changed proteins; testability is the fraction of all proteins meeting the cell minimum; realised false discovery rate is the mean false discovery proportion among rejections, reported with the mean number of rejections and the fraction of replicates with none. The registered grid is *m* ∈ {1, 2, 4, 8, 16, 32, 64} cells per sample, *n* ∈ {3, 4, 6, 8, 12} samples per condition and fold changes 1.25, 1.5 and 2, with 10% of proteins changed, 200 replicates per point, for every preset and for *k* ∈ {1, 0.75, 0.5} on the TMT presets (945 points). Cells are laid out in balanced acquisition batches (plexes of 12 or 14 single-cell channels, plexDIA runs of two cells under an equal-condition rule, or one cell per run), runs are counted as occupied batches, and nominal instrument hours are runs times the published method length of each preset (95, 95, 30, 53 and 35 minutes; Table 1), excluding handling and turnover.

## Validation

Two acceptance gates were registered before the corresponding stage was run. The fidelity gate requires the simulator, run at each dataset's own design with no condition effect and 20 replicates, to reproduce four statistics of the source data within fixed tolerances: missing fraction per tier within 5 percentage points, proteins per cell within 10%, median per-protein cell SD per tier within 20%, and median batch SD per tier within 25%. The semi-synthetic validation gate splits each dataset's proteins into calibration and validation halves by a recorded seed, refits a preset from the calibration half only, and then, on real cells of the largest single cell line, repeatedly assigns cells at random to 2*n* pseudo-samples of *m* cells, assigns pseudo-samples to two conditions, injects a log2 effect into 10% of validation proteins in one condition and an independent normal donor effect per pseudo-sample and protein, and evaluates the primary endpoint. Predicted power from the half-preset, computed with the realised pseudo-sample-by-batch layout and hashed before any observation, must lie within 10 percentage points of observed power at 80% of eligible grid points (2*nm* not exceeding the cell count, *m* ≥ 3) on at least two datasets. The validation tests cell-level noise, missingness and the test; it cannot test the donor-variance transfer or *k*, which are injected or assumed.

## Software

scpPower is an R package (GPL-3) with three exported functions, scp_calibrate(), scp_simulate() and scp_power(), and YAML presets carrying per-protein parameter tables, detection parameters, *k*, cells per run and method length with their sources (Figure 1b). It is checked on every commit in the Bioconductor release container, ships an renv lockfile, and has 268 tests. Every reported number is regenerated from public data by scripts that record input hashes, seeds, software versions and wall time. The grid and the validation were run on a Slurm cluster (R 4.6.1, scp 1.22.0, limma 3.68.5).

# Results

## Calibration

The five series span the current SCP instrument classes (Table 1, Figure 2). Mean quantified proteins per cell in the endpoint matrix range from 474 (plexDIA on the Q Exactive) to 1,284 (label-free DIA on the timsTOF SCP), with missing fractions from 0.33 to 0.80. Median cell SD is 0.40 to 0.54 log2 units in the TMT series and 0.78 to 0.84 in the plexDIA series; in both TMT datasets the batch (plex) SD, 0.44 and 0.46, is as large as or larger than the cell SD, so plex layout matters at least as much as cell number. The mapped between-donor SD is 0.35 to 0.38 log2 units in every series, which is larger than a 1.25-fold change (0.32 log2). TMT detection is bimodal at the protein level (Figure 2c): most proteins are either identified in most plexes and then quantified in nearly all channels, or rarely identified, which is the signature of carrier-driven, plex-level identification. DIA detection rises smoothly with intensity (Figure 2d).

## Simulator fidelity

With per-protein resampling, all five series pass the registered fidelity gate (Figure 2b; Supporting Table S1). The largest errors are 1.2 percentage points in missing fraction per tier, 1.7% in proteins per cell, 4.8% in median cell SD and 11.3% in batch SD, against tolerances of 5 points, 10%, 20% and 25%. An earlier version of the simulator that gave every protein its tier's mean variance failed the gate on four of five series, because the mean of a right-skewed variance distribution overstates the median protein's noise by 10 to 30%; that failure and its diagnosis are retained in the repository.

## Registered design grid

Three results follow from the grid (Figure 3, Table 2). First, a 1.25-fold change is out of reach: the largest registered design, 12 samples per condition and 64 cells per sample, gives power 0.01 to 0.07 on every preset, because the mapped between-donor SD exceeds the effect and no number of cells removes the donor term. A 1.5-fold change reaches 0.20 to 0.55 at that design. Second, at a 2-fold change the limit is testability rather than detection. At 12 samples and 64 cells per sample, the ratio of unconditional power to testability is 0.93 to 1.02 on every preset, so the test finds almost every changed protein it can test; what it cannot test are the proteins that fail the three-cell minimum in at least one of 24 samples. For carrier TMT this fraction is set by plex-level identification, so testability is 0.37 (Specht) and 0.49 (Leduc) at 64 cells per sample (Figure 3c) and power is capped at 0.35 and 0.47 (Figure 3a). Adding cells per sample raises testability slowly; the number of plexes, and hence samples, is the lever. Third, carrier compression halves TMT power: at 2-fold, 12 samples and 64 cells, power falls from 0.35 to 0.29 and 0.14 (Specht) and from 0.47 to 0.40 and 0.24 (Leduc) as *k* goes from 1 to 0.75 and 0.5, while testability is unchanged (Figure 3d).

Only one registered design reaches 0.8 power with the fewest hours: plexDIA on the timsTOF SCP at a 2-fold change with 12 samples per condition and 32 cells per sample, power 0.865 at 192 nominal instrument hours (Figure 3e, circled). Label-free DIA reaches 0.80 only at 64 cells per sample (896 hours), and the carrier-TMT and Q Exactive plexDIA presets do not reach 0.8 anywhere in the grid. Realised false discovery rate is at most 0.054 wherever a design produces at least ten rejections per replicate; the few values up to 0.068 occur at zero-power designs with a handful of rejections.

## Validation on real cells

Label-free DIA passed the semi-synthetic gate on every eligible design (57 of 57 within 10 percentage points, maximum gap 5.4 points; Figure 4c). The Q Exactive plexDIA series passed on its six eligible designs, but with only 45 cells per line all observed powers are below 2%, so this agreement says nothing about useful designs (Figure 4d). The timsTOF series has 10 cells per line and is not evaluable.

The two carrier-TMT series passed the registered criterion at first attempt (67 of 72 and 66 of 69 designs) but with a systematic bias: every discrepant design had predicted above observed power, by up to 19.5 points, and the whole excess was in testability while conditional power agreed (Figure 4e). A pre-specified diagnostic traced the cause: the TMT detection parameters had been estimated from the raw reporter-ion PSM mask, whereas the endpoint evaluates the ratio-based protein matrix, which is sparser by 9.8 (Specht) and 7.1 (Leduc) percentage points because of peptide filtering upstream in the source pipelines. Under a rule fixed before the diagnostic was run, the detection parameters were re-estimated from the mask of the endpoint assay. The fidelity gate was rerun and passed, and the validation was rerun on a fresh protein split with predictions again frozen before observation: 72 of 72 and 69 of 69 designs within 10 points, maximum gaps 3.3 and 2.3 points, and mean signed gaps of −0.8 and +1.5 points among designs with power above 0.2 (Figure 4a,b). At the largest validated design (8 samples, 64 cells, 2-fold) Specht predicted 0.30 against 0.31 observed. The same diagnostic showed that real cells vary in depth far more than independent per-protein draws produce (SD of quantified proteins per cell 5.2 and 7.5 times the simulated value); this dependence is not modelled and did not affect the validation verdict.

## Worked designs

For a 2-fold change with 32 cells per sample, moving from 8 to 12 samples per condition raises power from 0.35 to 0.36 on the Leduc TMT preset (testability 0.40 and 0.39, conditional power 0.86 and 0.95, 59 and 87 hours), from 0.78 to 0.87 on the timsTOF plexDIA preset (testability 0.90 and 0.89, 128 and 192 hours), and from 0.65 to 0.67 on the label-free preset (testability 0.67 and 0.65, 299 and 448 hours). The TMT and label-free examples show the same thing from two sides: once conditional power is above 0.9, more samples cannot help, and the remaining gap to 0.8 is proteins that are never testable.

# Discussion

scpPower turns the design of a two-condition SCP study from precedent into a calculation whose inputs are visible. Three of its outputs change practice. Studies that hope to detect 1.25-fold changes across donors need far more biological samples than the twelve per condition that bound this grid, or a smaller donor variance than whole PBMC implies; neither more cells nor a better instrument substitutes. In carrier-TMT designs the number of plexes limits what can be tested, so samples and plex layout should be planned before cells per sample, and the reported testability should accompany any claimed power. Carrier compression is a design parameter with a large effect that no public dataset measures; a known-ratio two-proteome standard belongs in every carrier design, and the package accepts the measured *k*. The instrument-hour axis (Figure 3e) makes the trade-off between instrument classes explicit: at a 2-fold change the plexDIA timsTOF preset reaches 0.8 at a quarter of the hours the label-free preset needs, and the TMT presets never reach it within the grid.

The limitations are stated in the package and repeated here. The between-donor variance is borrowed from one bulk PBMC cohort; its transfer to single cells of any particular cell type is an assumption, and whole-PBMC variance includes composition differences that a sorted population would not show, so the grid is probably pessimistic for sorted cells. The compression factor is an assumption. The validation reaches observed power of about 0.35 on the TMT series and 0.6 on the label-free series because public datasets contain too few cells for larger designs; predictions above that range rest on the model. Cell-depth dependence is not modelled. The Q Exactive and timsTOF plexDIA presets could not be validated at useful power. The primary endpoint is one reasonable choice, per-sample means with limma-trend; cell-level mixed models or the scplainer linear model will give different numbers, and the package reports only the registered endpoint. Only two-condition designs are covered.

The package ships a pilot protocol that addresses the two assumptions a laboratory can test: four donors by two conditions by twenty cells, a two-proteome ratio standard in the carrier design, and blank wells in every batch, with the calibration call that turns the pilot into a preset. Future work is Bioconductor submission, additional presets as public multi-donor SCP datasets appear, and further endpoints.

# Data and software availability

scpPower version 1.0.0 is available at https://github.com/yhengaa/scpPower under GPL-3 [ZENODO DOI]. All calibration inputs are public: scpdata 1.20.0 (specht2019v3, leduc2022_pSCoPE, derks2022, brunner2022) and PRIDE PXD058573 (Walzik et al.). The repository contains the frozen mandate with its addenda, the reproduction scripts with input hashes and seeds, the full power grid, the fidelity and validation tables, and the figure script used here.

# Supporting Information

S1, fidelity tables for the five series under the final simulator and the earlier failed version. S2, semi-synthetic validation tables and plots, first attempt and rerun. S3, the TMT mask diagnostic. S4, the full registered power grid with false discovery rate, rejection counts and instrument hours. S5, worked-example manifests. S6, the pilot protocol. S7, the frozen mandate and addenda (the registered record).

# Acknowledgements

The package depends on datasets deposited by the Slavov, Mann, Walzik and Gatto laboratories and on the Bioconductor scp and scpdata infrastructure. Generative AI tools were used in this work and are disclosed as follows. The package code, tests and analysis scripts were written by GPT-6.1 Sol (OpenAI) from written directives; the directives, the review of each evidence bundle against the registered criteria, and the first draft of this manuscript were produced with Claude Fable 5.1 (Anthropic). The author defined the model, endpoints, datasets and acceptance criteria, froze them before each stage, verified every reported result against the reproduction scripts and the source publications, and edited the manuscript. No AI tool generated data, results or figures other than through the plotting script described in Methods and deposited in the repository. The author received no specific funding for this work.

**Author contributions.** Y.H. conceived the project, made all design decisions, verified all results and wrote the manuscript.

**Conflicts of interest.** None declared.

# References

1. Specht, H.; Emmott, E.; Petelski, A. A.; Huffman, R. G.; Perlman, D. H.; Serra, M.; Kharchenko, P.; Koller, A.; Slavov, N. Single-cell proteomic and transcriptomic analysis of macrophage heterogeneity using SCoPE2. *Genome Biol.* **2021**, *22*, 50. https://doi.org/10.1186/s13059-021-02267-5
2. Leduc, A.; Huffman, R. G.; Cantlon, J.; Khan, S.; Slavov, N. Exploring functional protein covariation across single cells using nPOP. *Genome Biol.* **2022**, *23*, 261. https://doi.org/10.1186/s13059-022-02817-5
3. Derks, J.; Leduc, A.; Wallmann, G.; Huffman, R. G.; Willetts, M.; Khan, S.; Specht, H.; Ralser, M.; Demichev, V.; Slavov, N. Increasing the throughput of sensitive proteomics by plexDIA. *Nat. Biotechnol.* **2023**, *41*, 50–59. https://doi.org/10.1038/s41587-022-01389-w
4. Brunner, A.-D.; Thielert, M.; Vasilopoulou, C. G.; Ammar, C.; Coscia, F.; Mund, A.; Hoerning, O. B.; Bache, N.; Apalategui, A.; Lubeck, M.; Richter, S.; Fischer, D. S.; Raether, O.; Park, M. A.; Meier, F.; Theis, F. J.; Mann, M. Ultra-high sensitivity mass spectrometry quantifies single-cell proteome changes upon perturbation. *Mol. Syst. Biol.* **2022**, *18*, e10798. https://doi.org/10.15252/msb.202110798
5. Rosenberger, F. A.; Thielert, M.; Strauss, M. T.; Schweizer, L.; Ammar, C.; Mädler, S. C.; Metousis, A.; Skowronek, P.; Wahle, M.; Madden, K.; Gote-Schniering, J.; Semenova, A.; Schiller, H. B.; Rodriguez, E.; Nordmann, T. M.; Mund, A.; Mann, M. Spatial single-cell mass spectrometry defines zonation of the hepatocyte proteome. *Nat. Methods* **2023**, *20*, 1530–1536. https://doi.org/10.1038/s41592-023-02007-6
6. Guise, A. J.; Misal, S. A.; Carson, R.; Chu, J.-H.; Boekweg, H.; Van Der Watt, D.; Welsh, N. C.; Truong, T.; Liang, Y.; Xu, S.; Benedetto, G.; Gagnon, J.; Payne, S. H.; Plowey, E. D.; Kelly, R. T. TDP-43-stratified single-cell proteomics of postmortem human spinal motor neurons reveals protein dynamics in amyotrophic lateral sclerosis. *Cell Rep.* **2024**, *43*, 113636. https://doi.org/10.1016/j.celrep.2023.113636
7. Hao, B.; Wei, J.; Xu, J.; Fu, Y.; Zou, Q.; Wang, X.; Wang, S.; Shi, Y.; Chen, Y.; Xie, K.; Huang, X.; Zhao, T.; Dong, Z.; Li, P.; Qin, G.; Amit, I.; Li, B. EasySCP unveils extensive liver zonation at single-cell proteomics resolution. *Nat. Commun.* **2026**, *17*, 7630. https://doi.org/10.1038/s41467-026-74525-8
8. Schmid, K. T.; Höllbacher, B.; Cruceanu, C.; Böttcher, A.; Lickert, H.; Binder, E. B.; Theis, F. J.; Heinig, M. scPower accelerates and optimizes the design of multi-sample single cell transcriptomic studies. *Nat. Commun.* **2021**, *12*, 6625. https://doi.org/10.1038/s41467-021-26779-7
9. Boekweg, H.; Guise, A. J.; Plowey, E. D.; Kelly, R. T.; Payne, S. H. Calculating sample size requirements for temporal dynamics in single-cell proteomics. *Mol. Cell. Proteomics* **2021**, *20*, 100085. https://doi.org/10.1016/j.mcpro.2021.100085
10. MacLean, D. peppwR: Estimate best fit distributions and do power analysis for hypothesis tests on phosphopeptide abundance data. R package version 0.1.0, CRAN, **2026**. https://doi.org/10.32614/CRAN.package.peppwR
11. Vanderaa, C.; Gatto, L. Replication of single-cell proteomics data reveals important computational challenges. *Expert Rev. Proteomics* **2021**, *18*, 835–843. https://doi.org/10.1080/14789450.2021.1988571
12. Vanderaa, C.; Gatto, L. The current state of single-cell proteomics data analysis. *Curr. Protoc.* **2023**, *3*, e658. https://doi.org/10.1002/cpz1.658
13. Vanderaa, C.; Gatto, L. scplainer: using linear models to understand mass spectrometry-based single-cell proteomics data. *Genome Biol.* **2025**, *26*, 237. https://doi.org/10.1186/s13059-025-03713-4
14. Cheung, T. K.; Lee, C.-Y.; Bayer, F. P.; McCoy, A.; Kuster, B.; Rose, C. M. Defining the carrier proteome limit for single-cell proteomics. *Nat. Methods* **2021**, *18*, 76–83. https://doi.org/10.1038/s41592-020-01002-5
15. Walzik, D.; Joisten, N.; Metcalfe, A. J.; Proschinger, S.; Schenk, A.; Wenzel, C.; et al. Acute exercise rewires the proteomic landscape of human immune cells. *Nat. Commun.* **2026**, *17*, 130. https://doi.org/10.1038/s41467-025-68101-9
16. Ritchie, M. E.; Phipson, B.; Wu, D.; Hu, Y.; Law, C. W.; Shi, W.; Smyth, G. K. limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Res.* **2015**, *43*, e47. https://doi.org/10.1093/nar/gkv007

# Tables

**Table 1. Calibrated series.** Cells, cell lines and source assay from the calibration record; proteins per cell and missing fraction describe the endpoint protein matrix; SDs are medians of per-protein estimates in log2 units; mapped donor SD covers uniquely mapped proteins only; cells per run and method minutes are the preset values with their source papers.

| Series | Instrument | Labelling | Cells | Lines | Proteins per cell | Missing fraction | Cell SD | Batch SD | Donor SD (n mapped) | Cells per run | Minutes per run |
| ----------- | ------------- | ------------ | ------: | -----: | --------: | --------: | -----: | ----------: | ----------: | -----: | --------: |
| Specht 2021 | Q Exactive | carrier TMT | 1,490 | 2 | 624 | 0.80 | 0.40 | 0.44 | 0.38 (1,820) | 12 | 95 |
| Leduc 2022 | Q Exactive | carrier TMT (TMTpro 18) | 1,556 | 2 | 1,004 | 0.71 | 0.54 | 0.46 | 0.36 (2,545) | 14 | 95 |
| Derks 2023 | timsTOF SCP | plexDIA | 30 | 3 | 1,015 | 0.33 | 0.84 | 0.47 | 0.36 (776) | 3 | 30 |
| Derks 2023 | Q Exactive | plexDIA | 125 | 3 | 474 | 0.60 | 0.78 | 0.35 | 0.38 (600) | 3 | 53 |
| Brunner 2022 | timsTOF SCP | label-free DIA | 434 | 1 | 1,284 | 0.48 | 0.51 | absorbed | 0.35 (1,539) | 1 | 35 |

**Table 2. Registered design summaries at *k* = 1.** Maximum power in the grid with its design and nominal hours, and the cheapest design reaching 0.8. TMT rows at *k* = 0.75 and 0.5 are in Supporting Table S4.

| Series | Fold change | Maximum power | Design (n × m) | Hours | Cheapest design reaching 0.8 |
| -------------- | ------: | -------: | --------- | -----: | ---------------- |
| Specht 2021 TMT | 1.25 | 0.01 | 12 × 64 | 203 | none in grid |
| | 1.5 | 0.20 | 12 × 64 | 203 | none in grid |
| | 2 | 0.35 | 12 × 64 | 203 | none in grid |
| Leduc 2022 TMT | 1.25 | 0.04 | 12 × 64 | 174 | none in grid |
| | 1.5 | 0.31 | 12 × 64 | 174 | none in grid |
| | 2 | 0.47 | 12 × 64 | 174 | none in grid |
| Derks 2023 timsTOF plexDIA | 1.25 | 0.02 | 12 × 64 | 384 | none in grid |
| | 1.5 | 0.54 | 12 × 64 | 384 | none in grid |
| | 2 | 0.93 | 12 × 64 | 384 | 12 × 32, power 0.87, 192 h |
| Derks 2023 Q Exactive plexDIA | 1.25 | 0.01 | 12 × 64 | 678 | none in grid |
| | 1.5 | 0.37 | 12 × 64 | 678 | none in grid |
| | 2 | 0.72 | 12 × 64 | 678 | none in grid |
| Brunner 2022 label-free | 1.25 | 0.07 | 12 × 64 | 896 | none in grid |
| | 1.5 | 0.55 | 12 × 64 | 896 | none in grid |
| | 2 | 0.80 | 12 × 64 | 896 | none in grid (0.796 at 12 × 64) |

```{=openxml}
<w:p><w:r><w:br w:type="page"/></w:r></w:p>
```

# Figures

![](figures/Fig1.png){width=6.5in}

**Figure 1. Model and workflow.** (a) The five-term model: protein mean, condition effect scaled by the assumed carrier compression *k*, biological-sample effect, acquisition-batch effect and cell residual, with the two missingness mechanisms. In carrier TMT a protein is identified per plex with probability *p*~g~ and then quantified per channel with probability *q*~g~; in DIA and label-free data each cell is quantified with a protein baseline rate *r*~g~ shifted by the intensity slope. Filled and empty squares are illustrative. (b) SCP data and bulk donor data feed scp_calibrate(); the preset feeds scp_simulate() and scp_power(), which simulates internally, applies the per-sample limma-trend endpoint and returns the power table with instrument hours.

```{=openxml}
<w:p><w:r><w:br w:type="page"/></w:r></w:p>
```

![](figures/Fig2.png){width=6.5in}

**Figure 2. Calibration and fidelity.** (a) Cell, batch and mapped donor SD by abundance tier for the five series (square roots of the registered tier variances; Brunner batch variance is absorbed into the cell term). (b) Missing fraction by tier, observed and simulated mean, with the registered ±5-point band. (c) Specht and Leduc per-protein plex identification *p*~g~ against channel quantification *q*~g~ (two-dimensional counts, logarithmic opacity). (d) Per-protein quantification rate *r*~g~ against mean log2 intensity for the DIA series with the registered logistic curve; the simulator anchors each protein at its own *r*~g~.

```{=openxml}
<w:p><w:r><w:br w:type="page"/></w:r></w:p>
```

![](figures/Fig3.png){width=5.5in}

**Figure 3. Registered design grid.** (a, b) Unconditional power against cells per sample (log2 axis) at 2-fold and 1.5-fold, *k* = 1, one panel per series, line type by samples per condition, target 0.8 dotted. (c) Testability at 12 samples per condition and 2-fold. (d) Carrier-TMT power at 2-fold and 12 samples per condition for *k* = 1, 0.75 and 0.5. (e) Power against nominal instrument hours (log axis) for every 2-fold grid point including TMT sensitivity values; the circled point is the minimum-hour design reaching 0.8.

```{=openxml}
<w:p><w:r><w:br w:type="page"/></w:r></w:p>
```

![](figures/Fig4.png){width=6.5in}

**Figure 4. Validation on real cells.** Predicted against observed unconditional power at counted designs, identity line and ±10-point band, fill by fold change. (a, b) Specht and Leduc after the detection parameters were re-estimated from the endpoint mask. (c) Brunner. (d) Q Exactive, axes limited to 0 to 0.05; all powers are below 2%. (e) Specht before the mask correction, same axes as (a). timsTOF is not evaluable with 10 cells per line. Designs with fewer than three cells per sample are excluded. These comparisons validate cell-level noise, missingness and the test within the power range public cells allow; they cannot validate the donor-variance transfer or *k*.

```{=openxml}
<w:p><w:r><w:br w:type="page"/></w:r></w:p>
```

![](figures/TOC.png){width=3.25in}

**Table of contents graphic.** Reduced Figure 3e with the question "How many samples and cells does a single-cell proteomics study need?".
