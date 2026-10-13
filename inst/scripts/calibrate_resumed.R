# Regenerate the M1 evidence from public inputs. Run fetch_sources.py,
# fetch_bulk_metadata.py, then inspect_scp.R first when the input cache is empty.
if (dir.exists(".local/R-library")) .libPaths(c(normalizePath(".local/R-library"), .libPaths()))
source("R/scp_calibrate.R")
options(timeout = 1800)
seed <- 20261013L
started <- proc.time()[["elapsed"]]
out <- "reports/M1-resumed"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
versions <- .scp_versions()
lock <- jsonlite::read_json("renv.lock")
reference_versions <- c(R = lock$R$Version, vapply(lock$Packages[c("scp","limma","edgeR")], function(z) z$Version, character(1)))
version_text <- paste(names(versions), versions, sep = "=", collapse = "; ")
write_csv <- function(d, path) {
    d$software_versions <- rep(version_text, nrow(d))
    utils::write.csv(d, file.path(out, path), row.names = FALSE, na = "NA")
}
input_hash <- function(p) digest::digest(file = p, algo = "sha256")

# Cohort selection follows the PI's rule. The bounded search found no qualifying
# downloadable MS cohort with more donors than the Walzik PBMC study.
# Use the author supplement, not the conflicting community SDRF donor labels.
map_path <- ".local/sources/exercise-supplements/41467_2025_68101_MOESM3_ESM.xlsx"
bulk_path <- ".local/sources/PXD058573-output.tsv"
if (!file.exists(map_path) || !file.exists(bulk_path))
    stop("Chosen bulk cohort lacks downloadable protein quantities or donor mapping; stop M1.")
map <- as.data.frame(readxl::read_excel(map_path, sheet = "Sample allocation"))
names(map) <- trimws(names(map))
map[] <- lapply(map, function(z) if (is.character(z)) trimws(z) else z)
has_file <- !is.na(map$`Pride file name`) & map$`Pride file name` != "N/A"
baseline <- map[!is.na(map$Timepoint) & map$Timepoint == "Baseline" & has_file, , drop = FALSE]
bulk_table <- utils::read.delim(bulk_path, check.names = FALSE,
                               na.strings = c("NA", "NaN", ""))
quantity_cols <- grep("\\.PG.Quantity$", names(bulk_table), value = TRUE)
files <- sub("^\\[[0-9]+\\] ", "", sub("\\.PG.Quantity$", "", quantity_cols))
matched <- match(baseline$`Pride file name`, files)
if (anyNA(matched)) stop("Bulk baseline samples cannot be mapped to quantified raw files; stop M1.")
bulk <- as.matrix(bulk_table[, quantity_cols[matched], drop = FALSE])
storage.mode(bulk) <- "double"
colnames(bulk) <- baseline$`Pride file name`
bulk_ids <- bulk_table$PG.ProteinGroups
# Keep distinct identified protein groups; ambiguous groups remain unmapped.
keep <- !is.na(bulk_ids) & nzchar(bulk_ids) & !duplicated(bulk_ids)
bulk <- bulk[keep, , drop = FALSE]
rownames(bulk) <- bulk_ids[keep]
bulk_meta <- data.frame(donor = baseline$ID, row.names = colnames(bulk))
if (length(unique(bulk_meta$donor)) < 4L)
    stop("Chosen bulk cohort has fewer than four mapped healthy donors; stop M1.")
write_csv(data.frame(raw_file = colnames(bulk), donor = bulk_meta$donor,
                     condition = baseline$Condition, timepoint = baseline$Timepoint),
          "bulk-baseline-donor-map.csv")

cohort <- list(accession = "PXD058573", title = "Acute exercise rewires the proteomic landscape of human immune cells",
    publication = "Walzik et al. 2026, Nature Communications, doi:10.1038/s41467-025-68101-9",
    primary_cell_population = "PBMC", healthy_donors = length(unique(bulk_meta$donor)),
    baseline_observations = ncol(bulk), input_hashes = c(protein_table = input_hash(bulk_path),
        donor_map = input_hash(map_path)),
    donor_map_exclusions = "Author supplement includes a participant with N/A proteomics files; only donors with real proteomics files belong to the selected MS cohort",
    processing = "Unimputed Spectronaut PG.Quantity; baseline visits only; log2; average within donor; unbiased variance across measured donor means; minimum four quantified donors per protein",
    alternatives = list(
        list(cohort = "Rieckmann et al. PXD004352", donors = 4, status = "PI default; fewer donors"),
        list(cohort = "In-Depth Proteomic Characterization of Classical and Non-Classical Monocyte Subsets, 2018", donors = 12,
            status = "Fewer donors; spectral-count protein quantities"),
        list(cohort = "Shchukina et al. healthy aging", status = "Excluded: proteomic assay is plasma, not primary monocytes"),
        list(cohort = "Regulatory T-cell PXD022095 and MAIT PXD052574", donors = 3,
            status = "Excluded: below four donors")),
    selection_scope = "Most donors among qualifying cohorts found in the documented bounded public search; not a claim of an exhaustive global maximum",
    transfer_assumption = "Bulk PBMC variance transfers to SCP only for the same PBMC population; matching IDs alone do not establish biological transfer to cell lines. Lab pilot required; not tested by semi-synthetic M4.",
    seed = seed, versions = versions)
jsonlite::write_json(cohort, file.path(out, "bulk-cohort-selection.json"),
                     pretty = TRUE, auto_unbox = TRUE, na = "null", digits = NA)

# Use the approved assays. Specht's supplied peptide matrix is already logged
# and row/column normalised by its source producer: retain that fact and flag it.
inspection <- jsonlite::read_json(file.path(out,'input-inspection.json'))
results <- list(); manifests <- list(); summaries <- list(); all_tiers <- list()
cv <- function(v) {v<-v[is.finite(v)]; if(length(v)<2 || mean(v)<=0) NA_real_ else sd(v)/mean(v)}
for(dataset in c('specht2019v3','leduc2022_pSCoPE','derks2022','brunner2022')) {
 object <- readRDS(paste0('.local/datasets/',dataset,'-raw.rds'))
 for(a in names(object$assays)) {
  t0 <- proc.time()[['elapsed']]
  label <- if(dataset=='derks2022') paste0(dataset,if(a=='tims_prec_extracted') '_timsTOF' else '_QExactive') else dataset
  cat('Calibrating',label,'\n');flush.console()
  s <- object$assays[[a]]; z <- s$matrix; rd<-s$row_metadata; md<-s$cell_metadata
  metadata_reordered <- 0L
  if(dataset=='derks2022') {
   # The hub metadata rows are permuted relative to their row names. The
   # unique physical file/label key restores identity, independently of values.
   full <- object$global_metadata; keys<-paste0(full$File.Name,'.',full$Label)
   if(anyDuplicated(keys)) stop('Derks physical file/label metadata key is not unique.')
   ix<-match(colnames(z),keys)
   if(anyNA(ix)) stop('Derks file/label identity cannot be restored.')
   metadata_reordered<-sum(rownames(full)[ix]!=colnames(z))
   md<-full[ix,,drop=FALSE];rownames(md)<-colnames(z)
   if(!all(paste0(md$File.Name,'.',md$Label)==colnames(z)))stop('Derks identity check failed.')
   single<-!is.na(md$Real_single_cell)&md$Real_single_cell
   z<-z[,single,drop=FALSE];md<-md[single,,drop=FALSE]
   expected<-if(a=='tims_prec_extracted')'timsTOFSCP' else 'Q-Exactive'
   if(!all(md$Instrument==expected))stop('Derks instrument identity check failed.')
  } else if(!ncol(md)) md<-object$global_metadata[colnames(z),,drop=FALSE]
  line_col<-if(dataset=='derks2022')'Celltype' else 'SampleType'
  batch_col<-if(dataset=='derks2022')'File.Name' else if(dataset=='brunner2022')'Run' else 'Set'
  raw_features<-nrow(z); removed<-0L
  ids<-switch(dataset,specht2019v3=rd$protein,leduc2022_pSCoPE=rd$Leading.razor.protein.id,derks2022=rd$Protein.Group,brunner2022=rd$Protein.Ids)
  if(dataset=='leduc2022_pSCoPE') {
   # Annotated reverse/contaminant matches are not human protein measurements.
   bad<-grepl('REV__|CON__',rd$Leading.razor.protein)
   removed<-sum(bad);z<-z[!bad,,drop=FALSE];ids<-ids[!bad]
  }
  if(anyNA(ids)||any(!nzchar(ids)))stop('Unmapped feature identifier; no invented mapping.')
  feature_log<-if(dataset=='specht2019v3') z else .scp_matrix(z,FALSE,'quantitative assay')
  protein_log<-if(dataset=='brunner2022')feature_log else .scp_aggregate(feature_log,ids)
  if(dataset=='brunner2022')rownames(protein_log)<-ids
  r<-scp_calibrate(protein_log,md,sample_col=line_col,condition_col=line_col,
    batch_col=batch_col,biological_unit='culture',normalize_cells=TRUE,
    bulk=bulk,bulk_metadata=bulk_meta,bulk_log2_input=FALSE,seed=seed,reference_versions=reference_versions)
  m<-r$manifest;m$dataset<-dataset;m$instrument<-if(dataset=='derks2022')expected else unique(as.character(md$MsInstrument))
  m$assay<-a;m$dataset_version<-as.character(packageVersion('scpdata'))
  m$source_object<-inspection[[dataset]]
  m$input_file_hash<-input_hash(paste0('.local/datasets/',dataset,'-raw.rds'))
  m$quantitative_assay_hash<-digest::digest(s$matrix,algo='sha256')
  m$protein_input_hash<-digest::digest(protein_log,algo='sha256')
  m$raw_features<-raw_features;m$removed_reverse_contaminant_features<-removed
  m$metadata_reordered_by_file_label<-metadata_reordered
  m$metadata_correction<-if(dataset=='derks2022')'Match unique File.Name + Label against actual quantitative column identifiers, replacing permuted metadata row-name alignment; verify instrument and single-cell flags' else 'Assay-local metadata where present; otherwise exact global cell identifiers'
  m$processing<-c(paste('Selected',a),if(dataset=='specht2019v3') 'Input already log2: no double-log. Source Peptides-raw.csv is cr_norm_log(t3), not raw abundance; upstream centering cannot be undone' else 'Zero linear quantities become NA; log2 positive quantities',
   if(dataset=='brunner2022')'No aggregation' else 'Per-cell median of measured log2 features per protein group; all-missing remains NA',
   'No new abundance/missingness filtering; no imputation or batch correction',
   'Per-cell median shift only for variance; missingness uses pre-shift protein means',
   'One sample per cell line, a culture and not a donor')
  m$peptides_log_comparison<-if(dataset=='leduc2022_pSCoPE')'Not used: peptides_log changes rows and cells and applies further normalisation; not merely log2(peptides)' else NULL
  m$upstream_normalization<-if(dataset=='specht2019v3')'Source row and column normalisation before import' else if(dataset=='leduc2022_pSCoPE')'Source reference-channel normalisation before import' else 'No additional upstream correction inferred from the selected extracted/raw assay'
  m$bulk_transfer_status<-'Whole PBMC proxy only; no sorted-subset proteomics available. Whole-PBMC variance includes donor composition differences; lineage-specific transfer to these lines not established.'
  # Prespecified implausibility criterion: significant negative quantification
  # slope. No refit or alternative model is used in response to this result.
  flags<-character()
  if(dataset=='specht2019v3') flags<-c(flags,'Specified peptides assay already row/column normalised and logged: raw abundance prerequisite not satisfied')
  if(is.finite(r$missingness$slope_high)&&r$missingness$slope_high<0)flags<-c(flags,'Significantly negative quantification-versus-intensity slope')
  m$status<-if(length(flags))'STOPPED DATASET: estimates preserved for diagnosis; not accepted calibration' else 'CALIBRATED; unavailable components explicitly identified'
  m$flags<-flags;m$wall_seconds<-proc.time()[['elapsed']]-t0
  # Coverage and CV comparisons use observed quantities only. Report both
  # raw linear protein CV across cells and residual-derived lognormal CV.
  protein_cv<-apply(2^protein_log,1,cv)
  medcellcv<-NA_real_; meancellcv<-NA_real_
  cellqc<-data.frame(cell=colnames(protein_log),proteins=rowSums(t(!is.na(protein_log))))
  if(dataset %in% c('specht2019v3','leduc2022_pSCoPE')) {
   lin<-2^feature_log;lin<-lin/rowMeans(lin,na.rm=TRUE)
   groups<-split(seq_len(nrow(lin)),ids)
   # Paper-compatible Specht >5 peptides; same conservative criterion is
   # reported for Leduc with its scope explicitly named, not silently equated.
   groups<-groups[lengths(groups)>5L]
   cvs<-vapply(groups,function(i)apply(lin[i,,drop=FALSE],2,function(v)if(sum(is.finite(v))>5L)cv(v) else NA_real_),numeric(ncol(lin)))
   cellqc$peptide_agreement_median_cv<-apply(cvs,1,median,na.rm=TRUE)
   cellqc$peptide_agreement_mean_cv<-rowMeans(cvs,na.rm=TRUE)
   medcellcv<-median(cellqc$peptide_agreement_median_cv,na.rm=TRUE)
   meancellcv<-median(cellqc$peptide_agreement_mean_cv,na.rm=TRUE)
  }
  r$proteins$observed_linear_cv_across_cells<-protein_cv
  r$proteins$residual_lognormal_cv<-sqrt(exp(log(2)^2*r$proteins$sigma2_e)-1)
  within_complete<-NA_real_
  if(dataset=='derks2022') {
   batches<-split(seq_len(ncol(protein_log)),md[[batch_col]])
   within_complete<-mean(vapply(batches,function(j){zz<-protein_log[,j,drop=FALSE];keep<-rowSums(!is.na(zz))>0;mean(!is.na(zz[keep,,drop=FALSE]))},numeric(1)))
  }
  summary<-data.frame(dataset=label,proteins=nrow(protein_log),cells=ncol(protein_log),samples=length(unique(md[[line_col]])),missing_fraction=mean(is.na(protein_log)),mean_proteins_per_cell=mean(colSums(!is.na(protein_log))),median_proteins_per_cell=median(colSums(!is.na(protein_log))),minimum_proteins_per_cell=min(colSums(!is.na(protein_log))),median_observed_protein_cv=median(protein_cv,na.rm=TRUE),median_cell_peptide_cv=medcellcv,median_cell_mean_peptide_cv=meancellcv,mean_within_plex_completeness=within_complete,bulk_overlap=m$bulk_overlap,bulk_variance_overlap=m$bulk_variance_overlap,quantification_slope=r$missingness$slope,slope_low=r$missingness$slope_low,slope_high=r$missingness$slope_high,sigma2_e_proteins=sum(is.finite(r$proteins$sigma2_e)),sigma2_b_proteins=sum(is.finite(r$proteins$sigma2_b)),boundary_fits=sum(r$proteins$fit_status=='boundary fit'),status=m$status)
  r$manifest<-m;r$tiers$dataset<-label;r$tiers$assay<-a;r$tiers$dataset_version<-m$dataset_version
  r$tiers$status<-m$status;r$tiers$sample_unit<-'cell-line culture, not donor'
  r$tiers$missingness_slope_low<-r$missingness$slope_low;r$tiers$missingness_slope_high<-r$missingness$slope_high
  r$tiers$missingness_raw_intercept<-r$missingness$uncentered_intercept
  r$tiers$batch_reason<-paste(unique(r$proteins$sigma2_b_reason),collapse='; ')
  r$tiers$cell_reason<-paste(unique(r$proteins$sigma2_e_reason),collapse='; ')
  r$tiers$compression_status<-r$compression$status
  write_csv(r$proteins,paste0(label,'-proteins.csv'));write_csv(r$curve,paste0(label,'-missingness-curve.csv'))
  write_csv(r$line_components,paste0(label,'-line-components.csv'));write_csv(cellqc,paste0(label,'-cell-qc.csv'))
  saveRDS(list(matrix=protein_log,cell_metadata=md),paste0('.local/datasets/',label,'-aggregated.rds'))
  results[[label]]<-r;manifests[[label]]<-m;summaries[[label]]<-summary;all_tiers[[label]]<-r$tiers
  cat(label,'missing',summary$missing_fraction,'mean proteins/cell',summary$mean_proteins_per_cell,'slope',summary$quantification_slope,m$status,'\n');flush.console()
 }
}
write_csv(do.call(rbind,all_tiers),'calibration-table.csv');write_csv(do.call(rbind,summaries),'dataset-summary.csv')
write_csv(results[[1]]$bulk_proteins,'bulk-protein-variance.csv')
saveRDS(results,file.path(out,'calibration-results.rds'))
jsonlite::write_json(list(mandate='1.1 with frozen addendum 1.1.1',seed=seed,versions=as.list(versions),cohort=cohort,datasets=manifests,wall_seconds=proc.time()[['elapsed']]-started,compression=list(k=1,status='registered assumption; not estimated',sensitivity=c(0.75,0.5),pilot='Two-proteome ratio standard required to estimate k'),prior_stop='Preserved in reports/M1; not overwritten',run_stopped=!any(vapply(manifests,function(m)!length(m$flags),logical(1)))),file.path(out,'run-manifest.json'),pretty=TRUE,auto_unbox=TRUE,na='null',digits=NA)
cat('Every registered dataset and both Derks instruments evaluated.\n')
