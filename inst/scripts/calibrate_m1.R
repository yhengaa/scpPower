# Regenerate the M1 evidence from public inputs. Run fetch_sources.py,
# fetch_bulk_metadata.py, then inspect_scp.R first when the input cache is empty.
.libPaths(c(normalizePath(".local/R-library"), .libPaths()))
source("R/scp_calibrate.R")
options(timeout = 1800)
seed <- 20261013L
started <- proc.time()[["elapsed"]]
out <- "reports/M1"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
versions <- .scp_versions()
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

inspection <- jsonlite::read_json(file.path(out, "dataset-inspection.json"), simplifyVector = FALSE)
datasets <- c("specht2019v3", "leduc2022_pSCoPE", "derks2022", "brunner2022")
all_tiers <- list(); manifests <- list(); results <- list(); stopped <- FALSE
for (dataset in datasets) {
    t0 <- proc.time()[["elapsed"]]
    path <- paste0(".local/datasets/", dataset, "-proteins.rds")
    object <- readRDS(path)
    assay <- if (dataset == "leduc2022_pSCoPE") "proteins_norm2" else "proteins"
    selected <- object$assays[[assay]]
    z <- selected$matrix
    md <- selected$cell_metadata
    rd <- selected$row_metadata
    if (!ncol(md)) md <- object$global_metadata[colnames(z), , drop = FALSE]
    if (dataset == "derks2022") {
        single <- !is.na(md$Real_single_cell) & md$Real_single_cell
        z <- z[, single, drop = FALSE]; md <- md[single, , drop = FALSE]
    }
    raw_proteins <- nrow(z)
    removed_decoys <- 0L
    if (dataset == "leduc2022_pSCoPE") {
        decoy <- grepl("^REV__|^CON__", rd$Leading.razor.protein)
        removed_decoys <- sum(decoy)
        z <- z[!decoy, , drop = FALSE]; rd <- rd[!decoy, , drop = FALSE]
    }
    ids <- if (dataset == "brunner2022") rd$Protein.Ids else
        if (dataset == "leduc2022_pSCoPE") rd$Leading.razor.protein.id else rownames(z)
    if (anyDuplicated(ids)) stop("Protein identifiers are not unique; stop M1 rather than silently collapse rows.")
    rownames(z) <- ids
    is_log <- dataset != "brunner2022"
    condition_col <- if (dataset == "brunner2022") "CellCycleStage" else
        if (dataset == "derks2022") "Celltype" else "SampleType"
    batch_col <- if (dataset == "brunner2022") "Run" else
        if (dataset == "derks2022") "File.Name" else "Set"
    # No verified independent donor/culture IDs are provided by these metadata.
    # Glass slides, digest IDs, plates, and cell types are not invented samples.
    unavailable <- dataset %in% c("specht2019v3", "derks2022")
    r <- scp_calibrate(z, md, log2_input = is_log,
        sample_col = NULL, condition_col = condition_col, batch_col = batch_col,
        bulk = bulk, bulk_metadata = bulk_meta, bulk_log2_input = FALSE,
        carrier_ratio = if (dataset %in% c("specht2019v3", "leduc2022_pSCoPE")) NULL else NULL,
        seed = seed)
    r$manifest$dataset <- dataset
    r$manifest$assay <- assay
    r$manifest$scpdata_version <- as.character(utils::packageVersion("scpdata"))
    r$manifest$source_object_sha256 <- inspection[[dataset]]$object_sha256
    r$manifest$protein_input_sha256 <- input_hash(path)
    r$manifest$raw_proteins <- raw_proteins
    r$manifest$removed_decoys <- removed_decoys
    r$manifest$original_scale <- if (is_log) "normalized log2" else "linear DIA-NN normalized quantity"
    r$manifest$bulk_cohort <- cohort$accession
    r$manifest$bulk_transfer_status <- "PBMC-to-cell-line transfer not established; mapped variance is a bulk proxy"
    r$manifest$processing <- c(paste("Select earliest documented protein assay:", assay),
        "Preserve missing values; never impute or batch-correct",
        if (dataset == "leduc2022_pSCoPE") "Remove annotated reverse and contaminant protein groups" else "No protein QC filtering added",
        if (dataset == "derks2022") "Exclude all non-single-cell columns using Real_single_cell" else "Use the selected protein assay's cells",
        "Map exact unique protein identifiers; do not split ambiguous protein groups",
        "No verified biological sample IDs: sample count and within-sample variance unavailable")
    r$manifest$assay_eligible <- !unavailable
    if (unavailable) {
        reason <- "Only processed complete protein assay supplied; original quantification mask unavailable and no finite missingness MLE; cannot certify no imputation or batch correction"
        r$missingness$reason <- reason
        r$manifest$status <- "cannot calibrate original missingness or variance from this protein assay"
        r$manifest$reason <- reason
        # Preserve the inspection and mapped cohort overlap, without claiming
        # that complete processed entries are observed measurements.
        r$proteins$quantified_cells <- NA_integer_
        r$proteins$missing_fraction <- NA_real_
        r$proteins$tier <- NA_integer_
        r$proteins$mean_log2 <- NA_real_
        r$tiers$proteins <- NA_integer_
        r$tiers$sigma2_s <- NA_real_
        r$tiers$sigma2_s_proteins <- NA_integer_
        r$tiers$bulk_overlap <- NA_integer_
        r$tiers$missing_fraction <- NA_real_
        r$curve <- data.frame(expected_log2 = numeric(), missing_probability = numeric())
    } else {
        r$manifest$status <- "partial calibration; biological units and applicable compression preset unavailable"
    }
    if (dataset %in% c("specht2019v3", "leduc2022_pSCoPE")) {
        r$compression <- list(carrier_ratio = NA_real_, k = NA_real_,
            source = "Cheung et al. 2021, Nature Methods, PMID33288958, source datasets unavailable through publisher download",
            status = "TMT external compression preset unavailable; not replaced with k=1")
    }
    r$tiers$dataset <- dataset
    r$tiers$assay <- assay
    r$tiers$dataset_version <- r$manifest$scpdata_version
    r$tiers$status <- r$manifest$status
    r$tiers$variance_source <- "PXD058573 baseline PBMC donor means; transfer to cell lines unverified"
    r$tiers$missingness_slope_low <- r$missingness$slope_low
    r$tiers$missingness_slope_high <- r$missingness$slope_high
    r$tiers$sample_reason <- "Verified independent biological sample IDs absent"
    r$tiers$batch_reason <- unique(r$proteins$sigma2_b_reason)[1]
    r$tiers$compression_status <- r$compression$status
    r$manifest$wall_seconds <- proc.time()[["elapsed"]] - t0
    write_csv(r$proteins, paste0(dataset, "-proteins.csv"))
    write_csv(r$curve, paste0(dataset, "-missingness-curve.csv"))
    all_tiers[[dataset]] <- r$tiers
    manifests[[dataset]] <- r$manifest
    results[[dataset]] <- r
    jsonlite::write_json(list(seed = seed, versions = versions, cohort = cohort,
        datasets = manifests, wall_seconds = proc.time()[["elapsed"]] - started),
        file.path(out, "run-manifest.json"), pretty = TRUE, auto_unbox = TRUE,
        na = "null", digits = NA)
    # Preserve an implausible negative quantification slope and stop.
    if (!unavailable && is.finite(r$missingness$slope_high) && r$missingness$slope_high < 0) {
        stopped <- TRUE
        cat("STOP: quantification probability decreases with expected intensity in",dataset,"\n")
        break
    }
}
table <- do.call(rbind, all_tiers)
write_csv(table, "calibration-table.csv")
saveRDS(results, file.path(out, "calibration-results.rds"))
write_csv(results[[1]]$bulk_proteins, "bulk-protein-variance.csv")
jsonlite::write_json(list(stopped_for_implausibility = stopped,
    completed_datasets = names(results), remaining_datasets = setdiff(datasets,names(results)),
    required_components_unavailable = TRUE,
    M1_complete = FALSE, reason = "Unavailable eligible protein assays, biological sample IDs, and externally sourced compression presets",
    seed = seed, versions = versions, wall_seconds = proc.time()[["elapsed"]]-started),
    file.path(out,"execution-status.json"),pretty=TRUE,auto_unbox=TRUE,na="null",digits=NA)
cat("M1 evidence saved.\n")
