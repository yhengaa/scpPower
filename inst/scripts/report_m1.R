# Render the stopped M1 run without fitting or changing any calibration model.
.libPaths(c(normalizePath(".local/R-library"), .libPaths()))
source("R/scp_calibrate.R")
out <- "reports/M1"
results <- readRDS(file.path(out, "calibration-results.rds"))
status <- jsonlite::read_json(file.path(out, "execution-status.json"), simplifyVector = TRUE)
manifest <- jsonlite::read_json(file.path(out, "run-manifest.json"), simplifyVector = FALSE)
inspection <- jsonlite::read_json(file.path(out, "dataset-inspection.json"), simplifyVector = FALSE)
versions <- .scp_versions()
version_text <- paste(names(versions), versions, sep = "=", collapse = "; ")
manifest$versions <- as.list(versions)
manifest$cohort$versions <- as.list(versions)
for (n in names(manifest$datasets)) manifest$datasets[[n]]$versions <- as.list(versions)
manifest$decision <- list(status = "STOPPED", trigger_dataset = "leduc2022_pSCoPE",
    trigger = "Significantly negative quantification slope in normalized protein assay",
    no_recalibration_after_stop = TRUE, M1_closed = FALSE)
manifest$unrun_datasets <- as.list(status$remaining_datasets)
manifest$dataset_inspection <- inspection
manifest$source_downloads <- jsonlite::read_json(".local/sources/download_manifest.json")
manifest$metadata_download <- jsonlite::read_json(".local/sources/bulk-metadata-download-manifest.json")
manifest$compression_downloads <- jsonlite::read_json(".local/sources/supplement-download-manifest.json")
manifest$engineering_checks <- list(local_R_CMD_check = "Status: OK",
    tests = "55 passed, 0 failures, 0 warnings, 0 skipped", test_seeds = c(304L, 508L),
    github_commit = "b3e397b9670d6e2c175c10670c56c101d231370e",
    github_check = "success",
    github_run_url = "https://github.com/yhengaa/scpPower/actions/runs/38358072045")
manifest$interrupted_inspection <- list(status = "superseded before calibration",
    reason = "Full QFeatures serialization was replaced with compact materialization; public inputs unchanged",
    preserved_log = "scp-inspection.log")
jsonlite::write_json(manifest, file.path(out, "run-manifest.json"), pretty = TRUE,
    auto_unbox = TRUE, na = "null", digits = NA)
status$versions <- as.list(versions)
status$reason <- "Significantly negative Leduc quantification slope: stop rule triggered. Other components remain unavailable."
jsonlite::write_json(status, file.path(out, "execution-status.json"), pretty = TRUE,
    auto_unbox = TRUE, na = "null", digits = NA)

table <- utils::read.csv(file.path(out, "calibration-table.csv"), check.names = FALSE)
table <- table[table$dataset %in% names(results), , drop = FALSE]
for (n in status$remaining_datasets) {
    blank <- table[rep(1L, 5L), , drop = FALSE]
    for (col in names(blank)) blank[[col]] <- if (is.numeric(blank[[col]])) NA_real_ else NA_character_
    object <- readRDS(paste0(".local/datasets/", n, "-proteins.rds"))
    blank$dataset <- n; blank$tier <- seq_len(5L)
    blank$assay <- "proteins"
    blank$dataset_version <- as.character(utils::packageVersion("scpdata"))
    blank$cells <- ncol(object$assays$proteins$matrix)
    blank$status <- "NOT RUN: calibration stopped at Leduc"
    blank$sample_reason <- "No calibration after registered stop"
    blank$batch_reason <- "No calibration after registered stop"
    blank$compression_status <- "No calibration after registered stop"
    blank$software_versions <- version_text
    table <- rbind(table, blank)
}
utils::write.csv(table, file.path(out, "calibration-table.csv"), row.names = FALSE, na = "NA")

dataset_order <- c("specht2019v3", "leduc2022_pSCoPE", "derks2022", "brunner2022")
caption <- function() graphics::mtext("Bulk variance is a PBMC proxy; cell and batch components are NA. No model was adjusted.",
    side = 1, outer = TRUE, line = 1, cex = 0.7)
grDevices::png(file.path(out, "variance-components.png"), width = 1600, height = 900, res = 160)
graphics::par(mfrow = c(2, 2), oma = c(3, 0, 2, 0), mar = c(4, 4, 2, 1))
for (n in dataset_order) {
    d <- table[table$dataset == n, ]
    graphics::plot(c(1, 5), c(0, 0.25), type = "n", xlab = "Mean log2 intensity quintile",
        ylab = expression(sigma^2), main = n, xaxt = "n")
    graphics::axis(1, at = 1:5)
    if (any(is.finite(d$sigma2_s))) {
        graphics::lines(d$tier, d$sigma2_s, col = "#247B83", lwd = 2, type = "b", pch = 19)
        graphics::legend("bottomleft", "Between PBMC donors (proxy)", col = "#247B83", lty = 1, pch = 19, bty = "n", cex = 0.7)
        graphics::text(3, 0.05, "Cell variance: NA\nBatch variance: NA", cex = 0.8)
    } else graphics::text(3, 0.12, if (n == "specht2019v3")
        "Protein assay is processed/imputed\nNo eligible abundance tiers" else "Not run after Leduc stop", cex = 0.8)
}
caption()
graphics::mtext(version_text, side = 1, outer = TRUE, line = 2, cex = 0.45)
grDevices::dev.off()

grDevices::png(file.path(out, "missingness-fit.png"), width = 1600, height = 900, res = 160)
graphics::par(mfrow = c(2, 2), oma = c(2, 0, 2, 0), mar = c(4, 4, 3, 1))
for (n in dataset_order) {
    if (n == "leduc2022_pSCoPE") {
        object <- readRDS(paste0(".local/datasets/", n, "-proteins.rds"))
        a <- object$assays$proteins_norm2
        z <- a$matrix[!grepl("^REV__|^CON__", a$row_metadata$Leading.razor.protein), , drop = FALSE]
        md <- a$cell_metadata
        points <- do.call(rbind, lapply(split(seq_len(ncol(z)), md$SampleType), function(i) {
            d <- z[, i, drop = FALSE]
            data.frame(intensity = rowMeans(d, na.rm = TRUE), missing = rowMeans(is.na(d)))
        }))
        graphics::plot(points$intensity, points$missing, pch = 16,
            col = grDevices::adjustcolor("#526777", alpha.f = 0.18), cex = 0.35,
            xlab = "Expected log2 intensity (protein / condition mean)",
            ylab = "Fraction missing", main = "Leduc: negative slope preserved; STOP")
        r <- results[[n]]
        graphics::lines(r$curve$expected_log2, r$curve$missing_probability, col = "#B32936", lwd = 2)
    } else {
        graphics::plot.new(); graphics::title(main = n)
        graphics::text(0.5, 0.5, if (n == "specht2019v3")
            "Original missingness mask unavailable\nNo finite logistic MLE" else "Not run after Leduc stop")
    }
}
graphics::mtext("Fitted probabilities are conditional on observable expected intensities; exclusions are recorded in the manifest.",
                side = 1, outer = TRUE, line = 0.5, cex = 0.7)
graphics::mtext(version_text, side = 1, outer = TRUE, line = 1.5, cex = 0.45)
grDevices::dev.off()
grDevices::png(file.path(out, "compression-status.png"), width = 1200, height = 500, res = 160)
graphics::plot.new(); graphics::title(main = "Carrier compression preset unavailable")
graphics::text(0.5, 0.55, paste("Cheung source spreadsheets could not be downloaded through the publisher.",
    "No coefficient, carrier-ratio series, or k = 1 replacement was invented.",
    "This is a status panel, not a compression-versus-ratio estimate.", sep = "\n"), cex = 0.9)
graphics::mtext(version_text, side = 1, line = 0.5, cex = 0.45)
grDevices::dev.off()

l <- results$leduc2022_pSCoPE
cohort <- manifest$cohort
frac <- mean(l$proteins$missing_fraction)
note <- c("# scpPower M1 plausibility note", "", "13 October 2026", "",
    "**Disposition: STOPPED. M1 is not closed. The model and registered M4 criterion were not changed.**", "",
    sprintf("**Leduc pSCoPE.** `proteins_norm2` retains missing values. After removing %d rows annotated as reverse/contaminant groups, %d proteins across %d cells were examined. The missing fraction is %.2f%%. The quantification slope is %.4f (95%% interval %.4f to %.4f), so the fitted missing fraction rises with expected intensity. This is flagged as implausible for an abundance-dependent detection preset and triggers the stop rule. Protein-wise normalization has removed the original between-protein abundance scale; the quintiles describe this supplied normalized matrix. The exclusion rule and result are preserved without another fit. Scpdata documents that `proteins_processed` is imputed and batch corrected; it was not used.",
        manifest$datasets$leduc2022_pSCoPE$removed_decoys, nrow(l$proteins), l$manifest$cells,
        100*frac, l$missingness$slope, l$missingness$slope_low, l$missingness$slope_high), "",
    "**Specht SCoPE2.** The downloaded `proteins` assay is complete and contains the dimensions used for the paper's imputed analysis. The paper explicitly describes missing quantifications and nearest-neighbor imputation. Treating the complete matrix as an original missingness estimate would disagree with that evidence. This assay is reported as uncalibratable for original missingness and variance; no substitute dataset or new protein assay was introduced.", "",
    "**Derks plexDIA and Brunner label-free DIA.** Both inputs were downloaded and inspected before the stop, but neither was calibrated afterward. Derks' supplied protein table is complete; its original quantification mask cannot be certified. The source paper's high within-plex completeness is not global matrix completeness. Brunner's protein-group matrix retains missing values and the source paper describes abundance-related detection, but no fitted estimate is claimed here. No cross-instrument estimates were pooled and no bulk Derks columns were calibrated.", "",
    sprintf("**Variance and biological units.** The selected cohort is Walzik et al., PXD058573: %d healthy PBMC donors with %d baseline observations, mapped using the checksum-verified author supplement. Leduc has %d exact identifier overlaps and %d proteins with variance estimable from at least four quantified donors. Its pooled mapped log2 donor variances range from %.4f to %.4f. These are PBMC proxies, not established donor variances for the SCP cell lines. The bulk paper's reported CVs use different processing and are not directly comparable to these unimputed log2 variances. Transfer requires a matching-cell-type lab pilot and cannot be validated by M4. Verified donor/culture IDs are absent from the SCP metadata, so cell and batch components and biological sample counts are NA; technical slides and runs were not promoted to biological samples.",
        cohort$healthy_donors, cohort$baseline_observations, l$manifest$bulk_overlap,
        l$manifest$bulk_variance_overlap, min(l$tiers$sigma2_s), max(l$tiers$sigma2_s)), "",
    "**Carrier effect.** Cheung's PMID resolves to Nature Methods, not J Proteome Res. Its paper demonstrates dependence on carrier level and ion sampling. Publisher spreadsheet requests returned browser challenges or failed connections, so no defensible per-ratio default k was obtained. Compression is NA with a reason; no single-carrier coefficient, correlation-as-compression proxy, or uncompressed default was substituted. Power comparisons belong to a later directive and were not run.", "",
    "Sources: [Specht](https://doi.org/10.1186/s13059-021-02267-5), [Leduc assay documentation](https://bioconductor.org/packages/release/data/experiment/html/scpdata.html), [Derks](https://doi.org/10.1038/s41587-022-01389-w), [Brunner](https://doi.org/10.15252/msb.202110798), [Walzik](https://doi.org/10.1038/s41467-025-68101-9), [Cheung](https://doi.org/10.1038/s41592-020-01002-5).", "",
    paste("Software:", version_text), "",
    "Reproduction: public downloads and hashes are handled by `fetch_sources.py`, `fetch_bulk_metadata.py`, and `inspect_scp.R`; `calibrate_m1.R` regenerates the preserved fit and stop; this script renders the note and plots without fitting any further models.")
writeLines(note, file.path(out, "plausibility-note.md"))
writeLines(c("# M1 evidence from a stopped calibration run", "",
    "The first package CI run is green. The scientific calibration is stopped at Leduc; M1 is not closed.", "",
    "- [Calibration table](calibration-table.csv): all registered datasets and quintiles; unrun or unidentifiable entries remain NA.",
    "- [Variance components](variance-components.png)",
    "- [Missingness fit](missingness-fit.png)",
    "- [Compression status](compression-status.png): no numerical compression series available.",
    "- [Plausibility note](plausibility-note.md)",
    "- [Run manifest](run-manifest.json)",
    "- [Bulk cohort and alternatives](bulk-cohort-selection.json)",
    "- [Dataset versions, assays, and hashes](dataset-inspection.json)", "",
    "No simulations, power calculations, or M4 evaluation were performed.", "",
    paste("Software:", version_text)), file.path(out, "README.md"))
cat("Stopped-run report and plots rendered; no calibration was rerun.\n")
