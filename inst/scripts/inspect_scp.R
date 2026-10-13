.libPaths(c(normalizePath(".local/R-library"), .libPaths()))
options(timeout = 1800)
dir.create(".local/datasets", recursive = TRUE, showWarnings = FALSE)
dir.create("reports/M1", recursive = TRUE, showWarnings = FALSE)
source("R/scp_calibrate.R")
datasets <- c("specht2019v3", "leduc2022_pSCoPE", "derks2022", "brunner2022")
records <- list()
for (name in datasets) {
    start <- proc.time()[["elapsed"]]
    cat("Fetching", name, "\n"); flush.console()
    record <- list(dataset = name, seed = 20261013L, versions = .scp_versions(),
                   scpdata_version = as.character(utils::packageVersion("scpdata")))
    tryCatch({
        x <- getExportedValue("scpdata", name)()
        record$object_sha256 <- digest::digest(x, algo = "sha256")
        record$assays <- names(x)
        selected <- grep("protein", names(x), value = TRUE)
        path <- paste0(".local/datasets/", name, "-proteins.rds")
        # Retain only the protein assays and global metadata used by M1.
        # Full source objects remain in ExperimentHub's public-input cache.
        materialized <- list(global_metadata = as.data.frame(SummarizedExperiment::colData(x)),
            assays = lapply(selected, function(a) list(
                matrix = as.matrix(SummarizedExperiment::assay(x[[a]])),
                cell_metadata = as.data.frame(SummarizedExperiment::colData(x[[a]])),
                row_metadata = as.data.frame(SummarizedExperiment::rowData(x[[a]])))))
        names(materialized$assays) <- selected
        saveRDS(materialized, path, compress = "gzip")
        record$serialized_sha256 <- digest::digest(file = path, algo = "sha256")
        record$protein_assays <- lapply(selected, function(a) {
            z <- SummarizedExperiment::assay(x[[a]])
            md <- as.data.frame(SummarizedExperiment::colData(x[[a]]))
            rd <- as.data.frame(SummarizedExperiment::rowData(x[[a]]))
            cat(name, a, "dimensions", dim(z), "missing", mean(is.na(z)), "\n")
            print(names(md)); print(utils::head(md, 3)); print(utils::head(rd, 3))
            list(assay = a, dimensions = dim(z), missing_fraction = mean(is.na(z)),
                range = range(z, na.rm = TRUE), cell_metadata_columns = names(md),
                row_metadata_columns = names(rd))
        })
        write.csv(as.data.frame(SummarizedExperiment::colData(x)),
                  paste0("reports/M1/", name, "-cell-metadata.csv"))
        record$status <- "downloaded"
    }, error = function(e) {
        record$status <<- "unavailable"
        record$error <<- conditionMessage(e)
        cat(name, conditionMessage(e), "\n")
    })
    record$wall_seconds <- proc.time()[["elapsed"]] - start
    records[[name]] <- record
    jsonlite::write_json(records, "reports/M1/dataset-inspection.json",
                         auto_unbox = TRUE, pretty = TRUE, na = "null")
    rm(list = intersect("x", ls())); gc()
}
