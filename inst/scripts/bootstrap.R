# Run from the package root in R 4.6 or a Bioconductor release container.
dir.create(".local/R-library", recursive = TRUE, showWarnings = FALSE)
.libPaths(c(normalizePath(".local/R-library"), .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 4L,
        timeout = 1200)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
cran <- c("renv", "jsonlite", "digest", "testthat", "roxygen2", "readxl",
          "lme4")
missing <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing)
bioc <- c("scp", "scpdata", "QFeatures", "limma", "edgeR")
missing <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) BiocManager::install(missing, ask = FALSE, update = FALSE)
# scplainer is a workflow implemented in scp, not assumed to be a separate
# installed package. Its publication and scp workflow version are recorded.
writeLines(capture.output(sessionInfo()), "reports/M1/bootstrap-session.txt")
