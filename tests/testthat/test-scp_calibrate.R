fixture <- function() {
    set.seed(304)
    z <- matrix(rnorm(12 * 80), 12, 80,
        dimnames = list(paste0("P", 1:12), paste0("C", 1:80)))
    z <- z + seq_len(12)
    meta <- data.frame(sample = rep(paste0("S", 1:8), each = 10),
        condition = rep(rep(c("A", "B"), each = 4), each = 10),
        batch = rep(rep(paste0("B", 1:5), each = 2), 8),
        row.names = colnames(z))
    z[sample(length(z), 120)] <- NA_real_
    list(x = z, meta = meta)
}

test_that("calibration retains real units and provenance without imputing", {
    f <- fixture()
    original <- f$x
    random_state <- .Random.seed
    r <- scp_calibrate(f$x, f$meta, sample_col = "sample", condition_col = "condition")
    expect_identical(f$x, original)
    expect_identical(.Random.seed, random_state)
    expect_equal(sum(r$tiers$proteins), nrow(f$x))
    expect_equal(r$tiers$cells, rep(ncol(f$x), 5))
    expect_equal(r$tiers$samples, rep(8, 5))
    expect_true(all(is.finite(r$proteins$sigma2_e)))
    expect_true(all(is.na(r$proteins$sigma2_s)))
    expect_equal(r$manifest$seed, 20261013L)
    expect_true(all(c("R", "scp", "limma", "edgeR", "scplainer") %in%
                    names(r$manifest$versions)))
    expect_true(all(nchar(r$manifest$input_hashes) == 64L))
    expect_gte(r$manifest$wall_seconds, 0)
})

test_that("independent bulk donor means map by protein identifier", {
    f <- fixture()
    ids <- c("P2", "P1", "external")
    bulk <- matrix(rep(c(-1.5, -0.5, 0.5, 1.5), each = 3), 3, 4,
        dimnames = list(ids, paste0("D", 1:4)))
    md <- data.frame(donor = colnames(bulk), row.names = colnames(bulk))
    r <- scp_calibrate(f$x, f$meta, bulk = bulk, bulk_metadata = md)
    expect_equal(r$manifest$bulk_overlap, 2)
    expect_equal(r$proteins$sigma2_s[1:2], rep(5/3, 2))
    expect_true(all(is.na(r$proteins$sigma2_s[-(1:2)])))
    expect_equal(sum(r$tiers$bulk_overlap), 2)
    expect_error(scp_calibrate(f$x, f$meta, bulk = bulk[, 1:3],
        bulk_metadata = md[1:3, , drop = FALSE]), "fewer than four")
    replicated <- bulk[, rep(1:4, each = 2)]
    colnames(replicated) <- paste0("O", 1:8)
    repmd <- data.frame(donor = rep(paste0("D", 1:4), each = 2),
                        row.names = colnames(replicated))
    r2 <- scp_calibrate(f$x, f$meta, bulk = replicated, bulk_metadata = repmd)
    expect_equal(r2$proteins$sigma2_s, r$proteins$sigma2_s)
    expect_equal(r2$proteins$bulk_donors[1:2], c(4, 4))
})

test_that("unidentified sample and aliased batches are explicit NAs", {
    f <- fixture()
    r <- scp_calibrate(f$x, f$meta, batch_col = "batch")
    expect_true(all(is.na(r$proteins$sigma2_e)))
    expect_match(r$proteins$sigma2_e_reason[1], "sample IDs unavailable")
    f$meta$aliased <- f$meta$sample
    r <- scp_calibrate(f$x, f$meta, sample_col = "sample", batch_col = "aliased")
    expect_true(all(is.na(r$proteins$sigma2_b)))
    expect_match(r$proteins$sigma2_b_reason[1], "aliased with biological sample")
    f$meta$run <- rownames(f$meta)
    r <- scp_calibrate(f$x, f$meta, sample_col = "sample", batch_col = "run")
    expect_match(r$proteins$sigma2_b_reason[1], "One cell per run")
})

test_that("crossed batches are fitted with the specified random effects", {
    f <- fixture()
    f$x <- f$x + rep(rep(c(-2, -1, 0, 1, 2), each = 2), 8)
    r <- scp_calibrate(f$x, f$meta, sample_col = "sample",
                      condition_col = "condition", batch_col = "batch")
    expect_true(all(is.finite(r$proteins$sigma2_b)))
    expect_true(all(r$proteins$sigma2_b >= 0))
    expect_true(all(is.finite(r$proteins$sigma2_e)))
})

test_that("logistic calibration recovers intensity dependence and missing rates", {
    set.seed(508)
    mu <- seq(1, 10, length.out = 150)
    z <- matrix(rnorm(150*150, sd = 0.05), 150, 150) + mu
    p <- plogis(-3 + 0.6 * mu)
    observed <- matrix(runif(length(z)), 150, 150) < p
    z[!observed] <- NA_real_
    dimnames(z) <- list(paste0("P", 1:150), paste0("C", 1:150))
    r <- scp_calibrate(z, data.frame(row.names = colnames(z)))
    expect_equal(r$missingness$slope, 0.6, tolerance = 0.06)
    expect_equal(r$missingness$uncentered_intercept, -3, tolerance = 0.2)
    expect_equal(r$missingness$mode, "logistic")
    expect_true(all(diff(r$curve$missing_probability) < 0))
    expect_equal(mean(r$proteins$missing_fraction), mean(!observed))
})

test_that("all-observed proteins and constant expectations are not forced into an MLE", {
    f <- fixture()
    z <- f$x
    z[is.na(z)] <- 0
    r <- scp_calibrate(z, f$meta)
    expect_identical(r$missingness$mode, "unavailable")
    expect_true(is.na(r$missingness$slope))
    z <- matrix(2, 10, 10, dimnames = list(paste0("P", 1:10), paste0("C", 1:10)))
    z[, 1] <- NA_real_
    r <- scp_calibrate(z, data.frame(row.names = colnames(z)))
    expect_identical(r$missingness$mode, "MCAR")
    expect_equal(r$curve$missing_probability, 0.1)
    expect_equal(sum(r$tiers$proteins), 10)
    expect_equal(sum(r$tiers$proteins > 0), 1)
    z[1, ] <- NA_real_
    r <- scp_calibrate(z, data.frame(row.names = colnames(z)))
    expect_true(is.na(r$proteins$tier[1]))
    expect_equal(r$manifest$never_quantified_proteins, 1)
})

test_that("compression is an externally cited preset and missing ratios remain NA", {
    f <- fixture()
    presets <- data.frame(carrier_ratio = 100, k = 0.7, source = "test-only reference")
    r <- scp_calibrate(f$x, f$meta, carrier_ratio = 100, compression_presets = presets)
    expect_equal(r$compression$k, 0.7)
    expect_identical(r$compression$status, "external preset")
    r <- scp_calibrate(f$x, f$meta, carrier_ratio = 200, compression_presets = presets)
    expect_true(is.na(r$compression$k))
    expect_match(r$compression$status, "unavailable")
    expect_error(scp_calibrate(f$x, f$meta, carrier_ratio = -1), "positive")
})

test_that("software mismatches warn and return a calibration", {
    f <- fixture()
    expect_warning(r <- scp_calibrate(f$x, f$meta, reference_versions = c(R = "0.0")),
                   "mismatch.*continuing")
    expect_type(r$tiers, "list")
})

test_that("metadata alignment and input scale are checked", {
    f <- fixture()
    a <- scp_calibrate(f$x, f$meta)
    b <- scp_calibrate(f$x, f$meta[nrow(f$meta):1, , drop = FALSE])
    expect_identical(a$proteins, b$proteins)
    expect_error(scp_calibrate(f$x, f$meta[-1, ]), "row count")
    invalid <- f$meta
    rownames(invalid)[1] <- "unknown-cell"
    expect_error(scp_calibrate(f$x, invalid), "IDs do not match")
    expect_error(scp_calibrate(f$x, f$meta, sample_col = "absent"), "not found")
    linear <- 2^f$x
    linear[1, 1] <- 0
    r <- scp_calibrate(linear, f$meta, log2_input = FALSE)
    expect_equal(r$proteins$quantified_cells[1], sum(!is.na(f$x[1, -1])))
    linear[1, 1] <- -1
    expect_error(scp_calibrate(linear, f$meta, log2_input = FALSE), "negative")
})

test_that("QFeatures requires an explicit eligible assay and preserves metadata", {
    f <- fixture()
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(intensity = f$x),
        colData = f$meta)
    qf <- QFeatures::QFeatures(list(proteins_raw = se))
    expect_error(scp_calibrate(qf), "explicitly")
    r <- scp_calibrate(qf, assay = "proteins_raw", sample_col = "sample")
    m <- scp_calibrate(f$x, f$meta, sample_col = "sample")
    expect_identical(r$proteins, m$proteins)
    expect_identical(r$manifest$assay, "proteins_raw")
})

 test_that("culture components identify noise without substituting donors", {
    f <- fixture()
    f$meta$line <- rep(c('A','B'), 40)
    r <- scp_calibrate(f$x, f$meta, sample_col='line', batch_col='batch', biological_unit='culture')
    expect_true(all(is.finite(r$proteins$sigma2_e)))
    expect_true(all(is.finite(r$proteins$sigma2_b)))
    expect_equal(r$tiers$samples, rep(2,5))
    expect_true(all(is.na(r$proteins$sigma2_s)))
    f$meta$batch <- f$meta$line
    r <- scp_calibrate(f$x, f$meta, sample_col='line', batch_col='batch', biological_unit='culture')
    expect_true(all(is.finite(r$proteins$sigma2_e)))
    expect_true(all(is.na(r$proteins$sigma2_b)))
    expect_match(r$line_components$sigma2_b_reason[1], 'aliased')
    f$meta$batch <- rownames(f$meta)
    r <- scp_calibrate(f$x, f$meta, sample_col='line', batch_col='batch', biological_unit='culture')
    expect_true(all(is.na(r$proteins$sigma2_e)))
    expect_match(r$proteins$sigma2_e_reason[1], 'zero residual')
 })

 test_that("variance normalisation never changes raw missingness predictors", {
    f <- fixture()
    a <- scp_calibrate(f$x, f$meta)
    b <- scp_calibrate(f$x, f$meta, normalize_cells=TRUE)
    expect_identical(a$missingness,b$missingness)
    expect_identical(a$proteins$mean_log2,b$proteins$mean_log2)
    expect_match(b$manifest$normalization,'variance only')
    c <- scp_calibrate(f$x,f$meta,condition_col='condition')
    expect_identical(a$missingness,c$missingness)
    expect_equal(a$compression$k,1)
    expect_match(a$compression$status,'assumption; not estimated')
 })

 test_that("quantitative features aggregate on log scale without filling missing entries", {
    z <- matrix(c(1,16,NA,4,NA,NA),3,2,dimnames=list(c('f1','f2','f3'),c('c1','c2')))
    se <- SummarizedExperiment::SummarizedExperiment(assays=list(intensity=z),
        rowData=data.frame(protein=c('P1','P1','P2')),
        colData=data.frame(row.names=colnames(z)))
    qf <- QFeatures::QFeatures(list(peptides=se))
    a <- scp_calibrate(qf,assay='peptides',log2_input=FALSE,feature_protein_col='protein')
    expect_equal(a$proteins$mean_log2[1],2)
    expect_equal(a$proteins$quantified_cells,c(2,0))
    expect_true(is.na(a$proteins$mean_log2[2]))
    expect_equal(a$manifest$never_quantified_proteins,1)
 })

 test_that("feature aggregation retains orientation for a single cell", {
    z <- matrix(c(1,16,NA),3,1,dimnames=list(c('f1','f2','f3'),'c1'))
    se <- SummarizedExperiment::SummarizedExperiment(assays=list(intensity=z),
        rowData=data.frame(protein=c('P1','P1','P2')),
        colData=data.frame(row.names='c1'))
    r <- scp_calibrate(QFeatures::QFeatures(list(peptides=se)),
        assay='peptides', log2_input=FALSE,feature_protein_col='protein')
    expect_equal(r$proteins$quantified_cells,c(1,0))
    expect_equal(r$proteins$mean_log2[1],2)
    expect_equal(r$manifest$cells,1)
 })

 test_that("amendment E separates the missingness assay from retained variance", {
    f <- fixture()
    set.seed(219)
    raw <- matrix(rnorm(12*80,sd=.1),12,80)+seq(5,16,length.out=12)
    raw[matrix(runif(length(raw)),12,80)>plogis(-8+.7*rowMeans(raw))] <- NA_real_
    dimnames(raw) <- dimnames(f$x)
    old <- scp_calibrate(f$x,f$meta,sample_col='sample')
    new <- scp_calibrate(f$x,f$meta,sample_col='sample',missingness_input=raw)
    expect_identical(old$proteins,new$proteins)
    expect_gt(new$missingness$slope,0)
    expect_identical(new$manifest$separate_missingness_assay,TRUE)
    expect_equal(new$missingness_proteins$mean_log2,unname(rowMeans(raw,na.rm=TRUE)))
 })

 test_that("amendment F identifies the combined cell and run variance", {
    z <- matrix(c(0,1,2,3, 5,4,3,2),2,4,byrow=TRUE,
        dimnames=list(c('P1','P2'),paste0('C',1:4)))
    md <- data.frame(line='HeLa',run=colnames(z),row.names=colnames(z))
    r <- scp_calibrate(z,md,sample_col='line',batch_col='run',biological_unit='culture',normalize_cells=TRUE,absorb_cell_batches=TRUE)
    normalized <- sweep(z,2,apply(z,2,median)-median(apply(z,2,median)),'-')
    expect_equal(r$proteins$sigma2_e,unname(apply(normalized,1,var)))
    expect_equal(r$proteins$sigma2_b,c(0,0))
    expect_equal(r$proteins$sigma2_b_reason,rep('absorbed into sigma2_e: one cell per run',2))
    expect_equal(r$line_components$sigma2_b,c(0,0))
    expect_true(all(is.finite(r$proteins$sigma2_e)))
 })
