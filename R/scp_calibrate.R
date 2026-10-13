#' Calibrate variance and quantification probabilities for SCP data
#'
#' Fits the calibration stage of frozen mandate version 1.1 with addendum 1.1.1. No values are
#' imputed. Variances use log2 units and abundance tiers are quintiles of each
#' protein's mean measured log2 intensity; ties remain in the same tier.
#'
#' @param x A numeric protein-by-cell matrix or a QFeatures object.
#' @param cell_metadata Data frame with one row per cell. Named rows must match
#'   matrix columns. Required for matrix input; defaults to assay colData for
#'   QFeatures input.
#' @param assay Quantitative assay name for QFeatures input. Selection is explicit;
#'   the caller must verify no imputation or batch correction was applied.
#' @param log2_input Whether values are already on the log2 scale. For linear
#'   input, zeros are recorded as missing and positive values are log2 transformed.
#' @param sample_col Column containing biological sample IDs. With
#'   biological_unit='culture', these are cell lines, not donors.
#' @param condition_col Column specifying condition or known cell population.
#'   NULL uses an intercept-only abundance model.
#' @param batch_col Column specifying acquisition plex or run. NULL records
#'   batch variance as unavailable.
#' @param bulk Numeric protein-by-observation bulk matrix with protein IDs as
#'   row names. Replicate observations are averaged per donor on the log2 scale.
#' @param bulk_metadata Data frame aligned to bulk columns, with a `donor`
#'   column. At least four distinct donors are required for a supplied cohort.
#' @param bulk_log2_input Whether bulk quantities are already log2 transformed.
#' @param protein_ids SCP protein identifiers, defaulting to matrix row names.
#'   Only exact one-to-one identifiers are mapped; ambiguous groups receive NA.
#' @param carrier_ratio Carrier-to-single-cell amount ratio for TMT, or NULL
#'   for a dataset without a carrier.
#' @param biological_unit Either 'sample' for independent biological sample
#'   random effects, or 'culture' for the registered cell-line fixed effects.
#' @param normalize_cells Apply a per-cell median shift only for variance
#'   estimation. Missingness and tiers always use the unnormalised protein input.
#' @param feature_protein_col Optional rowData column mapping quantitative
#'   peptide or precursor rows to protein identifiers for QFeatures input.
#'   Features are log2 transformed first and aggregated by per-cell median.
#' @param compression_k Registered compression assumption, default 1. It is
#'   not estimated; future sensitivity assumptions are 0.75 and 0.5.
#' @param compression_presets Legacy optional data frame with `carrier_ratio`, `k`, and `source`
#'   columns. Values must have an explicit external source. No coefficient is
#'   estimated from a single-carrier SCP matrix or invented for an absent ratio.
#' @param seed Integer seed recorded in the manifest. Calibration is deterministic
#'   and does not change the caller's random-number state.
#' @param reference_versions Named character vector of reference software
#'   versions. Differences produce warnings and never prevent execution.
#'
#' @details The missingness predictor is the protein's mean observed log2
#'   intensity across all quantified cells in the unnormalised input. Grouped
#'   binomial maximum likelihood fits quantified versus missing counts per
#'   protein, without protein random intercepts. A common predictor offset
#'   improves numerical conditioning; it does not centre individual proteins.
#'   Both the offset and raw-intensity intercept are returned. If the slope's
#'   Wald 95 percent interval includes zero, the fitted offset intercept is
#'   retained as an MCAR rate. Never-quantified proteins are excluded and counted.
#'
#'   Culture calibration fits cell line as a fixed effect and shared batch
#'   intercepts by REML, preserving the registered hierarchical model. The
#'   cell component is the residual variance conditional on line and batch;
#'   the batch component is the variance of latent batch means after line
#'   adjustment. Batches containing only one line are excluded from batch
#'   fitting under the explicit aliasing rule. If no mixed batches identify
#'   the component, within-line/batch residual sums of squares identify cell
#'   variance only when positive residual degrees of freedom exist. One cell
#'   per run has no such degrees of freedom and cannot separate the components.
#'   Per-line component status is returned separately. Culture samples are
#'   not independent donors and do not identify between-donor variance.
#'
#'   Independent-sample calibration retains sample and batch random effects,
#'   a fixed condition, and explicit checks for aliasing. Boundary and
#'   convergence results are preserved without adjustment to source papers.
#'
#'   Bulk variance is the unbiased variance of measured log2 donor means,
#'   requiring at least four quantified donors per protein. It is mapped by
#'   identifier before being pooled within SCP abundance tiers. Bulk-to-SCP
#'   transfer is an assumption requiring a matching cell type and a lab pilot;
#'   the semi-synthetic M4 protocol cannot validate that assumption.
#'
#' @return A list containing `proteins`, `tiers`, `missingness`, `curve`,
#'   `compression`, `bulk_proteins`, and `manifest`. The manifest records versions
#'   of R, scp, limma, edgeR, and the scplainer workflow, input hashes, seed,
#'   and elapsed wall time. Every unavailable variance has a reason column.
#' @references Cheung et al. (2021). Defining the carrier proteome limit for
#'   single-cell proteomics. Nature Methods 18, 76-83. PMID: 33288958.
#'   doi:10.1038/s41592-020-01002-5. The mandate's journal label is corrected
#'   here using the verified PMID, without changing its approved text.
#' @examples
#' z <- matrix(seq_len(40) / 10, nrow = 4,
#'             dimnames = list(paste0("P", 1:4), paste0("C", 1:10)))
#' md <- data.frame(sample = rep(c("S1", "S2"), each = 5),
#'                  row.names = colnames(z))
#' calibration <- scp_calibrate(z, md, sample_col = "sample")
#' calibration$tiers
#' @export
scp_calibrate <- function(x, cell_metadata = NULL, assay = NULL,
                          log2_input = TRUE, sample_col = NULL,
                          condition_col = NULL, batch_col = NULL,
                          bulk = NULL, bulk_metadata = NULL,
                          bulk_log2_input = TRUE, protein_ids = NULL,
                          carrier_ratio = NULL, compression_presets = NULL,
                          seed = 20261013L, reference_versions = NULL,
                          biological_unit = c("sample", "culture"),
                          normalize_cells = FALSE, feature_protein_col = NULL,
                          compression_k = 1) {
    biological_unit <- match.arg(biological_unit)
    feature_ids <- NULL
    started <- proc.time()[["elapsed"]]
    versions <- .scp_versions()
    if (!is.null(reference_versions)) {
        common <- intersect(names(versions), names(reference_versions))
        changed <- common[versions[common] != reference_versions[common]]
        if (length(changed)) warning("Software version mismatch: ",
            paste(changed, collapse = ", "), "; continuing.", call. = FALSE)
    }
    if (length(seed) != 1L || !is.numeric(seed) || !is.finite(seed))
        stop("seed must be one finite integer.")
    if (inherits(x, "QFeatures")) {
        if (!requireNamespace("SummarizedExperiment", quietly = TRUE))
            stop("SummarizedExperiment is needed to read QFeatures.")
        if (is.null(assay) || length(assay) != 1L || !assay %in% names(x))
            stop("Choose one existing quantitative assay explicitly.")
        selected <- x[[assay]]
        if (!is.null(feature_protein_col)) {
            rd <- as.data.frame(SummarizedExperiment::rowData(selected))
            if (!feature_protein_col %in% names(rd)) stop("Feature protein column not found.")
            feature_ids <- as.character(rd[[feature_protein_col]])
        }
        if (is.null(cell_metadata)) cell_metadata <-
            as.data.frame(SummarizedExperiment::colData(selected))
        x <- as.matrix(SummarizedExperiment::assay(selected))
    }
    x <- .scp_matrix(x, log2_input, "x")
    feature_hash <- digest::digest(x, algo = "sha256")
    if (!is.null(feature_ids)) x <- .scp_aggregate(x, feature_ids)
    cell_metadata <- .scp_metadata(cell_metadata, colnames(x), "cell_metadata")
    for (field in c(sample_col, condition_col, batch_col)) {
        if (!field %in% names(cell_metadata))
            stop("Metadata column not found: ", field)
        if (anyNA(cell_metadata[[field]]) ||
            any(!nzchar(as.character(cell_metadata[[field]]))))
            stop("Metadata column contains missing IDs: ", field)
    }
    if (is.null(protein_ids)) protein_ids <- rownames(x)
    if (is.null(protein_ids) || length(protein_ids) != nrow(x) ||
        anyNA(protein_ids) || any(!nzchar(protein_ids)))
        stop("Provide one nonmissing protein identifier per row.")
    if (anyDuplicated(protein_ids)) stop("SCP protein identifiers must be unique.")
    hashes <- c(matrix = digest::digest(x, algo = "sha256"),
                cell_metadata = digest::digest(cell_metadata, algo = "sha256"))
    means <- rowMeans(x, na.rm = TRUE)
    means[!is.finite(means)] <- NA_real_
    tier <- .scp_tiers(means)
    conditions <- if (is.null(condition_col)) rep("all", ncol(x)) else
        as.character(cell_metadata[[condition_col]])
    curve <- .scp_missingness(x, conditions)
    variance_x <- x
    shifts <- rep(0, ncol(x))
    if (isTRUE(normalize_cells)) {
        medians <- apply(x, 2L, stats::median, na.rm = TRUE)
        if (any(!is.finite(medians))) stop("Cannot median-normalise an entirely missing cell.")
        shifts <- medians - stats::median(medians)
        variance_x <- sweep(x, 2L, shifts, "-")
    }
    fit <- lapply(seq_len(nrow(x)), function(g) {
        if (biological_unit == "culture")
            .scp_culture_variance(variance_x[g, ], cell_metadata, sample_col, batch_col)
        else .scp_variance(variance_x[g, ], cell_metadata, sample_col, condition_col, batch_col)
    })
    variances <- do.call(rbind, fit)
    proteins <- data.frame(protein = protein_ids, tier = tier,
        mean_log2 = means, quantified_cells = rowSums(!is.na(x)),
        missing_fraction = rowMeans(is.na(x)), variances,
        stringsAsFactors = FALSE)
    bulk_proteins <- data.frame(protein = character(), sigma2_s = numeric(),
                               quantified_donors = integer())
    if (!is.null(bulk)) {
        bulk <- .scp_matrix(bulk, bulk_log2_input, "bulk")
        bulk_metadata <- .scp_metadata(bulk_metadata, colnames(bulk), "bulk_metadata")
        if (!"donor" %in% names(bulk_metadata) || anyNA(bulk_metadata$donor) ||
            any(!nzchar(as.character(bulk_metadata$donor))))
            stop("bulk_metadata requires nonmissing donor IDs.")
        if (length(unique(bulk_metadata$donor)) < 4L)
            stop("The supplied bulk cohort has fewer than four donors.")
        if (is.null(rownames(bulk)) || anyDuplicated(rownames(bulk)))
            stop("bulk requires unique protein identifiers as row names.")
        donor_means <- vapply(split(seq_len(ncol(bulk)), bulk_metadata$donor),
            function(i) rowMeans(bulk[, i, drop = FALSE], na.rm = TRUE),
            numeric(nrow(bulk)))
        donor_means[!is.finite(donor_means)] <- NA_real_
        n_donors <- rowSums(!is.na(donor_means))
        donor_var <- apply(donor_means, 1L, stats::var, na.rm = TRUE)
        donor_var[n_donors < 4L] <- NA_real_
        bulk_proteins <- data.frame(protein = rownames(bulk), sigma2_s = donor_var,
                                    quantified_donors = n_donors)
        hashes <- c(hashes, bulk = digest::digest(bulk, algo = "sha256"),
            bulk_metadata = digest::digest(bulk_metadata, algo = "sha256"))
    }
    matched <- match(protein_ids, bulk_proteins$protein)
    ambiguous <- grepl(";", protein_ids) | grepl(";", bulk_proteins$protein[matched])
    matched[ambiguous & !is.na(ambiguous)] <- NA_integer_
    proteins$sigma2_s <- bulk_proteins$sigma2_s[matched]
    proteins$bulk_donors <- bulk_proteins$quantified_donors[matched]
    proteins$sigma2_s_reason <- ifelse(is.finite(proteins$sigma2_s), "estimated",
        if (is.null(bulk)) "No bulk cohort supplied" else
            "No unique identifier overlap or fewer than four quantified donors")
    compression <- .scp_compression(carrier_ratio, compression_presets)
    if (is.null(compression_presets)) {
        if (length(compression_k) != 1L || !is.numeric(compression_k) ||
            !is.finite(compression_k) || compression_k <= 0 || compression_k > 1)
            stop("compression_k must be in (0, 1].")
        compression$k <- compression_k
        compression$status <- "registered user assumption; not estimated"
        compression$source <- "Frozen addendum 1.1.1, amendment D"
        compression$sensitivity_assumptions <- c(0.75, 0.5)
        compression$pilot <- "Include a two-proteome ratio standard to estimate k."
    }
    average <- function(v) if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
    tiers <- do.call(rbind, lapply(seq_len(5L), function(q) {
        sel <- which(tier == q)
        data.frame(tier = q, proteins = length(sel), cells = ncol(x),
            samples = if (is.null(sample_col)) NA_integer_ else
                length(unique(cell_metadata[[sample_col]])),
            sigma2_s = average(proteins$sigma2_s[sel]),
            sigma2_b = average(proteins$sigma2_b[sel]),
            sigma2_e = average(proteins$sigma2_e[sel]),
            sigma2_s_proteins = sum(is.finite(proteins$sigma2_s[sel])),
            sigma2_b_proteins = sum(is.finite(proteins$sigma2_b[sel])),
            sigma2_e_proteins = sum(is.finite(proteins$sigma2_e[sel])),
            bulk_overlap = sum(!is.na(matched[sel])),
            missing_fraction = average(proteins$missing_fraction[sel]),
            missingness_intercept = curve$parameters$intercept,
            missingness_slope = curve$parameters$slope,
            intensity_center = curve$parameters$intensity_center,
            effective_slope = curve$parameters$effective_slope,
            compression_factor = compression$k,
            stringsAsFactors = FALSE)
    }))
    manifest <- list(seed = as.integer(seed), versions = versions,
        input_hashes = hashes, wall_seconds = proc.time()[["elapsed"]] - started,
        assay = assay, log2_input = log2_input, bulk_log2_input = bulk_log2_input,
        cells = ncol(x), proteins = nrow(x), never_quantified_proteins = sum(is.na(means)),
        bulk_overlap = sum(!is.na(matched)),
        bulk_variance_overlap = sum(is.finite(proteins$sigma2_s)),
        transfer_assumption = paste("Bulk donor variance transfers only to SCP",
            "of the same cell type; requires lab-pilot validation; not validated by M4."),
        pooling = "Unweighted mean of identifiable protein variances within abundance quintiles",
        expected_intensity = "Mean observed log2 protein intensity across all cells, before per-cell normalisation",
        biological_unit = biological_unit,
        normalization = if (normalize_cells) "Per-cell median shift for variance only" else "none",
        cell_median_shifts = shifts, feature_input_hash = feature_hash,
        variance_input_hash = digest::digest(variance_x, algo = "sha256"),
        compression = compression,
        missingness_excluded_trials = curve$parameters$excluded_trials,
        sample_col = sample_col, condition_col = condition_col, batch_col = batch_col)
    line_components <- data.frame()
    if (biological_unit == "culture" && !is.null(sample_col)) {
        line_components <- do.call(rbind, lapply(unique(as.character(cell_metadata[[sample_col]])), function(line) {
            cells <- as.character(cell_metadata[[sample_col]]) == line
            if (is.null(batch_col)) mixed <- rep(FALSE, ncol(x)) else {
                counts <- tapply(as.character(cell_metadata[[sample_col]]), cell_metadata[[batch_col]], function(v) length(unique(v)))
                mixed <- counts[as.character(cell_metadata[[batch_col]])] > 1L
            }
            identifiable <- sum(cells & mixed) > 1L
            data.frame(protein = protein_ids, cell_line = line,
                sigma2_e = proteins$sigma2_e,
                sigma2_b = if (identifiable) proteins$sigma2_b else NA_real_,
                sigma2_b_reason = if (identifiable) proteins$sigma2_b_reason else
                    "Batch holds one cell line only; batch and line aliased for this line",
                quantified_cells = rowSums(!is.na(x[,cells,drop=FALSE])),
                component_scope = "Shared culture-model residual and mixed-batch variance; not donor variance")
        }))
    }
    list(proteins = proteins, tiers = tiers, missingness = curve$parameters,
         curve = curve$curve, compression = compression,
         bulk_proteins = bulk_proteins, line_components = line_components, manifest = manifest)
}

.scp_matrix <- function(x, log2_input, name) {
    if (!is.matrix(x) || !is.numeric(x) || !nrow(x) || !ncol(x))
        stop(name, " must be a nonempty numeric matrix.")
    if (any(is.infinite(x))) stop(name, " contains infinite values.")
    if (!isTRUE(log2_input)) {
        if (any(x < 0, na.rm = TRUE)) stop(name, " has negative linear quantities.")
        x[x == 0] <- NA_real_
        x <- log2(x)
    }
    x[is.nan(x)] <- NA_real_
    x
}

.scp_metadata <- function(meta, ids, name) {
    if (is.null(meta) || !is.data.frame(meta)) stop(name, " must be a data frame.")
    if (is.null(ids) || anyDuplicated(ids)) stop("Matrix columns must have unique cell IDs.")
    if (nrow(meta) != length(ids)) stop(name, " row count must match matrix columns.")
    named <- !identical(rownames(meta), as.character(seq_len(nrow(meta))))
    if (named) {
        if (!setequal(rownames(meta), ids)) stop(name, " IDs do not match matrix columns.")
        meta <- meta[ids, , drop = FALSE]
    }
    rownames(meta) <- ids
    meta
}

.scp_tiers <- function(x) {
    out <- rep(NA_integer_, length(x))
    ok <- is.finite(x)
    if (!any(ok)) return(out)
    breaks <- as.numeric(stats::quantile(x[ok], probs = seq(0, 1, length.out = 6L),
                                        names = FALSE, type = 7L))
    out[ok] <- vapply(x[ok], function(y) 1L + sum(y > breaks[2:5]), integer(1))
    out
}

.scp_versions <- function() {
    packages <- c("scp", "limma", "edgeR")
    versions <- vapply(packages, function(p)
        if (length(find.package(p, quiet = TRUE))) as.character(utils::packageVersion(p))
        else "not installed", character(1))
    c(R = as.character(getRversion()), versions,
      scplainer = paste0("workflow in scp ", versions[["scp"]], "; no separate package"))
}

.scp_missingness <- function(x, conditions) {
    d <- data.frame(expected = rowMeans(x, na.rm = TRUE),
        quantified = rowSums(!is.na(x)), missing = rowSums(is.na(x)))
    ok <- is.finite(d$expected)
    excluded <- sum(d$quantified[!ok] + d$missing[!ok])
    d <- d[ok, , drop = FALSE]
    unavailable <- function(reason) list(parameters = list(intercept = NA_real_,
        slope = NA_real_, slope_low = NA_real_, slope_high = NA_real_,
        intensity_center = NA_real_, uncentered_intercept = NA_real_,
        effective_slope = NA_real_, mode = "unavailable", reason = reason,
        excluded_trials = excluded, converged = FALSE),
        curve = data.frame(expected_log2 = numeric(), missing_probability = numeric()))
    if (!nrow(d) || !sum(d$missing) || !sum(d$quantified))
        return(unavailable("Finite logistic MLE does not exist: all quantified, all missing, or no expected intensity"))
    center <- stats::weighted.mean(d$expected, d$quantified + d$missing)
    d$expected <- d$expected - center
    if (length(unique(d$expected)) < 2L) {
        intercept <- stats::qlogis(sum(d$quantified) / sum(d$quantified + d$missing))
        return(list(parameters = list(intercept = intercept, slope = NA_real_,
            slope_low = NA_real_, slope_high = NA_real_, intensity_center = center,
            uncentered_intercept = intercept, effective_slope = 0,
            mode = "MCAR", reason = "Constant expected intensity; only intercept identifiable",
            excluded_trials = excluded, converged = TRUE),
            curve = data.frame(expected_log2 = center,
                               missing_probability = stats::plogis(-intercept))))
    }
    fit_warnings <- character()
    fit <- tryCatch(withCallingHandlers(stats::glm(
        cbind(quantified, missing) ~ expected, data = d, family = stats::binomial()),
        warning = function(w) {
            fit_warnings <<- c(fit_warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        }), error = function(e) NULL)
    if (is.null(fit) || !fit$converged || any(!is.finite(stats::coef(fit))))
        return(unavailable("Logistic maximum likelihood fit failed or did not converge"))
    co <- stats::coef(fit)
    se <- sqrt(diag(stats::vcov(fit)))[[2L]]
    ci <- co[[2L]] + c(-1, 1) * stats::qnorm(0.975) * se
    mcar <- ci[[1L]] <= 0 && ci[[2L]] >= 0
    effective <- if (mcar) 0 else co[[2L]]
    grid <- seq(min(d$expected) + center, max(d$expected) + center, length.out = 200L)
    list(parameters = list(intercept = unname(co[[1L]]), slope = unname(co[[2L]]),
        slope_low = unname(ci[[1L]]), slope_high = unname(ci[[2L]]),
        intensity_center = center, uncentered_intercept = co[[1L]] - co[[2L]] * center,
        effective_slope = effective, mode = if (mcar) "MCAR" else "logistic",
        reason = paste(unique(fit_warnings), collapse = "; "),
        excluded_trials = excluded, converged = fit$converged),
        curve = data.frame(expected_log2 = grid,
            missing_probability = stats::plogis(-(co[[1L]] + effective * (grid - center)))))
}

.scp_variance <- function(y, meta, sample_col, condition_col, batch_col) {
    out <- data.frame(sigma2_e = NA_real_, sigma2_b = NA_real_,
        sigma2_e_reason = "Biological sample IDs unavailable",
        sigma2_b_reason = "Batch IDs unavailable", fit_status = "not fit",
        stringsAsFactors = FALSE)
    if (is.null(sample_col)) {
        if (!is.null(batch_col)) out$sigma2_b_reason <-
            "Biological sample IDs unavailable; sample-batch aliasing cannot be checked"
        return(out)
    }
    ok <- !is.na(y)
    d <- data.frame(y = y[ok], sample = factor(meta[[sample_col]][ok]),
        condition = if (is.null(condition_col)) factor(rep("all", sum(ok))) else
            factor(meta[[condition_col]][ok]))
    if (nrow(d) < 3L || max(table(d$sample), 0L) < 2L) {
        out$sigma2_e_reason <- "Too few quantified cells within samples"
        return(out)
    }
    use_batch <- FALSE
    if (!is.null(batch_col)) {
        d$batch <- factor(meta[[batch_col]][ok])
        by_sample <- split(as.character(d$batch), d$sample)
        sample_alias <- all(vapply(by_sample, function(z) length(unique(z)), integer(1)) == 1L)
        by_batch <- split(as.character(d$condition), d$batch)
        condition_alias <- nlevels(d$condition) > 1L &&
            all(vapply(by_batch, function(z) length(unique(z)), integer(1)) == 1L)
        reason <- if (nlevels(d$batch) < 2L) "Only one quantified batch" else
            if (all(table(d$batch) <= 1L)) "One cell per run; batch and cell noise inseparable" else
            if (sample_alias) "Batch aliased with biological sample" else
            if (condition_alias) "Batch aliased with condition" else "estimated"
        out$sigma2_b_reason <- reason
        use_batch <- identical(reason, "estimated")
        if (!use_batch && nlevels(d$batch) > 1L) {
            out$sigma2_e_reason <- paste("Cannot separate cell and batch variance:", reason)
            return(out)
        }
    }
    fixed <- if (nlevels(d$condition) > 1L) "y ~ condition" else "y ~ 1"
    random <- character()
    if (nlevels(d$sample) > 1L) random <- "(1 | sample)"
    if (use_batch) random <- c(random, "(1 | batch)")
    caught <- character()
    fit <- tryCatch(withCallingHandlers({
        if (length(random)) lme4::lmer(stats::as.formula(paste(
            fixed, paste(random, collapse = " + "), sep = " + ")), data = d,
            REML = TRUE, control = lme4::lmerControl(check.conv.singular = "ignore"))
        else stats::lm(stats::as.formula(fixed), data = d)
    }, warning = function(w) {
        caught <<- c(caught, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = function(e) {
        caught <<- conditionMessage(e); NULL
    })
    if (is.null(fit)) {
        out$sigma2_e_reason <- paste("Fit unavailable:", paste(caught, collapse = "; "))
        if (use_batch) out$sigma2_b_reason <- out$sigma2_e_reason
        return(out)
    }
    out$sigma2_e <- stats::sigma(fit)^2
    out$sigma2_e_reason <- "estimated"
    if (use_batch) {
        vc <- as.data.frame(lme4::VarCorr(fit))
        out$sigma2_b <- vc$vcov[vc$grp == "batch" & is.na(vc$var2)]
    }
    out$fit_status <- if (length(caught)) paste(caught, collapse = "; ") else
        if (inherits(fit, "merMod") && lme4::isSingular(fit)) "boundary fit" else "converged"
    out
}

.scp_compression <- function(ratio, presets) {
    out <- list(carrier_ratio = if (is.null(ratio)) NA_real_ else ratio,
        k = NA_real_, source = NA_character_, status = if (is.null(ratio))
            "not applicable: no carrier" else "external preset unavailable; not estimated from this dataset")
    if (is.null(ratio)) return(out)
    if (length(ratio) != 1L || !is.numeric(ratio) || !is.finite(ratio) || ratio <= 0)
        stop("carrier_ratio must be one positive number.")
    if (is.null(presets)) return(out)
    if (!is.data.frame(presets) || !all(c("carrier_ratio", "k", "source") %in% names(presets)))
        stop("compression_presets requires carrier_ratio, k, source.")
    selected <- which(presets$carrier_ratio == ratio)
    if (!length(selected)) return(out)
    if (length(selected) != 1L) stop("Multiple compression presets for this carrier ratio.")
    k <- presets$k[selected]
    if (!is.numeric(k) || !is.finite(k) || k < 0 ||
        is.na(presets$source[selected]) || !nzchar(presets$source[selected]))
        stop("Compression preset must have finite nonnegative k and an explicit source.")
    out$k <- k; out$source <- presets$source[selected]; out$status <- "external preset"
    out
}

# Log-scale median aggregation preserves the original quantification mask.
.scp_aggregate <- function(x, ids) {
    if (length(ids) != nrow(x) || anyNA(ids) || any(!nzchar(ids)))
        stop("Every feature needs a nonmissing protein-group identifier.")
    groups <- split(seq_len(nrow(x)), ids)
    values <- vapply(groups, function(i) apply(x[i, , drop = FALSE], 2L,
        function(v) if (all(is.na(v))) NA_real_ else stats::median(v, na.rm = TRUE)),
        numeric(ncol(x)))
    z <- t(matrix(values, nrow = ncol(x), ncol = length(groups)))
    dimnames(z) <- list(names(groups), colnames(x))
    z
}

.scp_culture_variance <- function(y, meta, line_col, batch_col) {
    out <- data.frame(sigma2_e = NA_real_, sigma2_b = NA_real_,
        sigma2_e_reason = "Culture line IDs unavailable",
        sigma2_b_reason = "Batch IDs unavailable", fit_status = "not fit")
    if (is.null(line_col)) return(out)
    ok <- is.finite(y)
    d <- data.frame(y = y[ok], line = factor(meta[[line_col]][ok]))
    if (is.null(batch_col)) {
        if (nrow(d) <= nlevels(d$line)) {
            out$sigma2_e_reason <- "No residual degrees of freedom within cell lines"; return(out)
        }
        form <- if (nlevels(d$line) > 1L) y ~ line else y ~ 1
        out$sigma2_e <- stats::sigma(stats::lm(form, d))^2
        out$sigma2_e_reason <- "estimated"; out$fit_status <- "converged"
        return(out)
    }
    d$batch <- factor(meta[[batch_col]][ok])
    counts <- tapply(as.character(meta[[line_col]]), meta[[batch_col]], function(z) length(unique(z)))
    mixed <- counts[as.character(d$batch)] > 1L
    dm <- droplevels(d[mixed, , drop = FALSE])
    use_batch <- nrow(dm) > nlevels(dm$batch) && nlevels(dm$batch) >= 2L &&
        nrow(dm) > nlevels(dm$line) + 1L
    if (!use_batch) {
        group <- interaction(d$line, d$batch, drop = TRUE)
        df <- nrow(d) - nlevels(group)
        out$sigma2_b_reason <- if (all(table(d$batch) <= 1L))
            "One cell per run; batch and cell noise inseparable; batch holds one cell line" else
            "Batch holds one cell line only or fewer than two informative mixed batches; batch-line aliasing"
        if (df > 0L) {
            residual <- d$y - stats::ave(d$y, group, FUN = mean)
            out$sigma2_e <- sum(residual^2) / df
            out$sigma2_e_reason <- "estimated from within-line/batch residuals"
            out$fit_status <- "Within-line/batch residual fit; batch unavailable"
        } else out$sigma2_e_reason <- "No replicated cells within line/batch; zero residual degrees of freedom"
        return(out)
    }
    caught <- character()
    form <- if (nlevels(dm$line) > 1L) y ~ line + (1 | batch) else y ~ 1 + (1 | batch)
    fit <- tryCatch(withCallingHandlers(lme4::lmer(form, data = dm, REML = TRUE,
        control = lme4::lmerControl(check.conv.singular = "ignore")),
        warning = function(w) {caught <<- c(caught, conditionMessage(w)); invokeRestart("muffleWarning")}),
        error = function(e) {caught <<- conditionMessage(e); NULL})
    if (is.null(fit)) {
        out$sigma2_e_reason <- out$sigma2_b_reason <- paste("Fit failed:", paste(caught, collapse = "; "))
        return(out)
    }
    out$sigma2_e <- stats::sigma(fit)^2
    vc <- as.data.frame(lme4::VarCorr(fit))
    out$sigma2_b <- vc$vcov[vc$grp == "batch" & is.na(vc$var2)]
    out$sigma2_e_reason <- out$sigma2_b_reason <- "estimated"
    out$fit_status <- if (length(caught)) paste(caught, collapse = "; ") else
        if (lme4::isSingular(fit)) "boundary fit" else "converged"
    out
}
