#' Simulate hierarchical single-cell proteomics measurements
#'
#' Generates log2 measurements under frozen mandate addendum 1.1.3. Independent
#' biological-sample effects use the mapped bulk donor variance. Missingness is
#' drawn before observing cell residual noise. Registered presets jointly
#' resample calibrated protein parameters under amendment G.
#'
#' @param calibration A single scp_calibrate result, a list of dataset results
#'   selected by preset, or NULL when a complete preset is supplied.
#' @param n_per_condition Number of samples per condition. One integer uses
#'   equal counts for A/B; two integers allow unequal counts. A zero count is
#'   permitted only with beta=0, for the registered native-design fidelity gate.
#' @param m_cells_per_sample Cells per sample: one positive integer or a vector
#'   of positive integers in sample order (A samples followed by B samples).
#' @param beta Scalar log2 condition effect or one effect per simulated protein.
#'   Named vectors must match the returned protein identifiers.
#' @param preset Dataset preset name, YAML path, or a complete preset list with
#'   tiers, curve and optional protein_parameters, detection_mode,
#'   intensity_pool, label_free, plex_size and k.
#'   A preset overrides calibration parameters. YAML expressions are disabled.
#' @param layout balanced, confounded, or an integer sample-by-batch count
#'   matrix. Matrix row sums must equal m_cells_per_sample. Optional row and
#'   column names define sample and batch IDs. Balanced multiplexed batches
#'   contain equal numbers from each condition and require equal cell totals.
#' @param k Compression assumption in (0,1], default 1. An explicitly supplied
#'   value overrides the preset; otherwise a preset's registered k is used.
#' @param seed Nonnegative integer seed. The generator uses Mersenne-Twister,
#'   Inversion and Rejection and restores the caller's RNG kind and state.
#' @param changed_fraction Fraction of proteins selected for an injected effect,
#'   default 0.1. The selected count is round(fraction * number of proteins).
#'
#' @details Under amendment G, a simulated protein jointly resamples one row
#'   of protein_parameters within its calibrated tier: mean_log2, sigma2_s,
#'   sigma2_b, sigma2_e and detection parameters. Returned source_protein
#'   identifies the draw. Missing protein variance estimates use the tier
#'   component; separately unidentified batch variance is zero. Unmapped donor
#'   variance uses the tier proxy. A never-quantified protein has no empirical
#'   intensity; its tier mean is recorded as a fallback and its detection rate
#'   remains zero. No proteins are removed for these fallbacks.
#'
#'   TMT detection_mode="tmt_two_level" draws plex identification with p_g,
#'   then channel detection with q_g conditional on identification. Neither
#'   responds to the condition effect: pooled carrier identification makes
#'   testability depend on plex count and barely on effect size. DIA and
#'   label-free detection_mode="protein_rate" use logit(r_g) plus the registered
#'   curve slope times (expected intensity minus the resampled protein mean).
#'   Baseline rates zero and one remain exact; no clipping is applied.
#'
#'   Sample, batch and cell effects are independent normal draws with the
#'   resampled variances. The selected condition effect is k times beta.
#'   TMT baseline and detection use raw PSM intensities, while variance retains
#'   the accepted ratio scale. Custom tier-only presets remain supported for
#'   compatibility and use the earlier tier-normal and global-curve model;
#'   all registered presets use per-protein parameters.
#'
#'   Label-free presets use one acquisition batch per cell and a zero separate
#'   batch component under amendment F. Their cell variance already includes
#'   run noise. A custom layout must also have one cell per batch in that case.
#'   For balanced multiplexing with an odd capacity, one slot can remain empty.
#'
#'   The M2 native-design gate uses exact sample-by-batch counts, cell lines as
#'   samples and beta=0. This diagnostic does not establish donor transfer or
#'   change the registered M4 gate. Simulation and calibration are never tuned
#'   to make fidelity comparisons pass. Software mismatches warn and continue.
#'
#' @return List with matrix (protein by cell log2 values, NA when missing),
#'   cell_metadata (sample, condition, batch), protein_metadata (tier, mu,
#'   source_protein, resampled components and detection parameters, changed,
#'   beta_true, beta_observed), expected_log2, truth effects and manifest.
#'   The manifest records seed, RNG kind, software versions, hashes and wall time.
#' @examples
#' preset <- list(tiers = data.frame(tier = 1, proteins = 10,
#'     mean_log2 = 10, sigma2_s = 0.2, sigma2_b = 0.1, sigma2_e = 0.3),
#'     curve = list(intercept = 1, effective_slope = 0.5, intensity_center = 10),
#'     label_free = FALSE, plex_size = 4, k = 1)
#' simulated <- scp_simulate(NULL, 2, 3, beta = 1, preset = preset, seed = 19)
#' dim(simulated$matrix)
#' @export
scp_simulate <- function(calibration, n_per_condition, m_cells_per_sample,
                         beta, preset = NULL, layout = "balanced", k = 1,
                         seed = 20261013L, changed_fraction = 0.1) {
    started <- proc.time()[["elapsed"]]
    config <- .scp_simulation_config(calibration, preset)
    versions <- .scp_versions()
    if (!is.null(config$versions)) {
        ref <- unlist(config$versions)
        common <- intersect(names(ref), names(versions))
        changed <- common[ref[common] != versions[common]]
        if (length(changed)) warning("Software version mismatch: ",
            paste(changed, collapse = ", "), "; continuing.", call. = FALSE)
    }
    integer_values <- function(v, name, allow_zero = FALSE) {
        if (!is.numeric(v) || !length(v) || any(!is.finite(v)) ||
            any(v != floor(v)) || any(v < if (allow_zero) 0 else 1) ||
            any(v > .Machine$integer.max)) stop(name, " requires valid integers.")
        as.integer(v)
    }
    n <- integer_values(n_per_condition, "n_per_condition", TRUE)
    if (length(n) == 1L) n <- rep(n, 2L)
    if (length(n) != 2L || sum(n) < 1L) stop("Supply one count or two condition counts.")

    cells <- integer_values(m_cells_per_sample, "m_cells_per_sample")
    ns <- sum(n)
    if (length(cells) == 1L) cells <- rep(cells, ns)
    if (length(cells) != ns) stop("Provide one cell count per biological sample.")
    if (missing(k) && !is.null(config$k)) k <- config$k
    if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k <= 0 || k > 1)
        stop("k must be in (0, 1].")
    if (length(seed) != 1L || !is.numeric(seed) || !is.finite(seed) ||
        seed < 0 || seed != floor(seed) || seed > .Machine$integer.max)
        stop("seed must be a nonnegative integer.")
    if (!is.numeric(changed_fraction) || length(changed_fraction) != 1L ||
        !is.finite(changed_fraction) || changed_fraction < 0 || changed_fraction > 1)
        stop("changed_fraction must be in [0,1].")
    d <- config$tiers
    required <- c("tier", "proteins", "sigma2_s", "sigma2_b", "sigma2_e")
    if (!is.data.frame(d) || !all(required %in% names(d)))
        stop("Preset tiers require tier, proteins, sigma2_s, sigma2_b, sigma2_e.")
    counts <- integer_values(d$proteins, "Tier protein counts", TRUE)
    if (!sum(counts) || anyDuplicated(d$tier)) stop("Tiers must be unique with at least one protein.")
    use <- counts > 0L
    for (field in c("sigma2_s", "sigma2_b", "sigma2_e")) {
        if (!is.numeric(d[[field]]) || any(!is.finite(d[[field]][use])) ||
            any(d[[field]][use] < 0)) stop("Unidentified or invalid tier component: ", field)
    }
    curve <- config$curve
    slope <- if (!is.null(curve$effective_slope)) curve$effective_slope else curve$slope
    center <- if (is.null(curve$intensity_center)) 0 else curve$intensity_center
    if (is.null(curve$intercept) || is.null(slope) ||
        any(!is.finite(c(curve$intercept, slope, center))) ||
        length(c(curve$intercept, slope, center)) != 3L)
        stop("A finite fitted quantification curve is required.")
    if (isTRUE(config$label_free) && any(d$sigma2_b[use] != 0))
        stop("A label-free preset must absorb run noise into sigma2_e and set sigma2_b=0.")
    pools <- lapply(seq_len(nrow(d)), function(i) {
        pool <- config$intensity_pool[[as.character(d$tier[i])]]
        if (is.null(pool) && "mean_log2" %in% names(d)) pool <- d$mean_log2[i]
        pool <- as.numeric(unlist(pool))
        pool <- pool[is.finite(pool)]
        if (counts[i] > 0 && !length(pool)) stop("Tier lacks calibrated log2 means: ", d$tier[i])
        pool
    })
    parameters <- config$protein_parameters
    mode <- if (is.null(config$detection_mode)) "global_curve" else config$detection_mode
    if (!mode %in% c("global_curve", "protein_rate", "tmt_two_level"))
        stop("Unknown protein detection mode.")
    if (!is.null(parameters)) {
        fields <- c("protein", "tier", "mean_log2", "sigma2_s", "sigma2_b", "sigma2_e")
        if (!all(fields %in% names(parameters)) || anyNA(parameters$protein) ||
            anyDuplicated(parameters$protein)) stop("Invalid per-protein parameter identifiers or fields.")
        for (field in c("mean_log2", "sigma2_s", "sigma2_b", "sigma2_e")) {
            if (!is.numeric(parameters[[field]]) || any(!is.finite(parameters[[field]])) ||
                (field != "mean_log2" && any(parameters[[field]] < 0)))
                stop("Invalid per-protein component: ", field)
        }
        for (tier in d$tier[use]) {
            if (!any(parameters$tier == tier)) stop("Tier lacks calibrated proteins: ", tier)
        }
        detection_fields <- if (mode == "tmt_two_level") c("p_g", "q_g") else
            if (mode == "protein_rate") "r_g" else character()
        for (field in detection_fields) {
            x <- parameters[[field]]
            if (is.null(x) || !is.numeric(x) || any(!is.finite(x)) || any(x < 0 | x > 1))
                stop("Invalid protein detection probability: ", field)
        }
        if (isTRUE(config$label_free) && any(parameters$sigma2_b != 0))
            stop("Label-free per-protein sigma2_b must be zero.")
    } else if (mode != "global_curve") stop("Per-protein detection requires protein_parameters.")
    rows <- rep(seq_len(nrow(d)), counts)
    np <- length(rows)
    protein_ids <- config$protein_ids
    if (is.null(protein_ids)) protein_ids <- paste0("P", sprintf("%05d", seq_len(np)))
    protein_ids <- as.character(unlist(protein_ids))
    if (length(protein_ids) != np || anyNA(protein_ids) || any(!nzchar(protein_ids)) || anyDuplicated(protein_ids))
        stop("Preset protein identifiers do not match the tier counts.")
    if (!is.numeric(beta) || !length(beta) || any(!is.finite(beta)) ||
        !length(beta) %in% c(1L, np)) stop("beta must be finite and scalar or one value per protein.")
    if (length(beta) == np && !is.null(names(beta))) {
        if (!setequal(names(beta), protein_ids)) stop("Named beta identifiers do not match proteins.")
        beta <- beta[protein_ids]
    }
    beta <- rep(beta, length.out = np)
    if (any(n == 0L) && any(beta != 0)) stop("A single-condition design requires beta=0.")
    # Preserve both the RNG kind and absence/presence of .Random.seed.
    previous_kind <- RNGkind()
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) previous_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
        do.call(RNGkind, as.list(previous_kind))
        if (had_seed) assign(".Random.seed", previous_seed, envir = .GlobalEnv)
        else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(as.integer(seed))
    metadata <- .scp_simulation_layout(n, cells, layout,
        isTRUE(config$label_free), config$plex_size)
    sample_ids <- unique(metadata$sample)
    batch_ids <- unique(metadata$batch)
    nc <- nrow(metadata)
    if (!is.null(parameters)) {
        source_rows <- unlist(lapply(seq_len(nrow(d)), function(i) {
            if (!counts[i]) return(integer())
            pool <- which(parameters$tier == d$tier[i])
            pool[sample.int(length(pool), counts[i], replace = TRUE)]
        }), use.names = FALSE)
        drawn <- parameters[source_rows, , drop = FALSE]
        mu <- drawn$mean_log2
        sigma_s <- drawn$sigma2_s; sigma_b <- drawn$sigma2_b; sigma_e <- drawn$sigma2_e
    } else {
        mu <- unlist(lapply(seq_len(nrow(d)), function(i) {
            if (!counts[i]) return(numeric())
            pools[[i]][sample.int(length(pools[[i]]), counts[i], replace = TRUE)]
        }), use.names = FALSE)
        sigma_s <- d$sigma2_s[rows]; sigma_b <- d$sigma2_b[rows]; sigma_e <- d$sigma2_e[rows]
    }
    changed <- rep(FALSE, np)
    number_changed <- as.integer(round(changed_fraction * np))
    if (number_changed) changed[sample.int(np, number_changed)] <- TRUE
    effect <- beta * changed
    sample_effects <- matrix(stats::rnorm(np * length(sample_ids)), np) *
        sqrt(sigma_s)
    batch_effects <- matrix(stats::rnorm(np * length(batch_ids)), np) *
        sqrt(sigma_b)
    dimnames(sample_effects) <- list(protein_ids, sample_ids)
    dimnames(batch_effects) <- list(protein_ids, batch_ids)
    expected <- matrix(mu, np, nc) +
        sample_effects[, match(metadata$sample, sample_ids), drop = FALSE] +
        batch_effects[, match(metadata$batch, batch_ids), drop = FALSE] +
        outer(k * effect, as.integer(metadata$condition == "B"))
    value <- expected + matrix(stats::rnorm(np * nc), np) * sqrt(sigma_e)
    identified <- NULL
    if (mode == "tmt_two_level") {
        identified <- matrix(stats::runif(np * length(batch_ids)), np) < drawn$p_g
        dimnames(identified) <- list(protein_ids, batch_ids)
        detected <- matrix(stats::runif(np * nc), np) < drawn$q_g
        detected <- detected & identified[, match(metadata$batch, batch_ids), drop = FALSE]
    } else if (mode == "protein_rate") {
        probability <- stats::plogis(stats::qlogis(drawn$r_g) + slope * (expected - mu))
        detected <- matrix(stats::runif(np * nc), np) < probability
    } else {
        probability <- stats::plogis(curve$intercept + slope * (expected - center))
        detected <- matrix(stats::runif(np * nc), np) >= (1 - probability)
    }
    value[!detected] <- NA_real_
    dimnames(value) <- dimnames(expected) <- list(protein_ids, rownames(metadata))
    protein_metadata <- data.frame(protein = protein_ids, tier = d$tier[rows],
        mu = mu, targeted = changed, changed = effect != 0, beta_true = effect, beta_observed = k * effect)
    if (!is.null(parameters)) {
        protein_metadata$source_protein <- drawn$protein
        for (field in intersect(c("sigma2_s", "sigma2_b", "sigma2_e", "p_g", "q_g", "r_g"), names(drawn)))
            protein_metadata[[field]] <- drawn[[field]]
    }
    hashes <- c(calibration = digest::digest(calibration, algo = "sha256"),
        preset = digest::digest(config, algo = "sha256"),
        layout = digest::digest(layout, algo = "sha256"),
        matrix = digest::digest(value, algo = "sha256"),
        cell_metadata = digest::digest(metadata, algo = "sha256"))
    manifest <- list(seed = as.integer(seed), RNGkind = RNGkind(),
        versions = as.list(versions), input_hashes = hashes,
        n_per_condition = n, m_cells_per_sample = cells, k = k,
        changed_fraction = changed_fraction, targeted_proteins = number_changed,
        changed_proteins = sum(effect != 0),
        proteins = np, cells = nc, samples = ns, batches = length(batch_ids),
        dataset = config$dataset, carrier_ratio = config$carrier_ratio,
        layout = if (is.matrix(layout)) "custom sample-by-batch counts" else layout,
        label_free = isTRUE(config$label_free),
        baseline = if (is.null(parameters)) "Legacy tier parameters" else "Joint calibrated-protein row resampling within tier",
        missingness = mode,
        modelling_consequence = if (mode == "tmt_two_level")
            "Carrier pools both conditions: p_g and q_g do not respond to effect size; testability depends on plex count" else
            "Detection changes with expected intensity before residual noise",
        mandate_addendum = "1.1.3",
        provenance = config$provenance,
        wall_seconds = proc.time()[["elapsed"]] - started)
    list(matrix = value, cell_metadata = metadata, protein_metadata = protein_metadata,
        expected_log2 = expected, truth = list(sample_effects = sample_effects,
        batch_effects = batch_effects, plex_identified = identified), manifest = manifest)
}

.scp_simulation_config <- function(calibration, preset) {
    if (is.character(preset) && length(preset) == 1L &&
        is.list(calibration) && is.null(calibration$tiers) && preset %in% names(calibration))
        return(.scp_simulation_config(calibration[[preset]], NULL))
    if (is.character(preset)) {
        if (length(preset) != 1L) stop("Select one preset.")
        path <- preset
        if (!file.exists(path)) {
            path <- system.file("presets", paste0(preset, ".yaml"), package = "scpPower")
            if (!nzchar(path)) path <- file.path("inst", "presets", paste0(preset, ".yaml"))
        }
        if (!file.exists(path)) stop("Preset YAML not found: ", preset)
        selected <- yaml::read_yaml(path, eval.expr = FALSE)
    } else if (is.list(preset)) selected <- preset
    else if (is.null(preset) && is.list(calibration) && !is.null(calibration$tiers)) {
        selected <- list(tiers = calibration$tiers, curve = calibration$missingness,
            label_free = isTRUE(calibration$manifest$absorb_cell_batches),
            plex_size = 2L, k = calibration$compression$k,
            dataset = calibration$manifest$dataset, versions = calibration$manifest$versions,
            carrier_ratio = calibration$compression$carrier_ratio,
            provenance = calibration$manifest)
        p <- calibration$proteins
        raw_mean <- if ("missingness_mean_log2" %in% names(p)) p$missingness_mean_log2 else p$mean_log2
        if (!"missingness_mean_log2" %in% names(p) && !is.null(calibration$missingness_proteins)) {
            matched <- match(p$protein, calibration$missingness_proteins$protein)
            available <- is.finite(calibration$missingness_proteins$mean_log2[matched])
            raw_mean[available] <- calibration$missingness_proteins$mean_log2[matched][available]
        }
        tidx <- match(p$tier, selected$tiers$tier)
        table <- data.frame(protein = p$protein, tier = p$tier, mean_log2 = raw_mean,
            sigma2_s = p$sigma2_s, sigma2_b = p$sigma2_b, sigma2_e = p$sigma2_e)
        table$mean_fallback <- !is.finite(table$mean_log2)
        table$mean_log2[table$mean_fallback] <- selected$tiers$mean_log2[tidx[table$mean_fallback]]
        for (field in c("sigma2_s", "sigma2_b", "sigma2_e")) {
            missing_component <- !is.finite(table[[field]])
            table[[paste0(field, "_tier_fallback")]] <- missing_component
            table[[field]][missing_component] <- selected$tiers[[field]][tidx[missing_component]]
        }
        if (isTRUE(selected$label_free)) table$sigma2_b <- 0
        # A structurally unidentified batch effect is zero; an otherwise
        # unavailable protein estimate uses the tier estimate under G1.
        unidentified <- grepl("alias|unidentified|absorbed|one cell per", p$sigma2_b_reason, ignore.case = TRUE)
        unidentified[is.na(unidentified)] <- FALSE
        table$sigma2_b[unidentified] <- 0
        table$sigma2_b_unidentified <- unidentified
        table$r_g <- 1 - p$missing_fraction
        if (!is.null(calibration$missingness_proteins)) {
            matched <- match(p$protein, calibration$missingness_proteins$protein)
            table$r_g <- 1 - calibration$missingness_proteins$missing_fraction[matched]
        }
        selected$protein_parameters <- table[!is.na(table$tier), , drop = FALSE]
        selected$detection_mode <- "protein_rate"
        selected$intensity_pool <- stats::setNames(lapply(selected$tiers$tier,
            function(t) raw_mean[!is.na(p$tier) & p$tier == t & is.finite(raw_mean)]),
            as.character(selected$tiers$tier))
        selected$protein_ids <- unlist(lapply(selected$tiers$tier,
            function(t) as.character(p$protein[!is.na(p$tier) & p$tier == t])), use.names = FALSE)
    } else stop("Supply one calibration result or a complete preset.")
    if (!is.data.frame(selected$tiers)) {
        if (all(c("tier", "proteins", "sigma2_s", "sigma2_b", "sigma2_e") %in% names(selected$tiers)))
            selected$tiers <- as.data.frame(selected$tiers, stringsAsFactors = FALSE)
        else selected$tiers <- do.call(rbind,
            lapply(selected$tiers, function(row) as.data.frame(row, stringsAsFactors = FALSE)))
    }
    if (!is.null(selected$protein_parameters) && !is.data.frame(selected$protein_parameters)) {
        pp <- selected$protein_parameters
        if (all(c("protein", "tier", "mean_log2") %in% names(pp)))
            selected$protein_parameters <- as.data.frame(pp, stringsAsFactors = FALSE)
        else selected$protein_parameters <- do.call(rbind,
            lapply(pp, function(row) as.data.frame(row, stringsAsFactors = FALSE)))
    }
    if (is.null(selected$intensity_pool)) selected$intensity_pool <- list()
    if (is.null(selected$k)) selected$k <- 1
    if (is.null(selected$plex_size)) selected$plex_size <- 2L
    selected
}

.scp_simulation_layout <- function(n, cells, layout, label_free, capacity) {
    condition <- rep(c("A", "B"), n)
    samples <- paste0("S", seq_along(cells))
    if (is.matrix(layout)) {
        if (!is.numeric(layout) || nrow(layout) != length(cells) || !ncol(layout) ||
            any(!is.finite(layout)) || any(layout < 0) || any(layout != floor(layout)) ||
            any(rowSums(layout) != cells) || any(colSums(layout) == 0))
            stop("Custom layout must have integer counts, the required row sums, and occupied batches.")
        if (label_free && any(colSums(layout) != 1L))
            stop("Label-free custom layouts require one cell per batch.")
        if (!is.null(rownames(layout))) samples <- rownames(layout)
        batches <- if (is.null(colnames(layout))) paste0("B", seq_len(ncol(layout))) else colnames(layout)
        if (anyNA(samples) || any(!nzchar(samples)) || anyDuplicated(samples) ||
            anyNA(batches) || any(!nzchar(batches)) || anyDuplicated(batches))
            stop("Layout sample and batch identifiers must be nonmissing and unique.")
        parts <- lapply(seq_len(ncol(layout)), function(b) {
            s <- rep(seq_len(nrow(layout)), layout[, b])
            data.frame(sample = samples[s], condition = condition[s], batch = batches[b])
        })
        metadata <- do.call(rbind, parts)
    } else {
        if (!is.character(layout) || length(layout) != 1L ||
            !layout %in% c("balanced", "confounded")) stop("Unknown layout.")
        metadata <- data.frame(sample = rep(samples, cells),
            condition = rep(condition, cells), batch = NA_character_)
        if (label_free) metadata$batch <- paste0("B", seq_len(nrow(metadata)))
        else {
            if (length(capacity) != 1L || !is.numeric(capacity) || !is.finite(capacity) ||
                capacity < 2 || capacity != floor(capacity)) stop("Multiplexed capacity must be an integer >=2.")
            A <- which(metadata$condition == "A"); B <- which(metadata$condition == "B")
            if (layout == "balanced") {
                if (length(A) != length(B)) stop("Balanced layout requires equal total cells in both conditions.")
                A <- A[sample.int(length(A))]; B <- B[sample.int(length(B))]
                half <- floor(capacity / 2)
                batch <- paste0("B", ceiling(seq_along(A) / half))
                metadata$batch[A] <- metadata$batch[B] <- batch
            } else {
                b <- 0L
                for (index in list(A, B)) {
                    if (!length(index)) next
                    index <- index[sample.int(length(index))]
                    ids <- ceiling(seq_along(index) / capacity)
                    metadata$batch[index] <- paste0("B", b + ids)
                    b <- b + max(ids)
                }
            }
        }
    }
    rownames(metadata) <- paste0("C", sprintf("%06d", seq_len(nrow(metadata))))
    metadata
}
