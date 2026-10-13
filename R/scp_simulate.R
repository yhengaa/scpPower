#' Simulate hierarchical single-cell proteomics measurements
#'
#' Generates log2 measurements under frozen mandate addendum 1.1.2. Independent
#' biological-sample effects use the mapped bulk donor variance. Missingness is
#' evaluated at expected intensity, before cell residual noise is added.
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
#'   tiers, curve and optional intensity_pool, label_free, plex_size and k.
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
#' @details Protein baseline means are sampled with replacement from calibrated
#'   observed protein means within each retained tier. A tier with only its mean
#'   supplied uses a point mass at that mean. TMT presets use raw reporter means
#'   for baseline and missingness and retain ratio-based log2 variance; the two
#'   scales are recorded in preset provenance rather than interchanged.
#'
#'   For each protein, sample and batch effects are independent normal draws
#'   with the tier variances. The selected condition effect is k times beta.
#'   Independent cell residuals have the tier cell variance. Quantification
#'   probability is plogis(intercept + effective_slope * (expected - center));
#'   the missingness draw is independent of the realized residual conditional
#'   on expected intensity. The registered MCAR slope rule is retained.
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
#'   changed, beta_true, beta_observed), expected_log2, truth effects and manifest.
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
    mu <- unlist(lapply(seq_len(nrow(d)), function(i) {
        if (!counts[i]) return(numeric())
        # sample.int avoids R's scalar-numeric sample shortcut.
        pools[[i]][sample.int(length(pools[[i]]), counts[i], replace = TRUE)]
    }), use.names = FALSE)
    changed <- rep(FALSE, np)
    number_changed <- as.integer(round(changed_fraction * np))
    if (number_changed) changed[sample.int(np, number_changed)] <- TRUE
    effect <- beta * changed
    sample_effects <- matrix(stats::rnorm(np * length(sample_ids)), np) *
        sqrt(d$sigma2_s[rows])
    batch_effects <- matrix(stats::rnorm(np * length(batch_ids)), np) *
        sqrt(d$sigma2_b[rows])
    dimnames(sample_effects) <- list(protein_ids, sample_ids)
    dimnames(batch_effects) <- list(protein_ids, batch_ids)
    expected <- matrix(mu, np, nc) +
        sample_effects[, match(metadata$sample, sample_ids), drop = FALSE] +
        batch_effects[, match(metadata$batch, batch_ids), drop = FALSE] +
        outer(k * effect, as.integer(metadata$condition == "B"))
    value <- expected + matrix(stats::rnorm(np * nc), np) * sqrt(d$sigma2_e[rows])
    probability <- stats::plogis(-(curve$intercept + slope * (expected - center)))
    missing_values <- matrix(stats::runif(np * nc), np) < probability
    value[missing_values] <- NA_real_
    dimnames(value) <- dimnames(expected) <- list(protein_ids, rownames(metadata))
    protein_metadata <- data.frame(protein = protein_ids, tier = d$tier[rows],
        mu = mu, targeted = changed, changed = effect != 0, beta_true = effect, beta_observed = k * effect)
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
        baseline = "Empirical tier mean resampling; point mass when only tier mean supplied",
        missingness = "Curve evaluated at mu + sample + batch + condition effect, before cell residual",
        provenance = config$provenance,
        wall_seconds = proc.time()[["elapsed"]] - started)
    list(matrix = value, cell_metadata = metadata, protein_metadata = protein_metadata,
        expected_log2 = expected, truth = list(sample_effects = sample_effects,
        batch_effects = batch_effects), manifest = manifest)
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
