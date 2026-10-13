#' Estimate power for the registered single-cell proteomics endpoint
#'
#' Simulates independent biological samples and evaluates per-sample protein
#' means with limma-trend, BH 0.05 and injected-sign agreement. The registered
#' M3 grid uses all combinations of the supplied design vectors.
#'
#' @param preset Registered preset name, YAML path or complete preset list.
#' @param n_per_condition Positive integers, the sample counts per condition.
#'   Each value defines a separate design with equal numbers in A and B.
#' @param m_cells_per_sample Positive integers, cells per sample for each design.
#' @param beta Nonnegative scalar or vector of log2 effects defining the grid.
#' @param layout Balanced, confounded, or a custom sample-by-batch count matrix.
#'   Custom layouts require one n and one m.
#' @param k Compression assumption or vector of assumptions in (0,1].
#' @param replicates Independent simulations per design, default 200.
#' @param seed Nonnegative integer starting seed. Replicate seeds are seed plus
#'   100000 times the zero-based grid index plus the one-based replicate index.
#' @param target_power Target for the cheapest qualifying design, default 0.8.
#'
#' @details A protein is testable only when it has at least three quantified
#'   cells in every sample. Unweighted sample means of log2 measurements enter
#'   a design with intercept and condition; limma::eBayes uses trend=TRUE.
#'   A truly changed protein is detected only when its BH adjusted p-value is
#'   below 0.05 and the estimated fold change has the injected sign. Power
#'   includes untestable changed proteins in its denominator. Testability is
#'   the proportion of all proteins satisfying the cell minimum in every sample.
#'   Conditional power uses testable changed proteins as its denominator and
#'   is NA when that denominator is zero. It is reported beside unconditional
#'   power. The summary averages defined replicate conditional powers.
#'
#'   Realised FDR is the mean replicate false-discovery proportion among all
#'   BH rejections: truly unchanged proteins divided by all rejected proteins,
#'   with zero when nothing is rejected. Sign agreement affects power, not the
#'   null-protein definition; wrong-sign changed rejections are also reported.
#'   Mean rejections and the fraction of replicates with no rejections are
#'   reported beside FDR so sparse rejection counts remain visible.
#'   Simulation standard errors are replicate standard deviations divided by
#'   sqrt(replicates), and are NA with only one replicate. Power is NA when
#'   beta=0 because there are no truly changed proteins.
#'
#'   Runs are actual occupied batches in the simulated layout. Instrument hours
#'   are runs times the preset's cited minutes_per_run divided by 60. Balanced
#'   odd-capacity plexes retain the registered empty slot, so nominal capacity
#'   alone must not be used to undercount runs. Published method/gradient times
#'   may exclude unreported overhead; hours are nominal acquisition estimates.
#'   No batch term is added to the registered primary endpoint's design.
#'
#'   Presets remain assumptions about future designs, including bulk-to-single-
#'   cell donor transfer and TMT compression. This function does not perform M4.
#'   Version differences warn and execution continues. Each replicate records
#'   seed, matrix and metadata hashes and wall time; the manifest records shared
#'   inputs, versions, RNG kind and endpoint definitions.
#'
#' @return List with results (one row per design), cheapest (minimum instrument
#'   hours meeting target power for each beta and k, or an explicit unavailable
#'   row), replicate_metrics and manifest. No cells are imputed.
#' @examples
#' preset <- list(tiers=data.frame(tier=1, proteins=20, mean_log2=10,
#'     sigma2_s=.1, sigma2_b=0, sigma2_e=.2),
#'     curve=list(intercept=10, effective_slope=0, intensity_center=10),
#'     label_free=TRUE, plex_size=1, cells_per_run=1, minutes_per_run=35)
#' result <- scp_power(preset, 3, 4, beta=1, replicates=2, seed=23)
#' result$results
#' @export
scp_power <- function(preset, n_per_condition, m_cells_per_sample, beta,
                      layout = "balanced", k = 1, replicates = 200L,
                      seed = 20261013L, target_power = 0.8) {
    started <- proc.time()[["elapsed"]]
    config <- .scp_simulation_config(NULL, preset)
    integer_vector <- function(x, name, minimum = 1) {
        if (!is.numeric(x) || !length(x) || any(!is.finite(x)) ||
            any(x < minimum | x != floor(x) | x > .Machine$integer.max))
            stop(name, " requires valid integers.")
        as.integer(x)
    }
    n <- integer_vector(n_per_condition, "n_per_condition", 2)
    m <- integer_vector(m_cells_per_sample, "m_cells_per_sample")
    replicates <- integer_vector(replicates, "replicates")
    seed <- integer_vector(seed, "seed", 0)
    if (length(replicates) != 1L || length(seed) != 1L)
        stop("replicates and seed must each be scalar.")
    if (!is.numeric(beta) || !length(beta) || any(!is.finite(beta) | beta < 0))
        stop("beta must contain finite nonnegative log2 effects.")
    if (!is.numeric(k) || !length(k) || any(!is.finite(k) | k <= 0 | k > 1))
        stop("k must be in (0,1].")
    if (length(target_power) != 1L || !is.finite(target_power) ||
        target_power < 0 || target_power > 1) stop("target_power must be in [0,1].")
    for (field in c("cells_per_run", "minutes_per_run")) {
        x <- config[[field]]
        if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0)
            stop("Preset requires a positive cited ", field, ".")
    }
    if (is.matrix(layout) && (length(n) != 1L || length(m) != 1L))
        stop("Custom layouts require a single n and m.")
    grid <- expand.grid(n_per_condition=unique(n), m_cells_per_sample=unique(m),
        beta=unique(beta), k=unique(k), KEEP.OUT.ATTRS=FALSE)
    if (seed + 100000 * (nrow(grid) - 1) + replicates > .Machine$integer.max)
        stop("Seed range exceeds the integer limit.")
    versions <- .scp_versions()
    if (!is.null(config$versions)) {
        ref <- unlist(config$versions); common <- intersect(names(ref), names(versions))
        mismatch <- common[ref[common] != versions[common]]
        if (length(mismatch)) warning("Software version mismatch: ",
            paste(mismatch,collapse=", "), "; continuing.", call.=FALSE)
        # Warn once per scp_power run rather than once per simulation.
        config$versions <- as.list(versions)
    }
    vt <- paste(names(versions), versions, sep="=", collapse="; ")
    records <- vector("list", nrow(grid)); summaries <- vector("list", nrow(grid))
    for (i in seq_len(nrow(grid))) {
        point_started <- proc.time()[["elapsed"]]
        g <- grid[i, , drop=FALSE]; rr <- vector("list", replicates)
        for (j in seq_len(replicates)) {
            replicate_started <- proc.time()[["elapsed"]]
            simseed <- as.integer(seed + 100000 * (i - 1) + j)
            sim <- scp_simulate(NULL, g$n_per_condition, g$m_cells_per_sample,
                g$beta, preset=config, layout=layout, k=g$k, seed=simseed)
            endpoint <- .scp_primary_endpoint(sim)
            rr[[j]] <- data.frame(grid_index=i, replicate=j, seed=simseed,
                n_per_condition=g$n_per_condition, m_cells_per_sample=g$m_cells_per_sample,
                beta=g$beta, k=g$k, power=endpoint$power,
                conditional_power=endpoint$conditional_power,
                testable_changed=endpoint$testable_changed,
                testability=endpoint$testability, realised_fdr=endpoint$realised_fdr,
                changed=endpoint$changed, detected=endpoint$detected,
                testable=endpoint$testable, rejected=endpoint$rejected,
                false_discoveries=endpoint$false_discoveries,
                wrong_sign_changed=endpoint$wrong_sign_changed,
                runs=sim$manifest$batches,
                matrix_sha256=unname(sim$manifest$input_hashes["matrix"]),
                metadata_sha256=unname(sim$manifest$input_hashes["cell_metadata"]),
                wall_seconds=proc.time()[["elapsed"]]-replicate_started,
                software_versions=vt)
        }
        r <- do.call(rbind, rr); records[[i]] <- r
        avg <- function(x) if (all(is.na(x))) NA_real_ else mean(x,na.rm=TRUE)
        se <- function(x) if (sum(is.finite(x))>1L)
            stats::sd(x,na.rm=TRUE)/sqrt(sum(is.finite(x))) else NA_real_
        summaries[[i]] <- data.frame(g, power=avg(r$power),
            conditional_power=avg(r$conditional_power),
            conditional_power_se=se(r$conditional_power),
            conditional_power_defined_replicates=sum(is.finite(r$conditional_power)),
            testability=avg(r$testability), realised_fdr=avg(r$realised_fdr),
            mean_rejections=mean(r$rejected),
            zero_rejection_fraction=mean(r$rejected==0L),
            simulation_se=se(r$power), testability_se=se(r$testability),
            realised_fdr_se=se(r$realised_fdr), replicates=replicates,
            cells=2*g$n_per_condition*g$m_cells_per_sample,
            runs=r$runs[1], cells_per_run=config$cells_per_run,
            minutes_per_run=config$minutes_per_run,
            instrument_hours=r$runs[1]*config$minutes_per_run/60,
            wall_seconds=proc.time()[["elapsed"]]-point_started,
            software_versions=vt)
    }
    result <- do.call(rbind,summaries); replicate_metrics <- do.call(rbind,records)
    cheapest <- do.call(rbind,lapply(unique(result$beta),function(b)
        do.call(rbind,lapply(unique(result$k),function(kk) {
            d <- result[result$beta==b & result$k==kk &
                is.finite(result$power) & result$power>=target_power,,drop=FALSE]
            if (!nrow(d)) return(data.frame(beta=b,k=kk,target_power=target_power,
                n_per_condition=NA_integer_,m_cells_per_sample=NA_integer_,
                power=NA_real_,runs=NA_integer_,instrument_hours=NA_real_,
                status="No tested design reaches target",software_versions=vt))
            d <- d[order(d$instrument_hours,d$cells,d$n_per_condition),,drop=FALSE]
            data.frame(beta=b,k=kk,target_power=target_power,
                n_per_condition=d$n_per_condition[1],m_cells_per_sample=d$m_cells_per_sample[1],
                power=d$power[1],runs=d$runs[1],instrument_hours=d$instrument_hours[1],
                status="Target reached in registered grid",software_versions=vt)
        }))))
    manifest <- list(seed=seed, RNGkind=c("Mersenne-Twister","Inversion","Rejection"),
        versions=as.list(versions), mandate_addendum="1.1.3", dataset=config$dataset,
        input_hashes=list(preset=digest::digest(config,algo="sha256"),
            layout=digest::digest(layout,algo="sha256"),grid=digest::digest(grid,algo="sha256")),
        output_hashes=list(results=digest::digest(result,algo="sha256"),
            replicate_metrics=digest::digest(replicate_metrics,algo="sha256")),
        seed_formula="seed + 100000 * zero-based grid index + one-based replicate",
        replicates=replicates, endpoint=list(minimum_quantified_cells_per_sample=3,
            testability="minimum satisfied in every sample",design="~ condition",
            inference="limma-trend, eBayes(trend=TRUE)",BH_threshold=.05,
            power="all truly changed proteins; BH rejection plus sign agreement",
            conditional_power="detected changed / testable changed; NA when none are testable; mean of defined replicates",
            rejection_reporting="mean BH rejections and fraction of replicates with zero BH rejections",
            realised_FDR="mean FDP among BH rejections, null proteins/rejections; zero for no rejections",
            simulation_se="SD of replicate estimate / sqrt(replicates)"),
        instrument_time_source=config$minutes_per_run_source,
        instrument_time_definition="Actual occupied batches times cited method/gradient duration; excludes any unreported overhead",
        target_power=target_power, k_status="Registered assumptions, not calibrated compression",
        between_donor_transfer_not_validated=TRUE, M4_registration_unchanged=TRUE,
        wall_seconds=proc.time()[["elapsed"]]-started)
    list(results=result,cheapest=cheapest,replicate_metrics=replicate_metrics,manifest=manifest)
}

.scp_primary_endpoint <- function(sim) {
    y <- sim$matrix; md <- sim$cell_metadata; samples <- unique(md$sample)
    columns <- split(seq_len(ncol(y)),factor(md$sample,levels=samples))
    counts <- vapply(columns,function(j)rowSums(is.finite(y[,j,drop=FALSE])),numeric(nrow(y)))
    testable <- rowSums(counts >= 3L) == length(samples)
    changed <- sim$protein_metadata$changed
    nchanged <- sum(changed); ntestable <- sum(testable)
    ntestable_changed <- sum(changed & testable)
    rejected <- rep(FALSE,nrow(y)); estimated <- rep(NA_real_,nrow(y))
    if (ntestable) {
        means <- vapply(columns,function(j)rowMeans(y[testable,j,drop=FALSE],na.rm=TRUE),
            numeric(ntestable))
        if (!is.matrix(means)) means <- matrix(means,nrow=ntestable)
        condition <- factor(md$condition[match(samples,md$sample)],levels=c("A","B"))
        design <- stats::model.matrix(~ condition)
        fit <- limma::eBayes(limma::lmFit(means,design),trend=TRUE)
        padj <- stats::p.adjust(fit$p.value[,2],method="BH")
        rejected[testable] <- is.finite(padj) & padj < .05
        estimated[testable] <- fit$coefficients[,2]
    }
    correct_sign <- is.finite(estimated) &
        sign(estimated)==sign(sim$protein_metadata$beta_observed)
    detected <- sum(rejected & changed & correct_sign)
    nr <- sum(rejected); nf <- sum(rejected & !changed)
    list(power=if(nchanged)detected/nchanged else NA_real_,
        conditional_power=if(ntestable_changed)detected/ntestable_changed else NA_real_,
        testable_changed=ntestable_changed,
        testability=ntestable/nrow(y),realised_fdr=if(nr)nf/nr else 0,
        changed=nchanged,detected=detected,testable=ntestable,rejected=nr,
        false_discoveries=nf,wrong_sign_changed=sum(rejected & changed & !correct_sign))
}
