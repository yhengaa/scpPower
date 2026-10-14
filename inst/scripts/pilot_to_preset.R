# Calibration utility: converts an existing calibration without fitting a new model.
# Source this file from system.file('scripts', 'pilot_to_preset.R', package='scpPower').
pilot_to_preset <- function(calibration, endpoint, metadata, batch_col,
                            cells_per_run, minutes_per_run, carrier = FALSE,
                            instrument_source = 'Measured pilot acquisition duration') {
    build <- getFromNamespace('.scp_simulation_config', 'scpPower')
    preset <- build(calibration, NULL)
    stopifnot(is.matrix(endpoint), ncol(endpoint) == nrow(metadata),
              identical(colnames(endpoint), rownames(metadata)),
              batch_col %in% names(metadata), cells_per_run > 0, minutes_per_run > 0)
    p <- preset$protein_parameters
    index <- match(p$protein, rownames(endpoint))
    if (anyNA(index)) stop('Endpoint must contain every calibrated protein identifier.')
    mask <- is.finite(endpoint[index, , drop = FALSE])
    p$r_g <- rowMeans(mask)
    if (carrier) {
        batches <- split(seq_len(ncol(mask)), metadata[[batch_col]])
        identified <- vapply(batches, function(j) rowSums(mask[, j, drop=FALSE]) > 0,
                             logical(nrow(mask)))
        if (!is.matrix(identified)) identified <- matrix(identified,nrow=nrow(mask))
        p$p_g <- rowMeans(identified)
        den <- rowSums(vapply(seq_along(batches), function(b)
            identified[, b] * length(batches[[b]]), numeric(nrow(mask))))
        p$q_g <- ifelse(den > 0, rowSums(mask) / den, 0)
        preset$detection_mode <- 'tmt_two_level'
    } else preset$detection_mode <- 'protein_rate'
    for (tier in preset$tiers$tier) {
        rows <- p$tier == tier
        preset$tiers$missing_fraction[preset$tiers$tier == tier] <- mean(!mask[rows,,drop=FALSE])
    }
    preset$protein_parameters <- p
    preset$plex_size <- as.integer(cells_per_run)
    preset$cells_per_run <- as.integer(cells_per_run)
    preset$minutes_per_run <- minutes_per_run
    preset$minutes_per_run_source <- instrument_source
    preset$provenance$endpoint_mask_hash <- digest::digest(mask, algo='sha256')
    preset$provenance$endpoint_input_hash <- digest::digest(endpoint, algo='sha256')
    preset$provenance$pilot_conversion <- 'G1 fallback and amendment-I endpoint mask; no new fit'
    preset
}
