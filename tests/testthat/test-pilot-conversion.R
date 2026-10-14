test_that('pilot conversion uses the endpoint mask and retains mapped variances', {
    set.seed(304)
    z <- matrix(rnorm(12 * 80),12,80,
        dimnames=list(paste0('P',1:12),paste0('C',1:80))) + seq_len(12)
    meta <- data.frame(sample=rep(paste0('S',1:8),each=10),
        condition=rep(rep(c('A','B'),each=4),each=10),
        batch=rep(rep(paste0('B',1:5),each=2),8),row.names=colnames(z))
    z[sample(length(z),120)] <- NA_real_
    f <- list(x=z,meta=meta)
    bulk <- f$x[, c(1,11,21,31),drop=FALSE]
    colnames(bulk) <- paste0('D',seq_len(ncol(bulk)))
    md <- data.frame(donor=colnames(bulk),row.names=colnames(bulk))
    cal <- scp_calibrate(f$x, f$meta, sample_col='sample', condition_col='condition',
        batch_col='batch', bulk=bulk, bulk_metadata=md, carrier_ratio=100)
    utility <- new.env(parent=globalenv())
    sys.source(system.file('scripts','pilot_to_preset.R',package='scpPower'),utility)
    endpoint <- f$x
    endpoint[1,] <- NA_real_
    preset <- utility$pilot_to_preset(cal,endpoint,f$meta,'batch',16,95,carrier=TRUE)
    expect_equal(preset$protein_parameters$r_g,
        unname(rowMeans(is.finite(endpoint))[match(preset$protein_parameters$protein,rownames(endpoint))]))
    i <- match('P1',preset$protein_parameters$protein)
    expect_equal(preset$protein_parameters$p_g[i],0)
    expect_equal(preset$protein_parameters$q_g[i],0)
    expected <- cal$proteins$sigma2_s[
        match(preset$protein_parameters$protein,cal$proteins$protein)]
    unavailable <- !is.finite(expected)
    expected[unavailable] <- cal$tiers$sigma2_s[
        match(preset$protein_parameters$tier[unavailable],cal$tiers$tier)]
    expect_equal(preset$protein_parameters$sigma2_s,expected)
    expect_equal(preset$protein_parameters$sigma2_s_tier_fallback,unavailable)
    simulated <- scp_simulate(NULL,3,4,1,preset=preset,seed=913)
    expect_equal(ncol(simulated$matrix),24)
    planning <- scp_power(preset,3,4,1,replicates=2,seed=913)
    expect_equal(planning$results$instrument_hours,planning$results$runs*95/60)
})
