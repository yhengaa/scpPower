simulation_fixture <- function() list(
    tiers = data.frame(tier = c(1,2), proteins = c(20,20), mean_log2 = c(8,12),
        sigma2_s = c(.2,.3), sigma2_b = c(.1,.2), sigma2_e = c(.3,.4)),
    curve = list(intercept=1,effective_slope=.6,intensity_center=10),
    plex_size=8L,label_free=FALSE,k=1)

 test_that("seeded simulation is reproducible and preserves caller RNG state", {
    p <- simulation_fixture();set.seed(721);before <- .Random.seed;kind <- RNGkind()
    a <- scp_simulate(NULL,2,7,beta=1,preset=p,seed=91)
    b <- scp_simulate(NULL,2,7,beta=1,preset=p,seed=91)
    expect_identical(a$matrix,b$matrix)
    expect_identical(a$cell_metadata,b$cell_metadata)
    expect_identical(a$protein_metadata,b$protein_metadata)
    expect_identical(a$truth,b$truth)
    expect_identical(.Random.seed,before);expect_identical(RNGkind(),kind)
    expect_equal(dim(a$matrix),c(40,28))
    expect_equal(as.integer(table(a$cell_metadata$sample)),rep(7L,4))
    expect_equal(sum(a$protein_metadata$changed),4)
    expect_equal(a$manifest$seed,91L)
    expect_true(all(c('R','scp','limma','edgeR','scplainer')%in%names(a$manifest$versions)))
    expect_true(all(nchar(a$manifest$input_hashes)==64L))
 })

 test_that("missingness follows expected intensity, independently of cell residual", {
    p <- simulation_fixture();p$tiers$proteins<-c(100,100);p$tiers$mean_log2<-10
    p$tiers$sigma2_s<-p$tiers$sigma2_b<-0;p$tiers$sigma2_e<-25
    p$curve<-list(intercept=0,effective_slope=1,intensity_center=10)
    a<-scp_simulate(NULL,1,100,beta=0,preset=p,seed=93)
    expect_equal(range(a$expected_log2),c(10,10))
    expect_false(any(a$protein_metadata$changed))
    expect_equal(mean(is.na(a$matrix)),.5,tolerance=.02)
    q<-a$matrix;q[is.na(q)]<-a$expected_log2[is.na(q)]
    # Direct expected-intensity semantics are also checked through low/high
    # expected intensities while the realized residual SD remains large.
    p$tiers$mean_log2<-c(7,13)
    b<-scp_simulate(NULL,1,100,beta=0,preset=p,seed=93)
    expect_gt(mean(is.na(b$matrix[1:100,])),.9)
    expect_lt(mean(is.na(b$matrix[101:200,])),.1)
    p$curve$intercept<-4
    c<-scp_simulate(NULL,1,100,beta=0,preset=p,seed=93)
    expect_lt(mean(is.na(c$matrix)),mean(is.na(b$matrix)))
 })

 test_that("balanced and confounded layouts change batch-induced contrasts", {
    p<-simulation_fixture();p$plex_size<-20;p$tiers$proteins<-c(200,200)
    p$tiers$sigma2_s<-p$tiers$sigma2_e<-0;p$tiers$sigma2_b<-9
    p$curve<-list(intercept=40,effective_slope=0,intensity_center=10)
    balanced<-scp_simulate(NULL,3,20,beta=0,preset=p,layout='balanced',seed=94)
    confounded<-scp_simulate(NULL,3,20,beta=0,preset=p,layout='confounded',seed=94)
    contrast<-function(r)rowMeans(r$matrix[,r$cell_metadata$condition=='B'])-rowMeans(r$matrix[,r$cell_metadata$condition=='A'])
    expect_lt(max(abs(contrast(balanced))),1e-12)
    expect_gt(sd(contrast(confounded)),1)
    expect_true(all(vapply(split(balanced$cell_metadata$condition,balanced$cell_metadata$batch),function(v)sum(v=='A')==sum(v=='B'),logical(1))))
    expect_true(all(vapply(split(confounded$cell_metadata$condition,confounded$cell_metadata$batch),function(v)length(unique(v))==1L,logical(1))))
 })

 test_that("sample effects are shared across cells and k scales only the condition effect", {
    p<-simulation_fixture();p$tiers$sigma2_b<-p$tiers$sigma2_e<-0
    p$curve<-list(intercept=40,effective_slope=0,intensity_center=10)
    a<-scp_simulate(NULL,2,10,beta=2,preset=p,k=1,seed=97,changed_fraction=1)
    b<-scp_simulate(NULL,2,10,beta=2,preset=p,k=.5,seed=97,changed_fraction=1)
    expect_equal(dim(a$truth$sample_effects),c(40,4))
    expect_true(all(vapply(split(seq_len(nrow(a$cell_metadata)),a$cell_metadata$sample),function(j)all(apply(a$matrix[,j],1,sd)==0),logical(1))))
    conditionB<-a$cell_metadata$condition=='B'
    expect_equal(unname(a$matrix[,conditionB]-b$matrix[,conditionB]),matrix(1,40,20))
    expect_identical(a$matrix[,!conditionB],b$matrix[,!conditionB])
 })

 test_that("custom layouts retain native unequal designs and label-free absorption", {
    p<-simulation_fixture()
    layout<-matrix(c(2,0,1, 0,3,1),3,2,dimnames=list(c('L1','L2','L3'),c('run1','run2')))
    r<-scp_simulate(NULL,c(2,1),c(2,3,2),beta=0,preset=p,layout=layout,seed=99)
    expect_equal(dim(r$matrix),c(40,7))
    expect_equal(unclass(table(factor(r$cell_metadata$sample,levels=rownames(layout)),factor(r$cell_metadata$batch,levels=colnames(layout)))),unname(layout),ignore_attr=TRUE)
    p$label_free<-TRUE;p$tiers$sigma2_b<-0
    native<-matrix(1,1,8,dimnames=list('HeLa',paste0('run',1:8)))
    r<-scp_simulate(NULL,c(1,0),8,beta=0,preset=p,layout=native,seed=99)
    expect_equal(length(unique(r$cell_metadata$batch)),8)
    expect_true(all(r$truth$batch_effects==0))
    expect_error(scp_simulate(NULL,c(1,0),8,beta=1,preset=p,layout=native),'beta=0')
    expect_error(scp_simulate(NULL,c(1,0),8,beta=0,preset=p,layout=matrix(8,1,1)),'one cell per batch')
    expect_error(scp_simulate(NULL,1,c(2,3),0,preset=simulation_fixture()),'equal total')
 })

 test_that("YAML user presets and calibration objects warn on version mismatch", {
    p<-simulation_fixture();path<-tempfile(fileext='.yaml');yaml::write_yaml(p,path)
    a<-scp_simulate(NULL,1,3,0,preset=p,seed=101)
    b<-scp_simulate(NULL,1,3,0,preset=path,seed=101)
    expect_identical(a$matrix,b$matrix)
    p$versions<-list(R='0.0')
    expect_warning(r<-scp_simulate(NULL,1,3,0,preset=p,seed=101),'mismatch.*continuing')
    expect_equal(dim(r$matrix),c(40,6))
    set.seed(882)
    z<-matrix(rnorm(12*30),12,30)+seq_len(12)
    z[matrix(runif(length(z)),12,30)<.15]<-NA_real_
    dimnames(z)<-list(paste0('P',1:12),paste0('C',1:30))
    md<-data.frame(sample=rep(paste0('S',1:3),each=10),batch=rep(rep(paste0('B',1:5),each=2),3),row.names=colnames(z))
    cal<-scp_calibrate(z,md,sample_col='sample',batch_col='batch')
    cal$tiers$sigma2_s<-0
    r<-scp_simulate(cal,1,5,0,seed=101)
    expect_equal(dim(r$matrix),c(12,10))
    multi<-scp_simulate(list(custom=cal),1,5,0,preset='custom',seed=101)
    expect_identical(r$matrix,multi$matrix)
    p$versions<-NULL
    p$tiers$sigma2_s[1]<-NA_real_
    expect_error(scp_simulate(NULL,1,3,0,preset=p),'Unidentified')
 })

protein_simulation_fixture <- function(mode = 'protein_rate') {
    p <- simulation_fixture()
    p$protein_parameters <- data.frame(
        protein=c('lowA','lowB','highA','highB'),tier=c(1,1,2,2),
        mean_log2=c(7,9,11,13),sigma2_s=c(.01,.02,.03,.04),
        sigma2_b=c(.1,.2,.3,.4),sigma2_e=c(.5,.6,.7,.8),
        r_g=c(.2,.4,.6,.8),p_g=c(.25,.5,.75,1),q_g=c(.8,.6,.4,.2))
    p$detection_mode <- mode
    p
}

test_that('G1 jointly resamples complete calibrated protein rows within tiers', {
    p <- protein_simulation_fixture()
    a <- scp_simulate(NULL,2,20,0,preset=p,seed=107)
    b <- scp_simulate(NULL,2,20,0,preset=p,seed=107)
    idx <- match(a$protein_metadata$source_protein,p$protein_parameters$protein)
    expect_false(anyNA(idx))
    expect_equal(a$protein_metadata$tier,p$protein_parameters$tier[idx])
    expect_equal(a$protein_metadata$mu,p$protein_parameters$mean_log2[idx])
    for(field in c('sigma2_s','sigma2_b','sigma2_e','r_g'))
        expect_equal(a$protein_metadata[[field]],p$protein_parameters[[field]][idx])
    expect_true(anyDuplicated(a$protein_metadata$source_protein)>0)
    expect_identical(a$matrix,b$matrix)
    expect_identical(a$truth,b$truth)
    expect_equal(a$manifest$missingness,'protein_rate')
    expect_equal(a$manifest$mandate_addendum,'1.1.3')
    # The tier summaries are reporting values rather than generating variances.
    altered <- p;altered$tiers$sigma2_s<-altered$tiers$sigma2_b<-altered$tiers$sigma2_e<-50
    z <- scp_simulate(NULL,2,20,0,preset=altered,seed=107)
    expect_identical(a$matrix,z$matrix)
    path <- tempfile(fileext='.yaml');yaml::write_yaml(p,path,precision=17)
    y <- scp_simulate(NULL,2,20,0,preset=path,seed=107)
    expect_identical(a$matrix,y$matrix)
})

test_that('G2 TMT identification removes entire plexes and is effect-independent', {
    p <- protein_simulation_fixture('tmt_two_level')
    p$protein_parameters$q_g <- 1
    a <- scp_simulate(NULL,2,100,0,preset=p,seed=109)
    b <- scp_simulate(NULL,2,100,3,preset=p,seed=109,changed_fraction=.1)
    expect_identical(is.na(a$matrix),is.na(b$matrix))
    expect_identical(a$truth$plex_identified,b$truth$plex_identified)
    for(batch in colnames(a$truth$plex_identified)) {
        j <- a$cell_metadata$batch==batch
        expect_equal(unname(rowSums(!is.na(a$matrix[,j,drop=FALSE]))),
            as.numeric(a$truth$plex_identified[,batch])*sum(j))
    }
    expect_match(a$manifest$modelling_consequence,'p_g and q_g do not respond')
    p$protein_parameters$p_g <- 1;p$protein_parameters$q_g <- .4
    c <- scp_simulate(NULL,1,1000,0,preset=p,seed=111)
    expect_true(all(c$truth$plex_identified))
    expect_equal(mean(!is.na(c$matrix)),.4,tolerance=.01)
    p$protein_parameters$p_g <- 0
    z <- scp_simulate(NULL,1,10,0,preset=p,seed=111)
    expect_true(all(is.na(z$matrix)))
})

test_that('G3 DIA baseline detection and slope respond before residual noise', {
    p <- protein_simulation_fixture()
    p$tiers$proteins<-c(100,100)
    p$protein_parameters$r_g<-.5
    p$protein_parameters$sigma2_s<-p$protein_parameters$sigma2_b<-0
    p$protein_parameters$sigma2_e<-30
    p$curve$effective_slope<-1
    a <- scp_simulate(NULL,1,300,2,preset=p,seed=113,changed_fraction=1)
    A<-a$cell_metadata$condition=='A';B<-!A
    expect_equal(mean(!is.na(a$matrix[,A])),.5,tolerance=.01)
    expect_equal(mean(!is.na(a$matrix[,B])),plogis(2),tolerance=.01)
    q<-p;q$protein_parameters$sigma2_e<-0
    b<-scp_simulate(NULL,1,300,2,preset=q,seed=113,changed_fraction=1)
    expect_identical(is.na(a$matrix),is.na(b$matrix))
    p$protein_parameters$r_g<-c(0,1,0,1)
    z<-scp_simulate(NULL,1,10,10,preset=p,seed=115,changed_fraction=1)
    zero<-z$protein_metadata$r_g==0;one<-!zero
    expect_true(all(is.na(z$matrix[zero,,drop=FALSE])))
    expect_true(all(is.finite(z$matrix[one,,drop=FALSE])))
    expect_null(z$truth$plex_identified)
})

test_that('invalid per-protein parameter tables cannot silently fall back', {
    p<-protein_simulation_fixture();p$protein_parameters$r_g[1]<-1.01
    expect_error(scp_simulate(NULL,1,5,0,preset=p),'detection probability')
    p<-protein_simulation_fixture();p$protein_parameters$sigma2_e[1]<-NA_real_
    expect_error(scp_simulate(NULL,1,5,0,preset=p),'per-protein component')
    p<-protein_simulation_fixture();p$protein_parameters$protein[1]<-p$protein_parameters$protein[2]
    expect_error(scp_simulate(NULL,1,5,0,preset=p),'identifiers')
    p<-protein_simulation_fixture();p$protein_parameters<-NULL
    expect_error(scp_simulate(NULL,1,5,0,preset=p),'requires protein_parameters')
})

test_that('G1 distinguishes unidentified batch variance from missing protein fits', {
    p <- protein_simulation_fixture()
    cal <- list(tiers=p$tiers,missingness=p$curve,
        proteins=data.frame(protein=c('x','y','z','w'),tier=c(1,1,2,2),
            mean_log2=c(7,9,11,13),sigma2_s=c(.01,NA,.03,.04),
            sigma2_b=c(NA,NA,.3,.4),sigma2_e=c(.5,NA,.7,.8),
            missing_fraction=c(.8,.6,.4,.2),
            sigma2_b_reason=c('batch-line aliasing','Fit failed: insufficient data','estimated','estimated')),
        manifest=list(absorb_cell_batches=FALSE),compression=list(k=1))
    config<-scpPower:::.scp_simulation_config(cal,NULL)
    expect_equal(config$protein_parameters$sigma2_b[1],0)
    expect_equal(config$protein_parameters$sigma2_b[2],p$tiers$sigma2_b[1])
    expect_equal(config$protein_parameters$sigma2_s[2],p$tiers$sigma2_s[1])
    expect_equal(config$protein_parameters$sigma2_e[2],p$tiers$sigma2_e[1])
    expect_equal(config$protein_parameters$sigma2_s[1],.01)
})
